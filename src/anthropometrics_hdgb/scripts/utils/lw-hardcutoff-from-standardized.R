# ------------------------------------------------------------
# MOBA anthropometrics – HARD CUTOFF (intermediate writer)
#
# Args:
#   1) tablesFolder (where child_anthropometrics_standardized.gz lives)
#   2) out_values_file (gz) -- INTERMEDIATE
#   3) out_change_log (csv)
#   4) out_timepoint_summary (csv)
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(readr)
  library(purrr)
  library(tibble)
})

args <- commandArgs(TRUE)
if (length(args) != 4) {
  stop("Usage: lw-hardcutoff-from-standardized.R <tablesFolder> <out_values_file.gz> <out_change_log.csv> <out_timepoint_summary.csv>")
}
tablesFolder          <- args[[1]]
out_values_file       <- args[[2]]
out_change_log        <- args[[3]]
out_timepoint_summary <- args[[4]]

infile <- file.path(tablesFolder, "child_anthropometrics_standardized.gz")

# --------- Load data ----------
dt <- data.table::fread(infile, na.strings = c("", "NA"))
stopifnot("child_sentrix_id" %in% names(dt))
setDT(dt)

# --------- Timepoints/ranges (days, [lower, upper)) ----------
tps <- c("birth","6w","3m","6m","8m","1y","16m","2y","3y","5y","7y","8y")
tp_ranges <- list(
  birth = c(0, 0),             # special handling
  `6w`  = c(21, 66),
  `3m`  = c(67, 136),
  `6m`  = c(137, 212),
  `8m`  = c(213, 303),
  `1y`  = c(304, 425),
  `16m` = c(426, 608),
  `2y`  = c(609, 912),
  `3y`  = c(913, 1460),
  `5y`  = c(1461, 2191),
  `7y`  = c(2192, 2738),
  `8y`  = c(2739, 3500)
)
tp_prev <- setNames(c(NA, tps[-length(tps)]), tps)
tp_next <- setNames(c(tps[-1], NA), tps)

age_col <- function(tp) sprintf("age_%s", tp)
len_col <- function(tp) sprintf("length_%s", tp)
wgt_col <- function(tp) sprintf("weight_%s", tp)

in_range      <- function(x, rng) !is.na(x) & x >= rng[1] & x < rng[2]
in_range_safe <- function(x, rng) if (any(is.na(rng))) rep(FALSE, length(x)) else in_range(x, rng)

# Ensure age columns exist (safety)
for (tp in tps) if (!age_col(tp) %in% names(dt)) dt[[age_col(tp)]] <- NA_real_

# --------- LOG scaffold ----------
log_rows <- list()
log_add <- function(tp, action, n) {
  log_rows[[length(log_rows) + 1]] <<- data.frame(
    timepoint = tp, action = action, n_changed = as.integer(n)
  )
}

num_equal <- function(a, b, tol = 1e-8) {
  both_na <- is.na(a) & is.na(b)
  both_num <- !is.na(a) & !is.na(b)
  eq_num <- rep(FALSE, length(a))
  eq_num[both_num] <- abs(a[both_num] - b[both_num]) <= tol
  both_na | eq_num
}

# Keep BEFORE snapshot (for summary)
before <- dt[, c("child_sentrix_id",
                 as.vector(outer(c("age","length","weight"), tps, paste, sep = "_"))),
             with = FALSE]

# Work copy
D <- copy(dt)

# ============================================================
# 0) MAKE AGES NON-NEGATIVE
# ============================================================
for (tp in tps) {
  ac <- age_col(tp)
  if (!ac %in% names(D)) next
  neg_idx <- which(!is.na(D[[ac]]) & D[[ac]] < 0)
  if (length(neg_idx)) {
    D[[ac]][neg_idx] <- abs(D[[ac]][neg_idx])
    log_add(tp, "age_negative_to_abs", length(neg_idx))
  }
}

compute_spans <- function(DT) {
  spans <- vector("list", length(tps)); names(spans) <- tps
  for (tp in tps) {
    if (tp == "birth") {
      spans[[tp]] <- list(span_down = rep(FALSE, nrow(DT)),
                          span_up = rep(FALSE, nrow(DT)))
      next
    }
    ac <- age_col(tp)
    rng     <- tp_ranges[[tp]]
    prevtp  <- tp_prev[[tp]]
    nexttp  <- tp_next[[tp]]
    prev_rng <- if (!is.na(prevtp)) tp_ranges[[prevtp]] else c(NA_real_, NA_real_)
    next_rng <- if (!is.na(nexttp)) tp_ranges[[nexttp]] else c(NA_real_, NA_real_)
    a <- DT[[ac]]
    inlier  <- in_range_safe(a, rng)
    outlier <- !is.na(a) & !inlier
    span_down <- outlier & in_range_safe(a, prev_rng)
    span_up   <- outlier & in_range_safe(a, next_rng)
    if (tp == "6w") span_down <- rep(FALSE, length(span_down))
    if (tp == "8y") span_up   <- rep(FALSE, length(span_up))
    spans[[tp]] <- list(span_down = span_down, span_up = span_up)
  }
  spans
}

# ============================================================
# STEP 1) RESCUE: move span → neighbor if neighbor has BOTH traits missing
# ============================================================
spans <- compute_spans(D)

for (tp in tps) {
  if (tp == "birth") next
  ac <- age_col(tp); lc <- len_col(tp); wc <- wgt_col(tp)

  # span UP → next neighbor
  nexttp <- tp_next[[tp]]
  if (!is.na(nexttp)) {
    nac <- age_col(nexttp); nlc <- len_col(nexttp); nwc <- wgt_col(nexttp)
    idx <- which(spans[[tp]]$span_up &
                   is.na(D[[nlc]]) & is.na(D[[nwc]]))
    if (length(idx)) {
      D[[nac]][idx] <- D[[ac]][idx]
      D[[nlc]][idx] <- D[[lc]][idx]
      D[[nwc]][idx] <- D[[wc]][idx]
      D[[ac]][idx] <- NA_real_; D[[lc]][idx] <- NA_real_; D[[wc]][idx] <- NA_real_
      log_add(tp, sprintf("reassign_span_up_to_%s_both_traits_missing", nexttp), length(idx))
    }
  }

  # span DOWN → previous neighbor
  prevtp <- tp_prev[[tp]]
  if (!is.na(prevtp)) {
    pac <- age_col(prevtp); plc <- len_col(prevtp); pwc <- wgt_col(prevtp)
    idx <- which(spans[[tp]]$span_down &
                   is.na(D[[plc]]) & is.na(D[[pwc]]))
    if (length(idx)) {
      D[[pac]][idx] <- D[[ac]][idx]
      D[[plc]][idx] <- D[[lc]][idx]
      D[[pwc]][idx] <- D[[wc]][idx]
      D[[ac]][idx] <- NA_real_; D[[lc]][idx] <- NA_real_; D[[wc]][idx] <- NA_real_
      log_add(tp, sprintf("reassign_span_down_to_%s_both_traits_missing", prevtp), length(idx))
    }
  }
}
spans <- compute_spans(D)

# ============================================================
# STEP 2) POOL CLOSE SPANS (<20d) with strict conditions
# ============================================================
pool_threshold <- 20
for (tp in tps) {
  if (tp == "birth") next
  ac <- age_col(tp); lc <- len_col(tp); wc <- wgt_col(tp)

  # span UP → pool into next if neighbor age inlier and both sides have both traits
  nexttp <- tp_next[[tp]]
  if (!is.na(nexttp)) {
    nac <- age_col(nexttp); nlc <- len_col(nexttp); nwc <- wgt_col(nexttp)
    next_rng <- tp_ranges[[nexttp]]
    idx <- which(spans[[tp]]$span_up &
                   !is.na(D[[lc]]) & !is.na(D[[wc]]) &
                   !is.na(D[[nlc]]) & !is.na(D[[nwc]]) &
                   in_range_safe(D[[nac]], next_rng))
    if (length(idx)) {
      has_both <- !is.na(D[[ac]][idx]) & !is.na(D[[nac]][idx])
      abs_diff <- ifelse(has_both, abs(D[[ac]][idx] - D[[nac]][idx]), NA_real_)
      pool_idx <- idx[which(has_both & abs_diff < pool_threshold)]
      if (length(pool_idx)) {
        pooled_len <- (D[[lc]][pool_idx] + D[[nlc]][pool_idx]) / 2
        pooled_wgt <- (D[[wc]][pool_idx] + D[[nwc]][pool_idx]) / 2
        D[[nlc]][pool_idx] <- pooled_len
        D[[nwc]][pool_idx] <- pooled_wgt
        D[[nac]][pool_idx] <- D[[ac]][pool_idx]
        D[[lc]][pool_idx] <- NA_real_; D[[wc]][pool_idx] <- NA_real_; D[[ac]][pool_idx] <- NA_real_
        log_add(tp, sprintf("pooled_span_up_to_%s_lt20days_neighbor_inlier", nexttp), length(pool_idx))
      }
    }
  }

  # span DOWN → pool into previous if neighbor age inlier and both sides have both traits
  prevtp <- tp_prev[[tp]]
  if (!is.na(prevtp)) {
    pac <- age_col(prevtp); plc <- len_col(prevtp); pwc <- wgt_col(prevtp)
    prev_rng <- tp_ranges[[prevtp]]
    idx <- which(spans[[tp]]$span_down &
                   !is.na(D[[lc]]) & !is.na(D[[wc]]) &
                   !is.na(D[[plc]]) & !is.na(D[[pwc]]) &
                   in_range_safe(D[[pac]], prev_rng))
    if (length(idx)) {
      has_both <- !is.na(D[[ac]][idx]) & !is.na(D[[pac]][idx])
      abs_diff <- ifelse(has_both, abs(D[[ac]][idx] - D[[pac]][idx]), NA_real_)
      pool_idx <- idx[which(has_both & abs_diff < pool_threshold)]
      if (length(pool_idx)) {
        pooled_len <- (D[[lc]][pool_idx] + D[[plc]][pool_idx]) / 2
        pooled_wgt <- (D[[wc]][pool_idx] + D[[pwc]][pool_idx]) / 2
        D[[plc]][pool_idx] <- pooled_len
        D[[pwc]][pool_idx] <- pooled_wgt
        D[[pac]][pool_idx] <- D[[ac]][pool_idx]
        D[[lc]][pool_idx] <- NA_real_; D[[wc]][pool_idx] <- NA_real_; D[[ac]][pool_idx] <- NA_real_
        log_add(tp, sprintf("pooled_span_down_to_%s_lt20days_neighbor_inlier", prevtp), length(pool_idx))
      }
    }
  }
}
spans <- compute_spans(D)

# ============================================================
# STEP 3) IMPUTE AGE_* TRUE OUTLIERS & NAs TO MEDIAN IN-RANGE
# ============================================================
for (tp in tps) {
  ac <- age_col(tp)
  a  <- D[[ac]]

  if (tp == "birth") {
    idx <- which(is.na(a) | a != 0)
    if (length(idx)) {
      D[[ac]][idx] <- 0
      log_add(tp, "impute_birth_age_to_0", length(idx))
    }
    next
  }

  rng     <- tp_ranges[[tp]]
  prevtp  <- tp_prev[[tp]]
  nexttp  <- tp_next[[tp]]
  prev_rng <- if (!is.na(prevtp)) tp_ranges[[prevtp]] else c(NA_real_, NA_real_)
  next_rng <- if (!is.na(nexttp)) tp_ranges[[nexttp]] else c(NA_real_, NA_real_)

  inlier_idx <- which(in_range_safe(a, rng))
  med_val <- if (length(inlier_idx)) median(a[inlier_idx], na.rm = TRUE) else mean(rng, na.rm = TRUE)
  if (is.na(med_val)) med_val <- mean(rng, na.rm = TRUE)

  inlier  <- in_range_safe(a, rng)
  outlier <- !is.na(a) & !inlier
  span_down_here <- outlier & in_range_safe(a, prev_rng)
  span_up_here   <- outlier & in_range_safe(a, next_rng)

  true_outlier_idx <- which(outlier & !span_down_here & !span_up_here)
  na_idx <- which(is.na(a))

  if (length(true_outlier_idx)) {
    D[[ac]][true_outlier_idx] <- med_val
    log_add(tp, "impute_age_true_outliers_to_median", length(true_outlier_idx))
  }
  if (length(na_idx)) {
    D[[ac]][na_idx] <- med_val
    log_add(tp, "impute_age_NA_to_median", length(na_idx))
  }
}

# ============================================================
# BEFORE/AFTER SUMMARY (untouched + coverage)
# ============================================================
after <- D[, c("child_sentrix_id",
               as.vector(outer(c("age","length","weight"), tps, paste, sep = "_"))),
           with = FALSE]

summ_rows <- lapply(tps, function(tp) {
  ac <- age_col(tp); lc <- len_col(tp); wc <- wgt_col(tp)
  a0 <- before[[ac]]; a1 <- after[[ac]]
  l0 <- before[[lc]]; l1 <- after[[lc]]
  w0 <- before[[wc]]; w1 <- after[[wc]]

  age_changed    <- !num_equal(a0, a1)
  length_changed <- !num_equal(l0, l1)
  weight_changed <- !num_equal(w0, w1)
  untouched <- !(age_changed | length_changed | weight_changed)

  tibble(
    timepoint = tp,
    n_total = nrow(D),
    n_untouched = sum(untouched, na.rm = TRUE),
    n_changed_any = n_total - n_untouched,
    n_age_changed = sum(age_changed, na.rm = TRUE),
    n_length_changed = sum(length_changed, na.rm = TRUE),
    n_weight_changed = sum(weight_changed, na.rm = TRUE),
    n_length_before = sum(!is.na(l0)),
    n_length_after  = sum(!is.na(l1)),
    n_weight_before = sum(!is.na(w0)),
    n_weight_after  = sum(!is.na(w1))
  )
})
timepoint_summary <- dplyr::bind_rows(summ_rows)
readr::write_csv(timepoint_summary, out_timepoint_summary)

# ============================================================
# WRITE INTERMEDIATE (FILTERED COLS) + RECOMPUTE BMI + ACTION LOG
#   - Drop any pre-existing bmi / z_* columns (input may contain old values)
#   - Recompute bmi_* from current length_* (cm) and weight_* (kg)
#   - Keep IDs/meta + all age_*, length_*, weight_*, bmi_* (drop z_* only)
# ============================================================

# 1) Drop any existing BMI and z-score columns from the working table
drop_mask <- grepl("^bmi($|_)", names(D)) |                 # matches 'bmi' and 'bmi_*'
             grepl("^z_length_", names(D)) |
             grepl("^z_weight_", names(D)) |
             grepl("^z_bmi_",    names(D))
if (any(drop_mask)) {
  D <- D[, !drop_mask, with = FALSE]
}

# 2) Recompute BMI columns (kg/m^2) from length_* (cm) and weight_* (kg)
bmi_col <- function(tp) sprintf("bmi_%s", tp)

for (tp in tps) {
  lc <- len_col(tp)
  wc <- wgt_col(tp)
  bc <- bmi_col(tp)

  # initialize target column
  D[[bc]] <- NA_real_

  # compute where both measures exist and length > 0
  valid <- !is.na(D[[wc]]) & !is.na(D[[lc]]) & (D[[lc]] > 0)
  D[[bc]][valid] <- D[[wc]][valid] / ((D[[lc]][valid] / 100)^2)
}

# 3) Select columns for output (keep bmi_*; drop all z_* just in case)
keep_cols <- names(D)[
  !grepl("^z_length_", names(D)) &
  !grepl("^z_weight_", names(D)) &
  !grepl("^z_bmi_",    names(D))
]
D_out <- D[, ..keep_cols]

# make sure folder exists
dir.create(dirname(out_values_file), recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(D_out, out_values_file, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")

log_df <- dplyr::bind_rows(log_rows)
if (nrow(log_df) == 0) {
  log_df <- tibble(timepoint = character(), action = character(), n_changed = integer())
}
totals <- log_df %>%
  group_by(action) %>%
  summarise(n_changed = sum(n_changed), .groups = "drop") %>%
  mutate(timepoint = "TOTAL") %>%
  select(timepoint, action, n_changed)

log_final <- bind_rows(log_df, totals) %>%
  arrange(timepoint, action)

readr::write_csv(log_final, out_change_log)

message("Wrote hardcutoff INTERMEDIATE: ", out_values_file)
message("Wrote action log:              ", out_change_log)
message("Wrote timepoint summary:       ", out_timepoint_summary)