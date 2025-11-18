#!/usr/bin/env Rscript

# ------------------------------------------------------------
# Length/Weight Hard-Cutoff + Z-Standardization (2 covariates)
# ------------------------------------------------------------
# Args:
#   1: tablesFolder
#   2: project_number (unused, for symmetry/logs)
#   3: docs_folder
#   4: pregnancy_file (tablesFolder/pregnancy.gz)
#   5: variables_mapping path (pipeline_HDGB_compatible/scripts/resources/variable_mapping)
#   6: output .gz path (child_anthropometrics_standardized_hardcutoff.gz)
#   7: change_log csv path
#   8: tp_summary csv path
#   9: qc_outdir (folder for extra logs)
#
# Output:
#   - child_anthropometrics_standardized_hardcutoff.gz
#   - hardcutoff_change_log.csv
#   - hardcutoff_timepoint_summary.csv
#   - standardization_hardcutoff.md
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(gamlss)
  library(glue)
  library(janitor)
  library(cli)
})

args <- commandArgs(TRUE)
tablesFolder   <- args[1]
project_number <- args[2]
docs_folder    <- args[3]
preg_file      <- args[4]
outfile        <- args[5]
logfile        <- args[6]
summaryfile    <- args[7]
qc_outdir      <- args[8]


dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
dir.create(qc_outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(docs_folder, recursive = TRUE, showWarnings = FALSE)

# ---------- Import vars mapping (provides length_columns, weight_columns, bmi_columns, age_columns, default_columns, etc.)
source("pipeline_HDGB_compatible/scripts/utils/variables_mapping.R")

# ---------- Load base data (like lw-standardization.R)
values <- read.table(
  file = file.path(tablesFolder, "child_anthropometrics.gz"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

pregnancy_values <- read.table(
  file = preg_file,
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
) %>% select(child_id, pregnancy_duration)

values <- values %>%
  left_join(pregnancy_values, by = "child_id")

# Prepare pregnancy_duration_1 as in original pipeline
median_duration <- median(values$pregnancy_duration, na.rm = TRUE)
values$pregnancy_duration_1 <- ifelse(is.na(values$pregnancy_duration), median_duration, values$pregnancy_duration)

# ---------- Timepoints & windows
tps <- gsub("^age_", "", age_columns)
tp_ranges <- list(
  birth = c(0, 0),
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
  `8y`  = c(2739, 3300)
)
tp_prev <- setNames(c(NA, tps[-length(tps)]), tps)
tp_next <- setNames(c(tps[-1], NA), tps)

age_col <- function(tp) paste0("age_", tp)
len_col <- function(tp) paste0("length_", tp)
wgt_col <- function(tp) paste0("weight_", tp)
bmi_col <- function(tp) paste0("bmi_", tp)

in_range <- function(x, rng) !is.na(x) & x >= rng[1] & x < rng[2]
in_range_safe <- function(x, rng) if (any(is.na(rng))) rep(FALSE, length(x)) else in_range(x, rng)

# ---------- Logging helpers
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

# ---------- BEFORE snapshot (for untouched + coverage)
before_cols <- c("child_sentrix_id",
                 as.vector(outer(c("age","length","weight"), tps, paste, sep = "_")))
before <- values[, before_cols, drop = FALSE]

# ---------- 0) Negative age fix: abs()
for (tp in tps) {
  ac <- age_col(tp)
  if (!ac %in% names(values)) values[[ac]] <- NA_real_
  neg_idx <- which(!is.na(values[[ac]]) & values[[ac]] < 0)
  if (length(neg_idx)) {
    values[[ac]][neg_idx] <- abs(values[[ac]][neg_idx])
    log_add(tp, "age_negative_to_abs", length(neg_idx))
  }
}

# ---------- Span computation
compute_spans <- function(DT) {
  spans <- vector("list", length(tps)); names(spans) <- tps
  for (tp in tps) {
    if (tp == "birth") {
      spans[[tp]] <- list(span_down = rep(FALSE, nrow(DT)),
                          span_up = rep(FALSE, nrow(DT)))
      next
    }
    ac <- age_col(tp)
    rng <- tp_ranges[[tp]]
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

# ---------- Step 1) Rescue only if neighbor totally missing BOTH traits
spans <- compute_spans(values)

for (tp in tps) {
  if (tp == "birth") next
  ac <- age_col(tp); lc <- len_col(tp); wc <- wgt_col(tp)
  
  # span UP -> nexttp
  nexttp <- tp_next[[tp]]
  if (!is.na(nexttp)) {
    nac <- age_col(nexttp); nlc <- len_col(nexttp); nwc <- wgt_col(nexttp)
    idx <- which(spans[[tp]]$span_up &
                   is.na(values[[nlc]]) & is.na(values[[nwc]]))  # BOTH missing
    if (length(idx)) {
      values[[nac]][idx] <- values[[ac]][idx]
      values[[nlc]][idx] <- values[[lc]][idx]
      values[[nwc]][idx] <- values[[wc]][idx]
      values[[ac]][idx] <- NA_real_
      values[[lc]][idx] <- NA_real_
      values[[wc]][idx] <- NA_real_
      log_add(tp, sprintf("reassign_span_up_to_%s_both_traits_missing", nexttp), length(idx))
    }
  }
  # span DOWN -> prevtp
  prevtp <- tp_prev[[tp]]
  if (!is.na(prevtp)) {
    pac <- age_col(prevtp); plc <- len_col(prevtp); pwc <- wgt_col(prevtp)
    idx <- which(spans[[tp]]$span_down &
                   is.na(values[[plc]]) & is.na(values[[pwc]]))  # BOTH missing
    if (length(idx)) {
      values[[pac]][idx] <- values[[ac]][idx]
      values[[plc]][idx] <- values[[lc]][idx]
      values[[pwc]][idx] <- values[[wc]][idx]
      values[[ac]][idx] <- NA_real_
      values[[lc]][idx] <- NA_real_
      values[[wc]][idx] <- NA_real_
      log_add(tp, sprintf("reassign_span_down_to_%s_both_traits_missing", prevtp), length(idx))
    }
  }
}
spans <- compute_spans(values)

# ---------- Step 2) Pool close spans (<20d) only if neighbor AGE inlier; both traits on both sides
pool_threshold <- 20
for (tp in tps) {
  if (tp == "birth") next
  ac <- age_col(tp); lc <- len_col(tp); wc <- wgt_col(tp)
  
  # span UP -> next
  nexttp <- tp_next[[tp]]
  if (!is.na(nexttp)) {
    nac <- age_col(nexttp); nlc <- len_col(nexttp); nwc <- wgt_col(nexttp)
    next_rng <- tp_ranges[[nexttp]]
    idx <- which(spans[[tp]]$span_up &
                   !is.na(values[[lc]]) & !is.na(values[[wc]]) &
                   !is.na(values[[nlc]]) & !is.na(values[[nwc]]) &
                   in_range_safe(values[[nac]], next_rng))
    if (length(idx)) {
      both_age <- !is.na(values[[ac]][idx]) & !is.na(values[[nac]][idx])
      abs_diff <- ifelse(both_age, abs(values[[ac]][idx] - values[[nac]][idx]), NA_real_)
      pool_idx <- idx[which(both_age & abs_diff < pool_threshold)]
      if (length(pool_idx)) {
        pooled_len <- (values[[lc]][pool_idx] + values[[nlc]][pool_idx]) / 2
        pooled_wgt <- (values[[wc]][pool_idx] + values[[nwc]][pool_idx]) / 2
        values[[nlc]][pool_idx] <- pooled_len
        values[[nwc]][pool_idx] <- pooled_wgt
        values[[nac]][pool_idx] <- values[[ac]][pool_idx] # neighbor age <- focal age
        values[[lc]][pool_idx] <- NA_real_
        values[[wc]][pool_idx] <- NA_real_
        values[[ac]][pool_idx] <- NA_real_
        log_add(tp, sprintf("pooled_span_up_to_%s_lt20days_neighbor_inlier", nexttp), length(pool_idx))
      }
    }
  }
  
  # span DOWN -> prev
  prevtp <- tp_prev[[tp]]
  if (!is.na(prevtp)) {
    pac <- age_col(prevtp); plc <- len_col(prevtp); pwc <- wgt_col(prevtp)
    prev_rng <- tp_ranges[[prevtp]]
    idx <- which(spans[[tp]]$span_down &
                   !is.na(values[[lc]]) & !is.na(values[[wc]]) &
                   !is.na(values[[plc]]) & !is.na(values[[pwc]]) &
                   in_range_safe(values[[pac]], prev_rng))
    if (length(idx)) {
      both_age <- !is.na(values[[ac]][idx]) & !is.na(values[[pac]][idx])
      abs_diff <- ifelse(both_age, abs(values[[ac]][idx] - values[[pac]][idx]), NA_real_)
      pool_idx <- idx[which(both_age & abs_diff < pool_threshold)]
      if (length(pool_idx)) {
        pooled_len <- (values[[lc]][pool_idx] + values[[plc]][pool_idx]) / 2
        pooled_wgt <- (values[[wc]][pool_idx] + values[[pwc]][pool_idx]) / 2
        values[[plc]][pool_idx] <- pooled_len
        values[[pwc]][pool_idx] <- pooled_wgt
        values[[pac]][pool_idx] <- values[[ac]][pool_idx]
        values[[lc]][pool_idx] <- NA_real_
        values[[wc]][pool_idx] <- NA_real_
        values[[ac]][pool_idx] <- NA_real_
        log_add(tp, sprintf("pooled_span_down_to_%s_lt20days_neighbor_inlier", prevtp), length(pool_idx))
      }
    }
  }
}
spans <- compute_spans(values)

# ---------- Step 3) Impute age_* true outliers & NAs to median in-range (birth=0)
for (tp in tps) {
  ac <- age_col(tp)
  a  <- values[[ac]]
  if (tp == "birth") {
    idx <- which(is.na(a) | a != 0)
    if (length(idx)) {
      values[[ac]][idx] <- 0
      log_add(tp, "impute_birth_age_to_0", length(idx))
    }
    next
  }
  rng     <- tp_ranges[[tp]]
  prevtp1 <- tp_prev[[tp]]
  nexttp1 <- tp_next[[tp]]
  prev_rng <- if (!is.na(prevtp1)) tp_ranges[[prevtp1]] else c(NA_real_, NA_real_)
  next_rng <- if (!is.na(nexttp1)) tp_ranges[[nexttp1]] else c(NA_real_, NA_real_)
  
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
    values[[ac]][true_outlier_idx] <- med_val
    log_add(tp, "impute_age_true_outliers_to_median", length(true_outlier_idx))
  }
  if (length(na_idx)) {
    values[[ac]][na_idx] <- med_val
    log_add(tp, "impute_age_NA_to_median", length(na_idx))
  }
}

# ---------- BEFORE/AFTER summary (untouched + coverage)
after_cols <- before_cols
after <- values[, after_cols, drop = FALSE]
tp_summary <- lapply(tps, function(tp) {
  ac <- age_col(tp); lc <- len_col(tp); wc <- wgt_col(tp)
  a0 <- before[[ac]]; a1 <- after[[ac]]
  l0 <- before[[lc]]; l1 <- after[[lc]]
  w0 <- before[[wc]]; w1 <- after[[wc]]
  age_changed    <- !num_equal(a0, a1)
  length_changed <- !num_equal(l0, l1)
  weight_changed <- !num_equal(w0, w1)
  untouched <- !(age_changed | length_changed | weight_changed)
  tibble::tibble(
    timepoint = tp,
    n_total = nrow(values),
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
}) %>% bind_rows()
readr::write_csv(tp_summary, summaryfile)

# ---------- Compute BMI from corrected L/W
for (i in seq_along(tps)) {
  tp <- tps[i]
  L <- len_col(tp); W <- wgt_col(tp); B <- bmi_col(tp)
  if (!L %in% names(values)) values[[L]] <- NA_real_
  if (!W %in% names(values)) values[[W]] <- NA_real_
  values[[B]] <- ifelse(!is.na(values[[L]]) & !is.na(values[[W]]),
                        10000 * values[[W]] / (values[[L]] * values[[L]]), NA_real_)
}

# ---------- Standardization (2 covariates: age_<tp> + pregnancy_duration_1)

# Helper: compute z for NO and LOGNO per observation given mu, sigma and y
z_from_mu_sigma <- function(y, mu, sigma, family_name) {
  if (family_name == "NO") {
    return( (y - mu) / sigma )
  } else if (family_name == "LOGNO") {
    return( (log(y) - mu) / sigma )
  } else {
    stop("Unsupported family for z_from_mu_sigma.")
  }
}

# Fit per sex, per timepoint (y ~ fp(age_tp) + fp(pregnancy_duration_1)), sigma ~ pregnancy_duration_1
standardize2D_by_sex <- function(values, y, zY, tp, family_name) {
  # family object
  fam <- if (family_name == "NO") NO else if (family_name == "LOGNO") LOGNO else stop("Unsupported family")
  ac <- age_col(tp)
  # Train data per sex (unrelated only)
  ref <- values %>% filter(unrelated_children == 1,
                           !is.na(sex), sex %in% c(1,2),
                           !is.na(.data[[y]]),
                           !is.na(.data[[ac]]),
                           !is.na(pregnancy_duration_1),
                           is.finite(.data[[y]]))
  if (nrow(ref) < 10) return(values) # too few to fit
  # prediction data
  pred <- values %>% filter(!is.na(sex), sex %in% c(1,2),
                            !is.na(.data[[y]]),
                            !is.na(.data[[ac]]),
                            !is.na(pregnancy_duration_1),
                            is.finite(.data[[y]]))
  if (nrow(pred) == 0) {
    values[[zY]] <- NA_real_
    return(values)
  }
  
  # model formulas
  # you can simplify fp() to linear if needed
  f_mu    <- as.formula(glue("{y} ~ fp({ac}) + fp(pregnancy_duration_1)"))
  f_sigma <- ~ pregnancy_duration_1
  
  # split by sex
  run_sex <- function(sex_code) {
    ref_sex  <- ref  %>% filter(sex == sex_code)
    pred_sex <- pred %>% filter(sex == sex_code)
    if (length(unique(ref_sex[[y]])) < 2) return(NULL)
    
    model <- tryCatch(
      gamlss(formula = f_mu, sigma.formula = f_sigma, family = fam, data = ref_sex),
      error = function(e) NULL
    )
    if (is.null(model)) {
      # fallback linear
      f_mu2 <- as.formula(glue("{y} ~ {ac} + pregnancy_duration_1"))
      model <- tryCatch(
        gamlss(formula = f_mu2, sigma.formula = f_sigma, family = fam, data = ref_sex),
        error = function(e) NULL
      )
      if (is.null(model)) return(NULL)
    }
    
    # predict mu, sigma at each observation's covariates
    mu_hat <- predict(model, what = "mu",    newdata = pred_sex, type = "response")
    sg_hat <- predict(model, what = "sigma", newdata = pred_sex, type = "response")
    
    z <- tryCatch(
      z_from_mu_sigma(pred_sex[[y]], mu_hat, sg_hat, family_name),
      error = function(e) rep(NA_real_, nrow(pred_sex))
    )
    tibble(child_id = pred_sex$child_id, z = z)
  }
  
  mz <- run_sex(1)
  fz <- run_sex(2)
  zdf <- bind_rows(mz, fz)
  
  values <- values %>%
    left_join(zdf %>% rename(!!zY := z), by = "child_id")
  
  # wipe extreme |z|>5 as original
  values[[zY]][!is.na(values[[zY]]) & abs(values[[zY]]) > 5] <- NA_real_
  
  values
}

# Run standardization for Length, Weight, BMI
for (phenoI in seq_along(length_columns)) {
  y  <- length_columns[phenoI]
  tp <- strsplit(y, "_")[[1]][2]
  zY <- paste0("z_", y)
  cli::cli_inform(glue("Standardizing {y} with covariates age_{tp} + pregnancy_duration_1"))
  values <- standardize2D_by_sex(values, y, zY, tp, "NO")
}
for (phenoI in seq_along(weight_columns)) {
  y  <- weight_columns[phenoI]
  tp <- strsplit(y, "_")[[1]][2]
  zY <- paste0("z_", y)
  cli::cli_inform(glue("Standardizing {y} with covariates age_{tp} + pregnancy_duration_1"))
  values <- standardize2D_by_sex(values, y, zY, tp, "NO")
}
# Ensure bmi_columns list exists for all tps; recreate from tps if not
if (!exists("bmi_columns")) {
  bmi_columns <- paste0("bmi_", tps)
}
for (phenoI in seq_along(bmi_columns)) {
  y  <- bmi_columns[phenoI]
  tp <- strsplit(y, "_")[[1]][2]
  zY <- paste0("z_", y)
  cli::cli_inform(glue("Standardizing {y} (LOGNO) with covariates age_{tp} + pregnancy_duration_1"))
  values <- standardize2D_by_sex(values, y, zY, tp, "LOGNO")
}

# ---------- Rank-based inverse normal transform (per trait & timepoint)
rinvnorm <- function(x) {
  idx <- which(!is.na(x))
  r   <- rank(x[idx], ties.method = "average")
  n   <- length(idx)
  z   <- rep(NA_real_, length(x))
  z[idx] <- qnorm((r - 0.5) / n)
  z
}
for (tp in tps) {
  L <- len_col(tp); W <- wgt_col(tp); B <- bmi_col(tp)
  if (L %in% names(values)) values[[paste0("length_ir_", tp)]] <- rinvnorm(values[[L]])
  if (W %in% names(values)) values[[paste0("weight_ir_", tp)]] <- rinvnorm(values[[W]])
  if (B %in% names(values)) values[[paste0("bmi_ir_", tp)]]    <- rinvnorm(values[[B]])
}

# ---------- Save logs
log_df <- bind_rows(log_rows)
if (nrow(log_df) == 0) {
  log_df <- tibble(timepoint = character(), action = character(), n_changed = integer())
}
totals <- log_df %>%
  group_by(action) %>%
  summarise(n_changed = sum(n_changed), .groups = "drop") %>%
  mutate(timepoint = "TOTAL") %>%
  select(timepoint, action, n_changed)
log_final <- bind_rows(log_df, totals) %>% arrange(timepoint, action)
write_csv(log_final, logfile)

# ---------- Docs stub
docs_file <- file.path(docs_folder, "standardization_hardcutoff.md")
writeLines(c(
  "# Standardization with hard-cutoff",
  "",
  "- Ages corrected by hard-cutoff (rescue, pooled <20d, median impute for true outliers/NA).",
  "- BMI recomputed from corrected length/weight.",
  "- Z-scores via GAMLSS (sex-specific) with covariates: age_<tp> + pregnancy_duration_1.",
  "- Added rank-based inverse normal (IR) columns for length/weight/BMI."
), con = docs_file)

# ---------- Export final table
# Keep default_columns + pregnancy_duration + all ages + L/W + BMI + z_* and *_ir_*
to_keep <- c(
  default_columns,
  "pregnancy_duration", "pregnancy_duration_1",
  age_columns,
  length_columns, weight_columns, bmi_columns,
  paste0("z_", length_columns), paste0("z_", weight_columns), paste0("z_", bmi_columns),
  paste0("length_ir_", tps), paste0("weight_ir_", tps), paste0("bmi_ir_", tps)
)
to_keep <- unique(to_keep[to_keep %in% names(values)])

values_to_export <- values %>% select(all_of(to_keep))

data.table::fwrite(values_to_export, outfile, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")

# Also write the BEFORE/AFTER summary produced earlier
readr::write_csv(tp_summary, summaryfile)

cli::cli_inform(glue("Wrote: {outfile}"))
cli::cli_inform(glue("Logs:  {logfile}, {summaryfile}"))
