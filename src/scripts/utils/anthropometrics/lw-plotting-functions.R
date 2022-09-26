##
#
# Functions used in the cleaning pipeline to plot length and weight values. 
#
##


#' Plots the weight against the length at each age for all kids before and after cleaning.
#' 
#' @param originalValues the original values
#' @param values the processed values
#' @param qcFolderLocal the folder where to save the results
#' 
#' @return returns the updated values
exportLenghtWeight <- function(
  originalValues,
  values,
  qcFolderLocal
) {
  
  timePoints <- c("Birth", "6 w", "3 m", "6 m", "8 m", "1 y", "1.5 y", "2 y", "3 y", "5 y", "7 y", "8 y", "14 y")
  
  for (ageI in 1:length(length_columns)) {
    
    lengthColumn <- length_columns[ageI]
    weightColumn <- weight_columns[ageI]
    timePoint <- timePoints[ageI]
    
    plotDF <- data.frame(
      length = c(originalValues[[lengthColumn]], values[[lengthColumn]]),
      weight = c(originalValues[[weightColumn]], values[[weightColumn]]),
      sex = c(originalValues$sex, values$sex),
      unrelated = c(originalValues$unrelated, values$unrelated),
      cleaningStep = c(rep("Before", nrow(originalValues)), rep("After", nrow(values))),
      stringsAsFactors = F
    ) %>% 
      filter(
        !is.na(length) | !is.na(weight)
    ) %>%
      mutate(
        length = ifelse(is.na(length), 0, length),
        weight = ifelse(is.na(weight), 0, weight),
        cleaningStep = factor(cleaningStep, levels = c("Before", "After")),
        sex = factor(sex, levels = c("Girl", "Boy")),
        unrelated = factor(unrelated, levels = c(0, 1))
      ) %>%
      arrange(
        unrelated
      )
    
    lwPlot <- ggplot(data = plotDF) + theme_bw() +
      geom_point(mapping = aes(x = length, y = weight, col = unrelated), alpha = 0.3) +
      facet_grid(sex ~ cleaningStep) +
      scale_color_manual(values = c("black", "blue3"), name = NULL, labels = c("Other", "Core")) +
      xlab("length [cm]") + ylab("weight [kg]") +
      ggtitle(timePoint) +
      theme(
        legend.position = "top"
      )
    
    plotFile <- file.path(qcFolderLocal, paste0(ageI, ".png"))
    png(plotFile, width = 1200, height = 1200)
    grid.draw(lwPlot)
    dummy <- dev.off()
    
  }
}


#' Exports the curves for all kids
#' 
#' @param originalValues the original values
#' @param values the processed values
#' @param bridgeDF the data frame to bridge identifiers
#' @param qcFolderLocal the folder where to save the results
#' 
#' @return returns the updated values
exportCurves <- function(
  originalValues,
  values,
  bridgeDF,
  qcFolderLocal
) {
  
  for (i in 1:nrow(values)) {
    
    exportCurvesForKid(
      originalValues = originalValues,
      i = i,
      values = values,
      bridgeDF = bridgeDF,
      qcFolderLocal = qcFolderLocal
    )
    
  }
}


#' Exports the groeth curve of a given kid.
#' 
#' @param originalValues the original values
#' @param i the line where the kid is found
#' @param values the processed values
#' @param bridgeDF the data frame to bridge identifiers
#' @param qcFolderLocal the folder where to save the results
#' 
#' @return returns the updated values
exportCurvesForKid <- function(
  originalValues,
  i,
  values,
  bridgeDF,
  qcFolderLocal
) {
  
  timePoints <- c("Birth", "6 w", "3 m", "6 m", "8 m", "1 y", "1.5 y", "2 y", "3 y", "5 y", "7 y", "8 y", "14 y")
  
  print(paste0(i, " in ", nrow(values)))
  
  logI <- values$log[i]
  
  if (logI != "") {
    
    child_idI <- originalValues$child_id[i]
    
    if (child_idI != values$child_id[i] || child_idI != bridgeDF$child_id[i]) {
      stop("id mismatch")
    }
    
    dummyIdI <- bridgeDF$dummyId[i]
    
    originalLength <- as.numeric(originalValues[i, length_columns])
    
    newLength <- as.numeric(values[i, length_columns])
    
    diff <- character(length(originalLength))
    diff[] <- "Unchanged"
    changedIs <- !is.na(originalLength) & is.na(newLength)
    diff[changedIs] <- "Outlier"
    changedIs <- is.na(originalLength) & !is.na(newLength)
    diff[changedIs] <- "Imputed"
    changedIs <- !is.na(originalLength) & !is.na(newLength) & originalLength != newLength
    diff[changedIs] <- "Corrected"
    
    pointsDF <- data.frame(
      timePoint = c(timePoints[diff == "Corrected"], timePoints[diff == "Outlier"], timePoints[diff != "Outlier"]),
      length = c(originalLength[diff == "Corrected"], originalLength[diff == "Outlier"], newLength[diff != "Outlier"]),
      category = c(rep("Original", sum(diff == "Corrected")), rep("Outlier", sum(diff == "Outlier")), diff[diff != "Outlier"]),
      stringsAsFactors = F
    ) %>%
      mutate(
        x = factor(timePoint, levels = timePoints),
        col = factor(category, levels = c("Unchanged", "Original", "Outlier", "Imputed", "Corrected"))
      )
    
    lengthDF <- data.frame(
      x = factor(timePoints, levels = timePoints),
      originalLength,
      newLength
    )
    
    colors <- c()
    if ("Unchanged" %in% pointsDF$col) {
      colors <- c(colors, "black")
    }
    if ("Original" %in% pointsDF$col) {
      colors <- c(colors, "grey")
    }
    if ("Outlier" %in% pointsDF$col) {
      colors <- c(colors, "red3")
    }
    if ("Imputed" %in% pointsDF$col) {
      colors <- c(colors, "green3")
    }
    if ("Corrected" %in% pointsDF$col) {
      colors <- c(colors, "blue3")
    }
    
    lengthPlot <- ggplot() + theme_bw() +
      geom_line(data = lengthDF, mapping = aes(x = x, y = originalLength, group = 1), col = "grey40", linetype = "dotted") +
      geom_line(data = lengthDF, mapping = aes(x = x, y = newLength, group = 1), col = "darkblue", linetype = "solid") +
      geom_point(data = pointsDF, mapping = aes(x = x, y = length, col = col, group = col)) +
      scale_color_manual(name = "Category", values = colors) +
      ylab("Length [cm]") +
      ggtitle(paste0(dummyIdI, " Log:", logI)) +
      theme(
        axis.title.x = element_blank()
      )
    
    originalWeight <- as.numeric(originalValues[i, weight_columns])
    
    newWeight <- as.numeric(values[i, weight_columns])
    
    diff <- character(length(originalWeight))
    diff[] <- "Unchanged"
    changedIs <- !is.na(originalWeight) & is.na(newWeight)
    diff[changedIs] <- "Outlier"
    changedIs <- is.na(originalWeight) & !is.na(newWeight)
    diff[changedIs] <- "Imputed"
    changedIs <- !is.na(originalWeight) & !is.na(newWeight) & originalWeight != newWeight
    diff[changedIs] <- "Corrected"
    
    pointsDF <- data.frame(
      timePoint = c(timePoints[diff == "Corrected"], timePoints[diff == "Outlier"], timePoints[diff != "Outlier"]),
      weight = c(originalWeight[diff == "Corrected"], originalWeight[diff == "Outlier"], newWeight[diff != "Outlier"]),
      category = c(rep("Original", sum(diff == "Corrected")), rep("Outlier", sum(diff == "Outlier")), diff[diff != "Outlier"]),
      stringsAsFactors = F
    ) %>%
      mutate(
        x = factor(timePoint, levels = timePoints),
        col = factor(category, levels = c("Unchanged", "Original", "Outlier", "Imputed", "Corrected"))
      )
    weightDF <- data.frame(
      x = factor(timePoints, levels = timePoints),
      originalWeight,
      newWeight
    )
    
    colors <- c()
    if ("Unchanged" %in% pointsDF$col) {
      colors <- c(colors, "black")
    }
    if ("Original" %in% pointsDF$col) {
      colors <- c(colors, "grey")
    }
    if ("Outlier" %in% pointsDF$col) {
      colors <- c(colors, "red3")
    }
    if ("Imputed" %in% pointsDF$col) {
      colors <- c(colors, "green3")
    }
    if ("Corrected" %in% pointsDF$col) {
      colors <- c(colors, "blue")
    }
    
    
    weightPlot <- ggplot() + theme_bw() +
      geom_line(data = weightDF, mapping = aes(x = x, y = originalWeight, group = 1), col = "grey40", linetype = "dotted") +
      geom_line(data = weightDF, mapping = aes(x = x, y = newWeight, group = 1), col = "blue3", linetype = "solid") +
      geom_point(data = pointsDF, mapping = aes(x = x, y = weight, col = col, group = col)) +
      scale_color_manual(name = "Category", values = colors) +
      ylab("Weight [kg]") +
      theme(
        axis.title.x = element_blank()
      )
    
    plotFile <- file.path(qcFolderLocal, paste0(dummyIdI, "_length.png"))
    png(plotFile, width = 900, height = 600)
    grid.draw(lengthPlot)
    dummy <- dev.off()
    
    plotFile <- file.path(qcFolderLocal, paste0(dummyIdI, "_weight.png"))
    png(plotFile, width = 900, height = 600)
    grid.draw(weightPlot)
    dummy <- dev.off()
    
  }
}

