#' Returns the number of children, mothers, and fathers for the given phenotype.
#' 
#' @param values the association data frame containing the phenotypes of all children
#' @param phenoName the name of the phenotype
#' 
#' @return The number of children, mothers, and fathers in a vector.
getNValues <- function(
    values,
    phenoName
) {
  
  result <- c(
    sum(!is.na(values$child_sentrix_id) & !is.na(values[[phenoName]])),
    sum(!is.na(values$mother_sentrix_id) & !is.na(values[[phenoName]])),
    sum(!is.na(values$father_sentrix_id) & !is.na(values[[phenoName]]))
  )
  
  return(result)
  
}


#' Builds a plot comparing two phenotypes and plots in the docs folder.
#' 
#' @param df the data frame containing the phenos
#' @param pheno1 the name of the column in df to plot on the x axis
#' @param pheno2 the name of the column in df to plot on the y axis
#' @param labelX the label for the x axis
#' @param labelY the label for the y axis
#' @param plotFile the file where to draw the plot
plotPhenos <- function(
    df,
    pheno1, 
    pheno2, 
    labelX, 
    labelY,
    plotFile
) {
  
  if (!pheno1 %in% names(df)) {
    stop(paste0("Pheno ", pheno1, " not found in values: ", paste(names(df), collapse = ", ")))
  }
  if (!pheno2 %in% names(df)) {
    stop(paste0("Pheno ", pheno2, " not found in values: ", paste(names(df), collapse = ", ")))
  }
  
  nX <- length(unique(df[[pheno1]][!is.na(df[[pheno1]]) & !is.infinite(df[[pheno1]])]))
  nY <- length(unique(df[[pheno2]][!is.na(df[[pheno2]]) & !is.infinite(df[[pheno2]])]))
  
  if (nX == 0 && nY == 0) {
    
    phenoError(
      message = glue("No values found for {pheno1} and {pheno2})"),
      plotFile = plotFile
    )
    
  } else if (nX == 0) {
    
    phenoError(
      message = glue("No values found for {pheno2})"),
      plotFile = plotFile
    )
    
  } else if (nY == 0) {
    
    phenoError(
      message = glue("No values found for {pheno2})"),
      plotFile = plotFile
    )
    
  } else if (nX > 10 && nY > 10) {
    
    phenoScatter(
      df = df,
      pheno1 = pheno1, 
      pheno2 = pheno2, 
      labelX = labelX, 
      labelY = labelY,
      plotFile = plotFile
    )
    
  } else if (nY > 10) {
    
    phenoViolinX(
      df = df,
      pheno1 = pheno1, 
      pheno2 = pheno2, 
      labelX = labelX, 
      labelY = labelY,
      plotFile = plotFile
    )
    
  } else if (nX > 10) {
    
    phenoViolinY(
      df = df,
      pheno1 = pheno1, 
      pheno2 = pheno2, 
      labelX = labelX, 
      labelY = labelY,
      plotFile = plotFile
    )
    
  } else {
    
    phenoDiscrete(
      df = df,
      pheno1 = pheno1, 
      pheno2 = pheno2, 
      labelX = labelX, 
      labelY = labelY,
      plotFile = plotFile
    )
    
  }
  
}


#' Builds a scatter plot comparing two phenotypes and plots in the docs folder.
#' 
#' @param df the data frame containing the phenos
#' @param pheno1 the name of the column in df to plot on the x axis
#' @param pheno2 the name of the column in df to plot on the y axis
#' @param labelX the label for the x axis
#' @param labelY the label for the y axis
#' @param plotFile the file where to draw the plot
phenoViolinX <- function(
    df,
    pheno1, 
    pheno2, 
    labelX, 
    labelY,
    plotFile
) {
  
  plotDF <- data.frame(
    x = df[[pheno1]],
    y = df[[pheno2]],
    stringsAsFactors = F
  ) %>%
    filter(
      !is.infinite(x) & !is.na(x) & !is.infinite(y) & !is.na(y)
    ) %>%
    mutate(
      x = factor(x)
    )
  
  levels(plotDF$x) <- round(as.numeric(levels(plotDF$x)), digits = 2)
  
  minY <- min(plotDF$y) - 0.05 * (max(plotDF$y) - min(plotDF$y))
  maxY <- max(plotDF$y) + 0.05 * (max(plotDF$y) - min(plotDF$y))
  
  # Build the scatter plot
  
  violinPlot <- ggplot(
    data = plotDF
  ) +
    geom_violin(
      mapping = aes(
        x = x,
        y = y
      ),
      fill = "grey80",
      alpha = 0.8
    ) +
    geom_boxplot(
      mapping = aes(
        x = x,
        y = y
      ),
      fill = "white",
      width = 0.1
    ) +
    scale_x_discrete(
      name = labelX
    ) +
    scale_y_continuous(
      name = labelY,
      limits = c(minY, maxY),
      expand = c(0, 0)
    ) +
    theme(
      legend.position = "none"
    )
  
  
  # Build the bar and density plot
  
  xBarPlot <- ggplot(
    data = plotDF
  ) + theme_minimal() + 
    geom_bar(
      mapping = aes(
        x = x
      ),
      fill = "black",
      col = "black",
      alpha = 0.1
    ) +
    scale_y_continuous(
      expand = c(0, 0)
    ) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank()
    )
  
  yDensityPlot <- ggplot(
    data = plotDF
  ) + theme_minimal() + 
    geom_density(
      mapping = aes(
        x = y
      ),
      fill = "black",
      alpha = 0.1
    ) +
    scale_x_continuous(
      limits = c(minY, maxY),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      expand = c(0, 0)
    ) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank()
    ) + 
    coord_flip()
  
  
  # Make grobs from plots
  
  violinGrob <- ggplotGrob(violinPlot)
  xBarGrob <- ggplotGrob(xBarPlot)
  yDensityGrob <- ggplotGrob(yDensityPlot)
  
  
  # Insert the bars and densities as new row and column in the scatter grob
  
  mergedGrob <- rbind(violinGrob[1:6, ], xBarGrob[7, ], violinGrob[7:nrow(violinGrob), ], size = "last")
  mergedGrob$heights[7] <- unit(0.15, "null")
  
  yDensityGrob <- gtable_add_rows(
    x = yDensityGrob, 
    heights = unit(rep(0, nrow(mergedGrob) - nrow(yDensityGrob)), "null"), 
    pos = 0
  )
  
  mergedGrob <- cbind(mergedGrob[, 1:5], yDensityGrob[, 5], mergedGrob[, 6:ncol(mergedGrob)], size = "first")
  mergedGrob$widths[6] <- unit(0.15, "null")
  
  
  # Plot
  
  png(
    filename = plotFile,
    width = 800,
    height = 600
  )
  grid.draw(mergedGrob)
  none <- dev.off()
  
}


#' Builds a scatter plot comparing two phenotypes and plots in the docs folder.
#' 
#' @param df the data frame containing the phenos
#' @param pheno1 the name of the column in df to plot on the x axis
#' @param pheno2 the name of the column in df to plot on the y axis
#' @param labelX the label for the x axis
#' @param labelY the label for the y axis
#' @param plotFile the file where to draw the plot
phenoViolinY <- function(
    df,
    pheno1, 
    pheno2, 
    labelX, 
    labelY,
    plotFile
) {
  
  plotDF <- data.frame(
    x = df[[pheno1]],
    y = df[[pheno2]],
    stringsAsFactors = F
  ) %>%
    filter(
      !is.infinite(x) & !is.na(x) & !is.infinite(y) & !is.na(y)
    ) %>%
    mutate(
      y = factor(y)
    )
  
  levels(plotDF$y) <- round(as.numeric(levels(plotDF$y)), digits = 2)
  
  minX <- min(plotDF$x) - 0.05 * (max(plotDF$x) - min(plotDF$x))
  maxX <- max(plotDF$x) + 0.05 * (max(plotDF$x) - min(plotDF$x))
  
  # Build the scatter plot
  
  violinPlot <- ggplot(
    data = plotDF
  ) +
    geom_violin(
      mapping = aes(
        x = x,
        y = y
      ),
      fill = "grey80",
      alpha = 0.8
    ) +
    geom_boxplot(
      mapping = aes(
        x = x,
        y = y
      ),
      fill = "white",
      width = 0.1
    ) +
    scale_x_continuous(
      name = labelX,
      limits = c(minX, maxX),
      expand = c(0, 0)
    ) +
    scale_y_discrete(
      name = labelY
    ) +
    theme(
      legend.position = "none"
    )
  
  
  # Build the bar and density plot
  
  xDensityPlot <- ggplot(
    data = plotDF
  ) + theme_minimal() + 
    geom_density(
      mapping = aes(
        x = x
      ),
      fill = "black",
      alpha = 0.1
    ) +
    scale_x_continuous(
      limits = c(minX, maxX),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      expand = c(0, 0)
    ) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank()
    )
  
  yBarPlot <- ggplot(
    data = plotDF
  ) + theme_minimal() + 
    geom_bar(
      mapping = aes(
        x = y
      ),
      fill = "black",
      col = "black",
      alpha = 0.1
    ) +
    scale_y_continuous(
      expand = c(0, 0)
    ) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank()
    ) + 
    coord_flip()
  
  
  # Make grobs from plots
  
  violinGrob <- ggplotGrob(violinPlot)
  xDensityGrob <- ggplotGrob(xDensityPlot)
  yBarGrob <- ggplotGrob(yBarPlot)
  
  
  # Insert the bars and densities as new row and column in the scatter grob
  
  mergedGrob <- rbind(violinGrob[1:6, ], xDensityGrob[7, ], violinGrob[7:nrow(violinGrob), ], size = "last")
  mergedGrob$heights[7] <- unit(0.15, "null")
  
  yBarGrob <- gtable_add_rows(
    x = yBarGrob, 
    heights = unit(rep(0, nrow(mergedGrob) - nrow(yBarGrob)), "null"), 
    pos = 0
  )
  
  mergedGrob <- cbind(mergedGrob[, 1:5], yBarGrob[, 5], mergedGrob[, 6:ncol(mergedGrob)], size = "first")
  mergedGrob$widths[6] <- unit(0.15, "null")
  
  
  # Plot
  
  png(
    filename = plotFile,
    width = 800,
    height = 600
  )
  grid.draw(mergedGrob)
  none <- dev.off()
  
}


#' Builds a scatter plot comparing two phenotypes and plots in the docs folder.
#' 
#' @param df the data frame containing the phenos
#' @param pheno1 the name of the column in df to plot on the x axis
#' @param pheno2 the name of the column in df to plot on the y axis
#' @param labelX the label for the x axis
#' @param labelY the label for the y axis
#' @param plotFile the file where to draw the plot
phenoDiscrete <- function(
    df,
    pheno1, 
    pheno2, 
    labelX, 
    labelY,
    plotFile
) {
  
  plotDF <- data.frame(
    x = df[[pheno1]],
    y = df[[pheno2]],
    stringsAsFactors = F
  ) %>%
    filter(
      !is.infinite(x) & !is.na(x) & !is.infinite(y) & !is.na(y)
    ) %>%
    mutate(
      x = factor(x),
      y = factor(y) 
    ) %>% 
    group_by(
      x, y
    ) %>% 
    mutate(
      n = n(),
      percent = round(100 * n()/nrow(df), digits = 2),
      label = paste0(n, "\n(", percent, " %)")
    ) %>% 
    ungroup()
  
  levels(plotDF$x) <- round(as.numeric(levels(plotDF$x)), digits = 2)
  levels(plotDF$y) <- round(as.numeric(levels(plotDF$y)), digits = 2)
  
  
  # Build the heatmap
  
  heatmapPlot <- ggplot(
    data = plotDF
  ) +
    geom_raster(
      mapping = aes(
        x = x,
        y = y,
        fill = percent
      )
    ) +
    geom_text(
      mapping = aes(
        x = x,
        y = y,
        label = label
      )
    ) +
    scale_fill_scico(
      palette = "grayC",
      begin = 0,
      end = 0.5
    ) +
    scale_x_discrete(
      name = labelX,
      expand = c(0, 0)
    ) +
    scale_y_discrete(
      name = labelY,
      expand = c(0, 0)
    ) +
    theme(
      legend.position = "none",
      panel.grid = element_blank()
    )
  
  
  # Build the bar plots
  
  xBarPlot <- ggplot(
    data = plotDF
  ) + theme_minimal() + 
    geom_bar(
      mapping = aes(
        x = x
      ),
      fill = "black",
      col = "black",
      alpha = 0.1
    ) +
    scale_y_continuous(
      expand = c(0, 0)
    ) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank()
    )
  
  yBarPlot <- ggplot(
    data = plotDF
  ) + theme_minimal() + 
    geom_bar(
      mapping = aes(
        x = y
      ),
      fill = "black",
      col = "black",
      alpha = 0.1
    ) +
    scale_y_continuous(
      expand = c(0, 0)
    ) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank()
    ) + 
    coord_flip()
  
  
  # Make grobs from plots
  
  heatmapGrob <- ggplotGrob(heatmapPlot)
  xBarGrob <- ggplotGrob(xBarPlot)
  yBarGrob <- ggplotGrob(yBarPlot)
  
  
  # Insert the bars and densities as new row and column in the scatter grob
  
  mergedGrob <- rbind(heatmapGrob[1:6, ], xBarGrob[7, ], heatmapGrob[7:nrow(heatmapGrob), ], size = "last")
  mergedGrob$heights[7] <- unit(0.15, "null")
  
  yBarGrob <- gtable_add_rows(
    x = yBarGrob, 
    heights = unit(rep(0, nrow(mergedGrob) - nrow(yBarGrob)), "null"), 
    pos = 0
  )
  
  mergedGrob <- cbind(mergedGrob[, 1:5], yBarGrob[, 5], mergedGrob[, 6:ncol(mergedGrob)], size = "first")
  mergedGrob$widths[6] <- unit(0.15, "null")
  
  
  # Plot
  
  png(
    filename = plotFile,
    width = 800,
    height = 600
  )
  grid.draw(mergedGrob)
  none <- dev.off()
  
}


#' Writes an error message in the docs folder.
#' 
#' @param message error message 
#' @param plotFile the file where to draw the plot
phenoError <- function(
    message,
    plotFile
) {
  
  messagePlot <- ggplot() +
    theme_void() +
    geom_text(
      mapping = aes(
        x = 0,
        y = 0,
        label = message
      ),
      col = "black"
    )
  
  png(
    filename = plotFile,
    width = 800,
    height = 600
  )
  grid.draw(messagePlot)
  none <- dev.off()
  
}


#' Builds a scatter plot comparing two phenotypes and plots in the docs folder.
#' 
#' @param df the data frame containing the phenos
#' @param pheno1 the name of the column in df to plot on the x axis
#' @param pheno2 the name of the column in df to plot on the y axis
#' @param labelX the label for the x axis
#' @param labelY the label for the y axis
#' @param plotFile the file where to draw the plot
phenoScatter <- function(
    df,
    pheno1, 
    pheno2, 
    labelX, 
    labelY,
    plotFile
) {
  
  plotDF <- data.frame(
    x = df[[pheno1]],
    y = df[[pheno2]],
    stringsAsFactors = F
  ) %>%
    filter(
      !is.infinite(x) & !is.na(x) & !is.infinite(y) & !is.na(y)
    ) %>%
    arrange(
      rev(abs(y - x))
    )
  
  minX <- min(plotDF$x) - 0.05 * (max(plotDF$x) - min(plotDF$x))
  maxX <- max(plotDF$x) + 0.05 * (max(plotDF$x) - min(plotDF$x))
  minY <- min(plotDF$y) - 0.05 * (max(plotDF$y) - min(plotDF$y))
  maxY <- max(plotDF$y) + 0.05 * (max(plotDF$y) - min(plotDF$y))
  
  # Build the scatter plot
  
  tryCatch(
    {
      
      scatterPlot <- ggplot(
        data = plotDF
      ) +
        geom_point(
          mapping = aes(
            x = x,
            y = y
          ),
          col = "black",
          alpha = 0.1
        ) + 
        geom_density_2d(
          mapping = aes(
            x = x,
            y = y
          ),
          col = "white",
        ) +
        scale_x_continuous(
          name = labelX,
          limits = c(minX, maxX),
          expand = c(0, 0)
        ) +
        scale_y_continuous(
          name = labelY,
          limits = c(minY, maxY),
          expand = c(0, 0)
        ) +
        theme(
          legend.position = "none"
        )
      
    }, error = function(error_condition) {
      
      scatterPlot <- ggplot(
        data = plotDF
      ) +
        geom_point(
          mapping = aes(
            x = x,
            y = y
          ),
          col = "black",
          alpha = 0.1
        ) + 
        scale_x_continuous(
          name = labelX,
          limits = c(minX, maxX),
          expand = c(0, 0)
        ) +
        scale_y_continuous(
          name = labelY,
          limits = c(minY, maxY),
          expand = c(0, 0)
        ) +
        theme(
          legend.position = "none"
        )
      
    }
  )
  
  
  # Build the density plots
  
  xDensityPlot <- ggplot(
    data = plotDF
  ) + theme_minimal() + 
    geom_density(
      mapping = aes(
        x = x
      ),
      fill = "black",
      alpha = 0.1
    ) +
    scale_x_continuous(
      limits = c(minX, maxX),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      expand = c(0, 0)
    ) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank()
    )
  
  yDensityPlot <- ggplot(
    data = plotDF
  ) + theme_minimal() + 
    geom_density(
      mapping = aes(
        x = y
      ),
      fill = "black",
      alpha = 0.1
    ) +
    scale_x_continuous(
      limits = c(minY, maxY),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      expand = c(0, 0)
    ) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank()
    ) + 
    coord_flip()
  
  
  # Make grobs from plots
  
  scatterGrob <- ggplotGrob(scatterPlot)
  xDensityGrob <- ggplotGrob(xDensityPlot)
  yDensityGrob <- ggplotGrob(yDensityPlot)
  
  
  # Insert the densities as new row and column in the scatter grob
  
  mergedGrob <- rbind(scatterGrob[1:6, ], xDensityGrob[7, ], scatterGrob[7:nrow(scatterGrob), ], size = "last")
  mergedGrob$heights[7] <- unit(0.15, "null")
  
  yDensityGrob <- gtable_add_rows(
    x = yDensityGrob, 
    heights = unit(rep(0, nrow(mergedGrob) - nrow(yDensityGrob)), "null"), 
    pos = 0
  )
  
  mergedGrob <- cbind(mergedGrob[, 1:5], yDensityGrob[, 5], mergedGrob[, 6:ncol(mergedGrob)], size = "first")
  mergedGrob$widths[6] <- unit(0.15, "null")
  
  
  # Plot
  
  png(
    filename = plotFile,
    width = 800,
    height = 600
  )
  grid.draw(mergedGrob)
  none <- dev.off()
  
}


#' Writes the documentation on a phenotype.
#' 
#' @param phenoName the name of the phenotype
#' @param zPhenoName the name of the standardized phenotype
#' @param phenoLabel the label to use for the phenotype in the docs
#' @param phenoUnit the unit to use for the phenotype in the docs
#' @param formula the formula to use for the model
#' @param sigmaFormula the formula to use for sigma
#' @param family the family of distribution fitted
#' @param values the phenotype data frame
#' @param docsFolder the folder where to write the docs
#' @param docsFile the file where to write the documentation on all phenotypes
#' @param phenoDocsFile the file where to write the documentation on this phenotype
writeDocs <- function(
    phenoName,
    zPhenoName,
    phenoLabel,
    phenoUnit,
    formula,
    sigmaFormula,
    family,
    values,
    docsFolder,
    docsFile
) {
  
  phenoDocsFolder <- file.path(docsFolder, zPhenoName)
  dir.create(phenoDocsFolder, showWarnings = F)
  phenoDocsFile <- file.path(phenoDocsFolder, paste0(zPhenoName, ".md"))
  phenoPlotsFolder <- file.path(phenoDocsFolder, "plots")
  dir.create(phenoPlotsFolder, showWarnings = F)
  
  write(x = paste0("- [", zPhenoName, "](", zPhenoName, "/", zPhenoName, ".md): ", phenoLabel, "\n\n"), file = docsFile, append = T)
  write(x = paste0("## ", phenoLabel, "\n"), file = phenoDocsFile, append = F)
  
  nValues <- getNValues(
    values = values,
    phenoName = phenoName
  )
  nZValues <- getNValues(
    values = values,
    phenoName = zPhenoName
  )
  
  write(x = paste0("| Name | # Children | # Mothers | # Fathers | # Total |"), file = phenoDocsFile, append = T)
  write(x = paste0("| ---- | ---------- | --------- | --------- | ------- |"), file = phenoDocsFile, append = T)
  write(x = paste0("| ", phenoName, " | ", nValues[1], " | ", nValues[2], " | ", nValues[3], " | ", sum(nValues), " |"), file = phenoDocsFile, append = T)
  write(x = paste0("| ", zPhenoName, " | ", nZValues[1], " | ", nZValues[2], " | ", nZValues[3], " | ", sum(nZValues), " |\n"), file = phenoDocsFile, append = T)
  
  write(x = paste0("- Formula: `", formula, "`"), file = phenoDocsFile, append = T)
  write(x = paste0("- Sigma formula: `", sigmaFormula, "`"), file = phenoDocsFile, append = T)
  write(x = paste0("- Distribution: `", family, "`"), file = phenoDocsFile, append = T)
  write(x = paste0("- Normalization: `centiles.pred` Z-scores"), file = phenoDocsFile, append = T)
  
  fileName <- paste0(zPhenoName, "_vs_", phenoName, "_child.png")
  plotFile <- file.path(phenoPlotsFolder, fileName)
  plotPhenos(
    df = values,
    pheno1 = phenoName,
    pheno2 = zPhenoName,
    labelX = paste0(phenoLabel, " [", phenoUnit, "]"),
    labelY = paste0(phenoLabel, " [Z-score]"),
    plotFile = plotFile
  )
  
  write(x = paste0("![](plots/", fileName, ")\n\n"), file = phenoDocsFile, append = T)
  
}


#' Writes the documentation on the comparison between phenotypes by plotting a reference phenotype against another phenotype.
#' 
#' @param values the phenotype data frame
#' @param refPhenoName the name of the reference phenotype
#' @param refZPhenoName the name of the standardized reference phenotype
#' @param refPhenoLabel the label of the refrence phenotype to use in the docs
#' @param refPhenoUnit the unit of the refrence phenotype
#' @param otherPhenoName the name of the phenotype to compare to
#' @param otherPhenoLabel the label of the refrence phenotype to use in the docs
#' @param otherPhenoUnit the unit of the refrence phenotype
writePhenocomparison <- function(
    values,
    refPhenoName,
    refZPhenoName,
    refPhenoLabel,
    refPhenoUnit,
    otherPhenoName,
    otherPhenoLabel,
    otherPhenoUnit,
    docsFolder
) {
  
  
  phenoDocsFolder <- file.path(docsFolder, refZPhenoName)
  phenoDocsFile <- file.path(phenoDocsFolder, paste0(refZPhenoName, ".md"))
  phenoPlotsFolder <- file.path(phenoDocsFolder, "plots")
  
  write(x = paste0("- ", refPhenoLabel, " vs. ", otherPhenoLabel, ":"), file = phenoDocsFile, append = T)
  
  fileName <- paste0(refPhenoName, "_vs_", otherPhenoName, ".png")
  plotFile <- file.path(phenoPlotsFolder, fileName)
  plotPhenos(
    df = values,
    pheno1 = otherPhenoName,
    pheno2 = refPhenoName,
    labelX = paste0(otherPhenoLabel, " [", otherPhenoUnit, "]"),
    labelY = paste0(refPhenoLabel, " [", refPhenoUnit, "]"),
    plotFile = plotFile
  )
  
  write(x = paste0("![](plots/", fileName, ")\n\n"), file = phenoDocsFile, append = T)
  
  fileName <- paste0(refZPhenoName, "_vs_", otherPhenoName, ".png")
  plotFile <- file.path(phenoPlotsFolder, fileName)
  plotPhenos(
    df = values,
    pheno1 = otherPhenoName,
    pheno2 = refZPhenoName,
    labelX = paste0(otherPhenoLabel, " [", otherPhenoUnit, "]"),
    labelY = paste0(refPhenoLabel, " [Z-score]"),
    plotFile = plotFile
  )
  
  write(x = paste0("![](plots/", fileName, ")\n\n"), file = phenoDocsFile, append = T)
  
}