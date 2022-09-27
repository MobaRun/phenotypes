##
#
# Functions used in the cleaning pipeline to process length and weight values. 
#
##


#' Returns the indexes where the length decreases or stagnates in a data frame. i for the index, j for the next non-NA index.
#' 
#' @param lengthValues the length values
#' 
#' @return the indexes of negative growth
getNegativeGrowthIndexes <- function(lengthValues) {
  
  negativeGrowthIndexesI <- c()
  negativeGrowthIndexesJ <- c()
  
  for(i in 1:(length(lengthValues)-1)) {
    
    if (!is.na(lengthValues[i])) {
      
      for (j in (i+1):length(lengthValues)) {
        
        if (!is.na(lengthValues[j])) {
          
          if (lengthValues[j] <= lengthValues[i]) {
            
            negativeGrowthIndexesI <- c(negativeGrowthIndexesI, i)
            negativeGrowthIndexesJ <- c(negativeGrowthIndexesJ, j)
            
          } else {
            
            break
            
          }
        }
      }
    }
  }
  
  result <- data.frame(
    i = negativeGrowthIndexesI, 
    j = negativeGrowthIndexesJ
  )
  
  return(result)
  
}


#' Returns a boolean indicating whether the given index is outlying the growth curve.
#' 
#' @param growthValues the normalized growth values
#' @param index the index on the curve
#' 
#' @return a boolean indicating whether the given index is outlying the growth curve
isOutlierGrowth <- function(
  growthValues, 
  index
) {
  
  if (index == length(growthValues) + 1) {
    
    index <- index - 1
    return(
      !is.na(growthValues[index]) && !is.na(growthValues[index - 1]) && 
        growthValues[index - 1] >= -1 && growthValues[index - 1] <= 1 && 
        (growthValues[index] > 2 || growthValues[index] < -2)
    )
    
  } else if (index == 1) {
    
    return(
      !is.na(growthValues[index]) && !is.na(growthValues[index + 1]) && 
        growthValues[index + 1] >= -1 && growthValues[index + 1] <= 1 && 
        (growthValues[index] > 2 || growthValues[index] < -2)
    )
  }
  
  return(
    !is.na(growthValues[index - 1]) && !is.na(growthValues[index]) && 
      (growthValues[index - 1] > 1 && growthValues[index] < -1
       || growthValues[index - 1] < -1 && growthValues[index] > 1)
  )
}


#' Interpolates a value based on a growth curve
#' 
#' @param originalIndex the original index to interpolate from
#' @param destinationIndex the destination index to interpolate to
#' @param measuredValues the measured values
#' @param referenceGrowthCurve the reference growth curve
#' 
#' @return returns the interpolated value
getInterpolationFromCurve <- function(
  originalIndex, 
  destinationIndex, 
  measuredValues, 
  referenceGrowthCurve
) {
  
  interpolatedValue <- measuredValues[originalIndex]
  
  if (destinationIndex > originalIndex) {
    
    for (i in originalIndex:(destinationIndex - 1)) {
      
      interpolatedValue <- interpolatedValue * (2 ^ referenceGrowthCurve[i])
      
    }
  } else if (destinationIndex < originalIndex) {
    
    for (i in (originalIndex - 1):destinationIndex) {
      
      interpolatedValue <- interpolatedValue / (2 ^ referenceGrowthCurve[i])
      
    }
  }
  
  return(interpolatedValue)
  
}


#' Imputes a value based on the interpolation from the nearest neighbors on the given growth curve
#' 
#' @param index the index where a value is needed
#' @param measuredValues the measured values
#' @param referenceGrowthCurve the reference growth curve
#' @param nNeighbors the number of neighbors
#' 
#' @return the imputed value
getImputationFromCurve <- function(
  index, 
  measuredValues, 
  referenceGrowthCurve, 
  nNeighbors
) {
  
  interpolatedValues <- c()
  weights <- c()
  nBelow <- floor(nNeighbors/2)
  
  if (index > 1) {
    
    for (i in (index - 1):1) {
      
      if (!is.na(measuredValues[i])) {
        
        interpolatedValue <-
          getInterpolationFromCurve(
            originalIndex = i, 
            destinationIndex = index, 
            measuredValues = measuredValues, 
            referenceGrowthCurve = referenceGrowthCurve
          )
        
        interpolatedValues <- c(interpolatedValues, interpolatedValue)
        
        weight <- 1 / (index - i)
        weights <- c(weights, weight)
        
        if (length(interpolatedValues) == nBelow) {
          
          break
          
        }
      }
    }
  }
  
  if (index < length(measuredValues)) {
    
    for (i in (index + 1):length(measuredValues)) {
      
      if (!is.na(measuredValues[i])) {
        interpolatedValue <-
          getInterpolationFromCurve(
            originalIndex = i, 
            destinationIndex = index, 
            measuredValues = measuredValues, 
            referenceGrowthCurve = referenceGrowthCurve)
        
        interpolatedValues <-
          c(interpolatedValues, interpolatedValue)
        
        weight <- 1 / (i - index)
        weights <- c(weights, weight)
        
        if (length(interpolatedValues) == nNeighbors) {
          
          break
          
        }
      }
    }
  }
  
  if (length(interpolatedValues) < nNeighbors && index > 1) {
    
    for (i in (index - 1):1) {
      
      if (!is.na(measuredValues[i])) {
        
        interpolatedValue <-
          getInterpolationFromCurve(
            originalIndex = i, 
            destinationIndex = index, 
            measuredValues = measuredValues, 
            referenceGrowthCurve = referenceGrowthCurve)
        
        interpolatedValues <- c(interpolatedValues, interpolatedValue)
        
        weight <- 1 / (index - i)
        weights <- c(weights, weight)
        
        if (length(interpolatedValues) == nNeighbors) {
          
          break
          
        }
      }
    }
  }
  
  interpolatedValue <- weighted.mean(x = interpolatedValues, w = weights, na.rm = T)
  
  return(interpolatedValue)
  
}


#' Imputes a value based on the interpolation from the nearest neighbors on the given growth curve
#' 
#' @param index the index where a value is needed
#' @param measuredValues the measured values
#' @param populationGrowth the population growth data frame
#' @param sex the sex of the individual
#' @param pheno the pheno to extract from the population growth data frame
#' @param nNeighbors the number of neighbors
#' 
#' @return the imputed value
getImputationFromCurveBySex <- function(
  index, 
  measuredValues, 
  populationGrowth,
  sex,
  pheno,
  longitudinalCategory,
  nNeighbors
) {
  
  if (
    !is.na(measuredValues[index]) |
    longitudinalCategory != "longitudinal" &
    (longitudinalCategory != "early" | index > 7)
  ) {
    return(measuredValues[index])
  }
  
  nBefore <- 0
  for (i in 1:index) {
    
    if(!is.na(measuredValues[i])) {
      
      nBefore <- nBefore + 1
      
      if (nBefore >= minValidValuesBeforeFirstImputed) {
        
        break
        
      }
    }
  }
  if (nBefore < minValidValuesBeforeFirstImputed) {
    return(NA)
  }
  
  nAfter <- 0
  for (i in index:length(measuredValues)) {
    
    if(!is.na(measuredValues[i])) {
      
      nAfter <- nAfter + 1
      
      if (nAfter >= minValidValuesAfterLastImputed) {
        
        break
        
      }
    }
  }
  if (nAfter < minValidValuesAfterLastImputed) {
    return(NA)
  }
  
  
  referenceGrowthCurve <- getPopulationCurve(
    populationDF = populationGrowth, 
    quantileType = "median", 
    ageIs = 0:10, 
    sex = sex, 
    pheno = pheno
  )
  
  return(
    getImputationFromCurve(
      index = index, 
      measuredValues = measuredValues, 
      referenceGrowthCurve = referenceGrowthCurve, 
      nNeighbors = nNeighbors
    )
  )
  
}

#' Creates an imputed curve based on the imputation of missing values using the 3 nearest neighbours with no missing values.
#' 
#' @param measuredValues the values
#' @param referenceGrowthCurve the reference growth curve
#' 
#' @return the imputed curve
getImputedCurve <- function(
  measuredValues, 
  referenceGrowthCurve
) {
  
  imputedCurve <- measuredValues
  
  for (i in 1:length(measuredValues)) {
    
    if (is.na(measuredValues[i])) {
      
      imputedValue <- getImputationFromCurve(
        index = i, 
        measuredValues = measuredValues, 
        referenceGrowthCurve = referenceGrowthCurve, 
        nNeighbors = 3
      )
      
      imputedCurve[i] <- imputedValue
      
    }
  }
  
  return(imputedCurve)
  
}


#' Smoothes a curve based on the nearest neighbours.
#' 
#' @param measuredValues the values
#' @param referenceGrowthCurve the reference growth curve
#' @param nNeighbors the number of neighbors to consider
#' 
#' @return returns the smoothed curve
getSmoothedCurve <- function(
  measuredValues, 
  referenceGrowthCurve, 
  nNeighbors
) {
  
  smoothedCurve <- numeric(length(measuredValues))
  
  for (i in 1:length(measuredValues)) {
    
    imputedValue <- getImputationFromCurve(
      index = i, 
      measuredValues = measuredValues, 
      referenceGrowthCurve = referenceGrowthCurve, 
      nNeighbors = 4
    )
    
    smoothedCurve[i] <- imputedValue
    
  }
  
  return(smoothedCurve)
  
}

#' Corrects a growth curve by outlier removal and inference
#' 
#' @param child_id the child id
#' @param lengthValues the length values
#' @param sex the sex of the individual
#' @param growthCurves the reference growth curves
#' @param populationAnchors the population anchor values
#' 
#' @return the corrected growth curve
correctGrowthCurve <- function(
  child_id, 
  lengthValues, 
  sex, 
  growthCurves, 
  populationAnchors
) {
  
  if (sum(!is.na(lengthValues)) <= 1) {
    return(NULL)
  }
  
  # Get indexes of negative growth
  outliers <- getNegativeGrowthIndexes(lengthValues)
  
  if (nrow(outliers) == 0) {
    return(NULL)
  }
  
  # Initiate new values
  
  newLengthValues <- lengthValues
  
  
  # get reference curves and length
  
  referenceGrowthCurve <- getPopulationCurve(
    populationDF = growthCurves, 
    quantileType = "median", 
    ageIs = 0:10, 
    sex = sex, 
    pheno = "growth"
  )
  lows <- getPopulationCurve(
    populationDF = populationAnchors, 
    quantileType = "low", 
    ageIs = 0:11, 
    sex = sex, 
    pheno = "length"
  )
  medians <- getPopulationCurve(
    populationDF = populationAnchors, 
    quantileType = "median", 
    ageIs = 0:11, 
    sex = sex, 
    pheno = "length"
  )
  highs <- getPopulationCurve(
    populationDF = populationAnchors, 
    quantileType = "high", 
    ageIs = 0:11, 
    sex = sex, 
    pheno = "length"
  )
  
  
  # Correct the decreasing values
  
  while (nrow(outliers) > 0) {
    
    tempLengthValues <- newLengthValues
    outlierFound <- FALSE
    
    for (k in 1:nrow(outliers)) {
      
      # Try to find reference values outside the incorrect areas
      i <- outliers$i[k]
      j <- outliers$j[k]
      refValuesAtI <- c()
      refValuesAtJ <- c()
      
      for (l in 1:length(newLengthValues)) {
        
        if (!is.na(newLengthValues[l])) {
          
          if (!l %in% outliers$i) {
            
            interpolatedValue <- getInterpolationFromCurve(
              originalIndex = l, 
              destinationIndex = i, 
              measuredValues = newLengthValues, 
              referenceGrowthCurve = referenceGrowthCurve
            )
            
            refValuesAtI <- c(refValuesAtI, interpolatedValue)
            
          }
          
          if (!l %in% outliers$j) {
            
            interpolatedValue <- getInterpolationFromCurve(
              originalIndex = l, 
              destinationIndex = j, 
              measuredValues = newLengthValues, 
              referenceGrowthCurve = referenceGrowthCurve
            )
            
            refValuesAtJ <- c(refValuesAtJ, interpolatedValue)
            
          }
        }
      }
      
      if (length(refValuesAtI) <= 2) {
        
        # Try to find reference values from the background distribution
        for (l in 1:length(newLengthValues)) {
          
          if (
            !is.na(newLengthValues[l]) && 
            newLengthValues[l] > lows[l] && 
            newLengthValues[l] < highs[l]
          ) {
            
            interpolatedValue <- getInterpolationFromCurve(
              originalIndex = l, 
              destinationIndex = i, 
              measuredValues = newLengthValues, 
              referenceGrowthCurve = referenceGrowthCurve
            )
            refValuesAtI <- c(refValuesAtI, interpolatedValue)
            
            interpolatedValue <- getInterpolationFromCurve(
              originalIndex = l, 
              destinationIndex = j, 
              measuredValues = newLengthValues, 
              referenceGrowthCurve = referenceGrowthCurve
            )
            refValuesAtJ <- c(refValuesAtJ, interpolatedValue)
            
          }
        }
      }
      
      if (length(refValuesAtI) == 0) {
        stop(paste("No reference found for", child_id))
      }
      
      # Infer the values of obvious outliers from the reference values
      interpolatedValue <- getInterpolationFromCurve(
        originalIndex = j, 
        destinationIndex = i, 
        measuredValues = newLengthValues, 
        referenceGrowthCurve = referenceGrowthCurve
      )
      refValuesAtI <- c(refValuesAtI, interpolatedValue)
      
      interpolatedValue <- getInterpolationFromCurve(
        originalIndex = i, 
        destinationIndex = j, 
        measuredValues = newLengthValues, 
        referenceGrowthCurve = referenceGrowthCurve
      )
      refValuesAtJ <- c(refValuesAtJ, interpolatedValue)
      
      maxAtI <- max(refValuesAtI)
      
      if (is.na(maxAtI)) {
        stop(paste("No max value found at", i, "for", child_id))
      }
      
      if (newLengthValues[i] > maxAtI) {
        
        tempLengthValues[i] <- median(refValuesAtI)
        outlierFound <- TRUE
        
      } 
      
      minAtJ <- min(refValuesAtJ)
      
      if (is.na(minAtJ)) {
        stop(paste("No min value found at", i, "for", child_id))
      }
      
      if (newLengthValues[j] < minAtJ) {
        
        tempLengthValues[j] <- median(refValuesAtJ)
        outlierFound <- TRUE
        
      }
    }
    
    if (!outlierFound) {
      
      # Perform a curve smoothing
      
      for (k in 1:nrow(outliers)) {
        
        # Try to find reference values outside the incorrect areas
        i <- outliers$i[k]
        j <- outliers$j[k]
        refValuesAtI <- c()
        refValuesAtJ <- c()
        
        for (l in 1:length(newLengthValues)) {
          
          if (!is.na(newLengthValues[l]) && !l %in% outliers$i && !l %in% outliers$j) {
            
            interpolatedValue <- getInterpolationFromCurve(
              originalIndex = l, 
              destinationIndex = i, 
              measuredValues = newLengthValues, 
              referenceGrowthCurve = referenceGrowthCurve
            )
            refValuesAtI <- c(refValuesAtI, interpolatedValue)
            
            interpolatedValue <- getInterpolationFromCurve(
              originalIndex = l, 
              destinationIndex = j, 
              measuredValues = newLengthValues, 
              referenceGrowthCurve = referenceGrowthCurve)
            refValuesAtJ <- c(refValuesAtJ, interpolatedValue)
            
          }
        }
        
        if (length(refValuesAtI) == 0) {
          
          # Try to find reference values from the background distribution
          for (l in 1:length(newLengthValues)) {
            
            if (!is.na(newLengthValues[l]) && newLengthValues[l] > lows[l] && newLengthValues[l] < highs[l]) {
              
              interpolatedValue <- getInterpolationFromCurve(
                originalIndex = l, 
                destinationIndex = i, 
                measuredValues = newLengthValues, 
                referenceGrowthCurve = referenceGrowthCurve
              )
              refValuesAtI <- c(refValuesAtI, interpolatedValue)
              
              interpolatedValue <- getInterpolationFromCurve(
                originalIndex = l, 
                destinationIndex = j, 
                measuredValues = newLengthValues, 
                referenceGrowthCurve = referenceGrowthCurve
              )
              refValuesAtJ <- c(refValuesAtJ, interpolatedValue)
              
            }
          }
        }
        
        if (length(refValuesAtI) == 0) {
          stop(paste("No reference found for curve correction for", child_id))
        }
        
        # Smooth the most extreme values
        while (tempLengthValues[j] <= tempLengthValues[i]) {
          
          medianAtI <- median(refValuesAtI)
          medianAtJ <- median(refValuesAtJ)
          
          if (medianAtJ < medianAtI) {
            stop(paste("Decreasing reference curve for", child_id))
          }
          
          deviationAtI <- tempLengthValues[i] - medianAtI
          deviationAtJ <- tempLengthValues[j] - medianAtJ
          
          if (abs(deviationAtI) < 0.0001 && abs(deviationAtJ) < 0.0001) {
            stop(paste("Smoothing does not converge for:", child_id))
          }
          
          if (abs(deviationAtI) > abs(deviationAtJ)) {
            
            tempLengthValues[i] <- medianAtI + 0.9 * deviationAtI
            
          } else {
            
            tempLengthValues[j] <- medianAtJ + 0.9 * deviationAtJ
            
          }
        }
      }
    }
    
    if (all(newLengthValues == tempLengthValues)) {
      
      stop(paste("No correction performed to the curve for", child_id))
      
    }
    
    newLengthValues <- tempLengthValues
    
    outliers <- getNegativeGrowthIndexes(newLengthValues)
    
  }
  
  return(newLengthValues)
  
}

#' Corrects points of a growth curve presenting decrease similar to the growth curve.
#' 
#' @param oldLengthValues the old length values
#' @param newLengthValues the new length values
#' @param weightValues the weight values
#' @param sex the sex of the individual
#' @param growthCurves the reference growth curves
#' 
#' @return the corrected weight curve
correctWeightIncreaseCurve <- function(
  oldLengthValues, 
  newLengthValues,
  weightValues,
  sex, 
  growthCurves
) {
  
  referenceGrowthCurve <- getPopulationCurve(
    populationDF = growthCurves, 
    quantileType = "median", 
    ageIs = 0:10, 
    sex = sex, 
    pheno = "growth"
  )
  referenceWeightIncreaseCurve <- getPopulationCurve(
    populationDF = growthCurves, 
    quantileType = "median", 
    ageIs = 0:10, 
    sex = sex, 
    pheno = "weightIncrease"
  )
  
  newWeightValues <- weightValues
  
  for (i in 1:length(oldLengthValues)) {
    
    if (
      !is.na(oldLengthValues[i]) &&
      !is.na(newLengthValues[i]) &&
      !is.na(weightValues[i]) &&
      oldLengthValues[i] != newLengthValues[i]
    ) {
      
      imputedLengthValue <- getImputationFromCurve(
        index = i, 
        measuredValues = oldLengthValues, 
        referenceGrowthCurve = referenceGrowthCurve, 
        nNeighbors = 4
      )
      
      oldLengthD <- oldLengthValues[i] - imputedLengthValue
      newLengthD <- newLengthValues[i] - imputedLengthValue
      
      if (abs(newLengthD) < abs(oldLengthD)) {
        
        imputedWeightValue <- getImputationFromCurve(
          index = i, 
          measuredValues = weightValues, 
          referenceGrowthCurve = referenceWeightIncreaseCurve, 
          nNeighbors = 4
        )
        
        weightD <- weightValues[i] - imputedWeightValue
        
        if (sign(weightD * oldLengthD) == 1) {
          
          newWeightValues[i] <- imputedWeightValue + weightD * newLengthD / oldLengthD
          
        }
      }
    }
  }
  
  return(newWeightValues)
  
}


#' Returns the low, median, and high growth and weight increase at the different time points in a data frame.
#' 
#' @param methodValues the values to inspect
#' @param phenoColumns the columns containing the phenotypes from the variables mapping
#' @param categoryColumns the columns containing the category for the phenos
#' @param ages the ages to get values from
#' 
#' @return returns the population quantiles
getPopulationQuantiles <- function(
  methodValues, 
  phenoColumns, 
  categoryColumns,
  ages
) {
  
  k <- 1
  n <- 3 * 3 * length(phenoColumns) * length(ages)
  
  quantilesDF <- data.frame(
    quantile = numeric(n),
    quantileType = character(n),
    age = numeric(n),
    sex = character(n),
    pheno = character(n),
    stringsAsFactors = F
  )
  
  for (phenoI in 1:length(phenoColumns)) {
    
    pheno <- names(phenoColumns)[phenoI]
    phenotypes <- phenoColumns[[phenoI]]
    categoryColumn <- categoryColumns[phenoI]
    
    for (ageI in ages) {

      phenoName <- phenotypes[ageI]

      for (sex in c("Girl", "Boy", "All")) {
        
        is <- methodValues$child_core == 1 & methodValues$pregnancy_duration_37w == 1 & methodValues[[categoryColumn]] == "longitudinal"
        
        if (sex == "Girl") {
          
          is <- is & methodValues$sex == 2
          
        } else if (sex == "Boy") {
          
          is <- is & methodValues$sex == 1
          
        }
        
        quantiles <- quantile(methodValues[is, phenoName], c(anchorDown, anchorCenter, anchorUp), na.rm = T, names = F)
        
        low <- quantiles[2] - margin99 * (quantiles[2] - quantiles[1])
        median <- quantiles[2]
        high <- quantiles[2] + margin99 * (quantiles[3] - quantiles[2])
        
        quantilesDF$quantile[k] <- low
        quantilesDF$quantileType[k] <- "low"
        quantilesDF$age[k] <- ageI
        quantilesDF$sex[k] <- sex
        quantilesDF$pheno[k] <- pheno
        k <- k + 1
        
        quantilesDF$quantile[k] <- median
        quantilesDF$quantileType[k] <- "median"
        quantilesDF$age[k] <- ageI
        quantilesDF$sex[k] <- sex
        quantilesDF$pheno[k] <- pheno
        k <- k + 1
        
        quantilesDF$quantile[k] <- high
        quantilesDF$quantileType[k] <- "high"
        quantilesDF$age[k] <- ageI
        quantilesDF$sex[k] <- sex
        quantilesDF$pheno[k] <- pheno
        k <- k + 1
        
      }
    }
  }
  
  return(quantilesDF)
  
}


#' Returns the low, median, and high length and weight at the different time points in a data frame.
#' 
#' @param methodValues the values to inspect
#' 
#' @return returns the population length weight summary statistics
getPopulationLengthWeight <- function(methodValues) {
  
  populationLengthWeight <- getPopulationQuantiles(
    methodValues = methodValues,
    phenoColumns = list(
      length = length_columns, 
      weight = weight_columns
    ),
    categoryColumns = c("lengthLongitudinalCategory", "weightLongitudinalCategory"),
    ages = c(0:11)
  )
  
  return(populationLengthWeight)
  
}


#' Returns the low, median, and high growth and weight increase at the different time points in a data frame.
#' 
#' @param methodValues the values to inspect
#' 
#' @return returns the population growth and weight increase summary statistics
getPopulationGrowth <- function(methodValues) {
  
  populationGrowth <- getPopulationQuantiles(
    methodValues = methodValues,
    phenoColumns = list(
      growth = growth_columns, 
      weight_gain = weight_gain_columns
        ),
    categoryColumns = c("lengthLongitudinalCategory", "weightLongitudinalCategory"),
    ages = c(0:10)
  )
  
  return(populationGrowth)
  
}


#' Returns the population value from the given data frame.
#' 
#' @param populationDF the population values as produced by getPopulationQuantiles
#' @param quantileType the type of quantile (low, median, high)
#' @param ageI the age index
#' @param sex the sex if available
#' @param pheno the pheno
#' 
#' @return returns the population value
getPopulationValue <- function(
  populationDF, 
  quantileType, 
  ageI, 
  sex, 
  pheno
) {
  
  if (sex != "Girl" && sex != "Boy") {
    sex <- "All"
  }
  
  result <- populationDF$quantile[
    populationDF$quantileType == quantileType &
      populationDF$age == ageI &
      populationDF$sex == sex &
      populationDF$pheno == pheno
    ]
  
  if (length(result) != 1) {
    stop(paste0("No or multiple population value found.\nageI: ", ageI, "\nsex: ", sex, "\npheno: ", pheno, "\nquantileType: ", quantileType))
  }
  
  return(result[1])
  
}


#' Returns the population curve from the given data frame.
#' 
#' @param populationDF the population values as produced by getPopulationQuantiles
#' @param quantileType the type of quantile (low, median, high)
#' @param ageIs the age indexes
#' @param sex the sex or All for all
#' @param pheno the pheno
#' 
#' @return returns the population curve
getPopulationCurve <- function(
  populationDF, 
  quantileType, 
  ageIs, 
  sex, 
  pheno
) {
  
  result <- numeric(length(ageIs))
  
  for (i in 1:length(ageIs)) {
    
    ageI <- ageIs[i]
    result[i] <- getPopulationValue(
      populationDF = populationDF, 
      quantileType = quantileType, 
      ageI = ageI, 
      sex = sex, 
      pheno = pheno
    )
    
  }
  
  return(result)
  
}


#' Updates the growth and weight gain values.
#' 
#' @param methodValues the values to inspect
#' 
#' @return returns the updated values
getGrowth <- function(methodValues) {
  
  # Get growth
  for (ageI in 1:(length(length_columns) - 1)) {
    
    # Length
    
    colName <- growth_columns[ageI]
    phenoNameAtI <- length_columns[ageI]
    phenoNameAtNextI <- length_columns[ageI + 1]
    
    if (!phenoNameAtI %in% names(methodValues)) {
      stop(paste0("Column ", phenoNameAtI, " not found in method values."))
    }
    if (!phenoNameAtNextI %in% names(methodValues)) {
      stop(paste0("Column ", phenoNameAtNextI, " not found in method values."))
    }
    if (sum(!is.na(methodValues[[phenoNameAtI]]) & methodValues[[phenoNameAtI]] <= 0)) {
      stop(paste0("Null or negative length."))
    }
    if (sum(!is.na(methodValues[[phenoNameAtNextI]]) & methodValues[[phenoNameAtNextI]] <= 0)) {
      stop(paste0("Null or negative length."))
    }
    
    methodValues[[colName]] <- ifelse(
      !is.na(methodValues[[phenoNameAtI]]) & !is.na(methodValues[[phenoNameAtNextI]]),
      log2(methodValues[[phenoNameAtNextI]] / methodValues[[phenoNameAtI]]),
      NA
    )
    
    
    # Weight
    
    colName <- weight_gain_columns[ageI]
    phenoNameAtI <- weight_columns[ageI]
    phenoNameAtNextI <- weight_columns[ageI + 1]
    
    if (!phenoNameAtI %in% names(methodValues)) {
      stop(paste0("Column ", phenoNameAtI, " not found in method values."))
    }
    if (!phenoNameAtNextI %in% names(methodValues)) {
      stop(paste0("Column ", phenoNameAtNextI, " not found in method values."))
    }
    if (sum(!is.na(methodValues[[phenoNameAtI]]) & methodValues[[phenoNameAtI]] <= 0)) {
      stop(paste0("Null or negative length."))
    }
    if (sum(!is.na(methodValues[[phenoNameAtNextI]]) & methodValues[[phenoNameAtNextI]] <= 0)) {
      stop(paste0("Null or negative length."))
    }
    
    methodValues[[colName]] <- ifelse(
      !is.na(methodValues[[phenoNameAtI]]) & !is.na(methodValues[[phenoNameAtNextI]]),
      log2(methodValues[[phenoNameAtNextI]] / methodValues[[phenoNameAtI]]),
      NA
    )
    
  }
  
  return(methodValues)
  
}


#' Updates the growth and weight gain values.
#' 
#' @param methodValues the values to inspect
#' @param populationGrowth the reference population growth metrics
#' 
#' @return returns the updated values
getNormalizedGrowth <- function(
  methodValues, 
  populationGrowth
) {
  
  phenotypes <- list(
    growth = growth_columns,
    weightIncrease = weight_gain_columns
  )
    
    for (phenoI in 1:length(phenotypes)) {
      
      pheno = names(phenotypes)[phenoI]
      phenoNames <- phenotypes[[phenoI]]
  
  for (ageI in 1:length(phenoNames)) {
      
      phenoName <- phenoNames[ageI]
      normalizedPhenoName <- paste0(phenoName, "_normalized")
      
      methodValues[[normalizedPhenoName]] <- NA
      
      for (sex in unique(methodValues$sex)) {
        
        lowValue <- getPopulationValue(
          populationDF = populationGrowth, 
          quantileType = "low", 
          ageI = ageI, 
          sex = sex, 
          pheno = pheno
        )
        medianValue <- getPopulationValue(
          populationDF = populationGrowth, 
          quantileType = "median", 
          ageI = ageI, 
          sex = sex, 
          pheno = pheno
        )
        highValue <- getPopulationValue(
          populationDF = populationGrowth, 
          quantileType = "high", 
          ageI = ageI, 
          sex = sex, 
          pheno = pheno
        )
        
        rows <- !is.na(methodValues[[phenoName]]) & methodValues$sex == sex & methodValues[[phenoName]] >= medianValue
        methodValues[rows, normalizedPhenoName] <- (methodValues[rows, phenoName] - medianValue)/(highValue - medianValue)
        
        rows <- !is.na(methodValues[[phenoName]]) & methodValues$sex == sex & methodValues[[phenoName]] < medianValue
        methodValues[rows, normalizedPhenoName] <- (methodValues[rows, phenoName] - medianValue)/(medianValue - lowValue)
        
      }
    }
    }
  
  test <- data.frame(
    a = c(1, 2, 3),
    b = c(2, NA, 4),
    c = c(3, 4, 5),
    d = c("a", "b", "c"),
    stringsAsFactors = F
  )
  
  test$median <- apply(
    X = test[, c("a", "b", "c")], 
    MARGIN = 1,
    FUN = median,
    na.rm = T
    )
  
  # Normalize each growth curve
  normalized_growth_columns <- paste0(growth_columns, "_normalized")
  methodValues$normalized_growth_median = apply(
    X = methodValues[, normalized_growth_columns], 
    MARGIN = 1,
    FUN = median,
    na.rm = T
  )
  
  centered_normalized_growth_columns <- paste0(growth_columns, "_normalized_centered")
  
  for (ageI in 1:length(normalized_growth_columns)) {
    
    methodValues[[centered_normalized_growth_columns[ageI]]] <- methodValues[[normalizedGrowth]] - methodValues$normalized_growth_median
    
  }
  
  normalized_weight_gain_columns <- paste0(weight_gain_columns, "_normalized")
  methodValues$normalized_weight_gain_median = apply(
    X = methodValues[, normalized_weight_gain_columns], 
    MARGIN = 1,
    FUN = median,
    na.rm = T
  )
  
  centered_normalized_weight_gain_columns <- paste0(weight_gain_columns, "_normalized_centered")
  
  for (ageI in 1:length(normalized_weight_gain_columns)) {
    
    methodValues[[centered_normalized_weight_gain_columns[ageI]]] <- methodValues[[normalized_weight_gain_columns]] - methodValues$normalized_weight_gain_median
    
  }
  
  return(methodValues)
  
}


#' Removes growth outliers from the given data frame and returns it.
#' 
#' @param methodValues the values to process
#' 
#' @return the updated data frame
cleanGrowthOutliers <- function(
  methodValues
) {
  
  for (index in 1:length(length_columns)) {
    
    colName <- length_columns[index]
    oldValues <- methodValues[[colName]]
    
    centered_normalized_growth_columns <- paste0(growth_columns, "_normalized_centered")
    
    outliers = apply(
      X = methodValues[, centered_normalized_growth_columns], 
      MARGIN = 1,
      FUN = isOutlierGrowth,
      index = index
    )
    
    if (sum(outliers) > 0) {
      
      methodValues[outliers, colName] <- NA
      
      newValues <- methodValues[[colName]]
      
      changedIs <- is.na(newValues) & !is.na(oldValues)
      
      if (sum(changedIs) > 0) {
        
        log <- paste0("outlier@", colName)
        methodValues[changedIs, "log"] <- ifelse(methodValues[changedIs, "log"] == "", log, paste(methodValues[changedIs, "log"], log, sep = " "))
        
      }
    }
  }
  
  return(methodValues)
}


#' Removes weight increase outliers from the given data frame and returns it.
#' 
#' @param methodValues the values to process
#' 
#' @return the updated data frame
cleanWeightIncreaseOutliers <- function(
  methodValues
) {
  
  for (index in 1:length(length_columns)) {
    
    colName <- weight_columns[index]
    oldValues <- methodValues[[colName]]
    
    centered_normalized_weight_gain_columns <- paste0(weight_gain_columns, "_normalized_centered")
    
    outliers = apply(
      X = methodValues[, centered_normalized_weight_gain_columns], 
      MARGIN = 1,
      FUN = isOutlierGrowth,
      index = index
    )
    
    if (sum(outliers) > 0) {
      
      methodValues[outliers, colName] <- NA
      
      newValues <- methodValues[[colName]]
      
      changedIs <- is.na(newValues) & !is.na(oldValues)
      
      if (sum(changedIs) > 0) {
        
        log <- paste0("outlier@", colName)
        methodValues[changedIs, "log"] <- ifelse(methodValues[changedIs, "log"] == "", log, paste(methodValues[changedIs, "log"], log, sep = " "))
        
      }
    }
  }
  
  return(methodValues)
}


#' Interpolates and smoothes the growth curves along a reference growth curve so that the growth is always positive.
#' 
#' @param methodValues the values to process
#' @param growthCurves the growth curves to use as reference
#' @param populationAnchors reference values for the population length
#' 
#' @return the corrected data frame
correctNegativeGrowth <- function(
  methodValues, 
  growthCurves, 
  populationAnchors
) {
  
  for (i in 1:nrow(methodValues)) {
    
    lengthValues <- as.numeric(methodValues[i, length_columns])
    
    newLengthValues <- correctGrowthCurve(
      child_id = methodValues[i, "child_id"], 
      sex = methodValues[i, "sex"], 
      lengthValues = lengthValues, 
      growthCurves = growthCurves, 
      populationAnchors = populationAnchors
    )
    
    if (!is.null(newLengthValues)) {
      
      for (j in 1:length(length_columns)) {
        
        methodValues[i, length_columns[j]] <- newLengthValues[j]
        
      }
      
      changedIs <- which(lengthValues != newLengthValues)
      log <- paste0("CorrectedGrowth@", paste(changedIs, collapse = ""))
      methodValues[i, "log"] <- ifelse(methodValues[i, "log"] == "", log, paste(methodValues[i, "log"], log, sep = " "))
      
      weightValues <- as.numeric(methodValues[i, weight_columns])
      
      newWeightValues <- correctWeightIncreaseCurve(
        oldLengthValues = lengthValues, 
        newLengthValues = newLengthValues,
        weightValues = weightValues,
        sex = methodValues[i, "sex"], 
        growthCurves = growthCurves
      )
      
      changedIs <- which(weightValues != newWeightValues)
      
      if (length(changedIs) > 0) {
        
        log <- paste0("CorrectedWeightIncrease@", paste(changedIs, collapse = ""))
        methodValues[i, "log"] <- ifelse(methodValues[i, "log"] == "", log, paste(methodValues[i, "log"], log, sep = " "))
        
        for (j in 1:length(weight_columns)) {
          
          methodValues[i, weight_columns[j]] <- newWeightValues[j]
          
        }
      }
    }
  }
  
  return(methodValues)
  
}


#' Imputes the length missing values in the given data frame and returns an updated version.
#' 
#' @param methodValues the values to process
#' @param growthCurves the growth curves to use as reference
#' 
#' @return returns the updated values
imputeLengthMissingValues <- function(
  methodValues, 
  growthCurves
) {
  
  for (ageI in minValidValuesBeforeFirstImputed:(11-minValidValuesAfterLastImputed)) {
    
    colName <- length_columns[ageI]
    
    oldValues <- methodValues[[colName]]
    
    outliers = apply(
      X = methodValues[, centered_normalized_growth_columns], 
      MARGIN = 1,
      FUN = isOutlierGrowth,
      index = index
    )
    
    methodValues %>%
      mutate(
        !!colName := 
          as.numeric(
            pmap(
              list(
                sex,
                length_birth,
                length_6w,
                length_3m,
                length_6m,
                length_8m,
                length_1y,
                length_16m,
                length_2y,
                length_3y,
                length_5y,
                length_7y,
                length_8y,
                lengthLongitudinalCategory
              ),
              ~getImputationFromCurveBySex(
                index = ageI + 1, 
                measuredValues = c(..2, ..3, ..4, ..5, ..6, ..7, ..8, ..9, ..10, ..11, ..12, ..13), 
                populationGrowth = growthCurves,
                sex = ..1,
                pheno = "growth",
                longitudinalCategory = ..14,
                nNeighbors = 3
              )
            )
          )
      ) -> methodValues
    
    newValues <- methodValues[[colName]]
    
    changedIs <- !is.na(newValues) & is.na(oldValues)
    log <- paste0("imputed@", colName)
    methodValues[changedIs, "log"] <- ifelse(methodValues[changedIs, "log"] == "", log, paste(methodValues[changedIs, "log"], log, sep = " "))
    
  }
  
  return(methodValues)
  
}


#' Imputes the weight missing values in the given data frame and returns an updated version.
#' 
#' @param methodValues the values to process
#' @param growthCurves the growth curves to use as reference
#' 
#' @return returns the updated values
imputeWeightMissingValues <- function(
  methodValues, 
  growthCurves
) {
  
  for (ageI in minValidValuesBeforeFirstImputed:(11-minValidValuesAfterLastImputed)) {
    
    colName <- weight_columns[ageI]
    
    oldValues <- methodValues[[colName]]
    
    methodValues %>%
      mutate(
        !!colName := 
          as.numeric(
            pmap(
              list(
                sex,
                weight_birth,
                weight_6w,
                weight_3m,
                weight_6m,
                weight_8m,
                weight_1y,
                weight_16m,
                weight_2y,
                weight_3y,
                weight_5y,
                weight_7y,
                weight_8y,
                weightLongitudinalCategory
              ),
              ~getImputationFromCurveBySex(
                index = ageI + 1, 
                measuredValues = c(..2, ..3, ..4, ..5, ..6, ..7, ..8, ..9, ..10, ..11, ..12, ..13), 
                populationGrowth = growthCurves,
                sex = ..1,
                pheno = "weightIncrease",
                longitudinalCategory = ..14,
                nNeighbors = 3
              )
            )
          )
      ) -> methodValues
    
    newValues <- methodValues[[colName]]
    
    changedIs <- !is.na(newValues) & is.na(oldValues)
    log <- paste0("imputed@", colName)
    methodValues[changedIs, "log"] <- ifelse(methodValues[changedIs, "log"] == "", log, paste(methodValues[changedIs, "log"], log, sep = " "))
    
  }
  
  return(methodValues)
  
}


#' Sets the height and weight longitudinal categories.
#' 
#' @param methodValues the values to process
#' 
#' @return returns the updated values
setLongitudinalCategories <- function(methodValues) {
  
  columnsBefore2 <- c(
    "length_birth",
    "length_6w",
    "length_3m",
    "length_6m",
    "length_8m",
    "length_1y",
    "length_16m"
  )
  columnsAfter2 <- c(
    "length_2y",
    "length_3y",
    "length_5y",
    "length_7y",
    "length_8y",
    "length_14y"
  )
  methodValues$lengthLongitudinalCategory <- "sparse"
  is <- rowSums(!is.na(methodValues[, columnsBefore2])) >= minValidValuesBefore2
  methodValues$lengthLongitudinalCategory[is] <- "early"
  is <- rowSums(!is.na(methodValues[, columnsBefore2])) >= minValidValuesBefore2 &
    rowSums(!is.na(methodValues[, columnsAfter2])) >= minValidValuesAfter2
  methodValues$lengthLongitudinalCategory[is] <- "longitudinal"
  
  columnsBefore2 <- c(
    "weight_birth",
    "weight_6w",
    "weight_3m",
    "weight_6m",
    "weight_8m",
    "weight_1y",
    "weight_16m"
  )
  columnsAfter2 <- c(
    "weight_2y",
    "weight_3y",
    "weight_5y",
    "weight_7y",
    "weight_8y",
    "weight_14y"
  )
  methodValues$weightLongitudinalCategory <- "sparse"
  methodValues$weightLongitudinalCategory <- "sparse"
  is <- rowSums(!is.na(methodValues[, columnsBefore2])) >= minValidValuesBefore2
  methodValues$weightLongitudinalCategory[is] <- "early"
  is <- rowSums(!is.na(methodValues[, columnsBefore2])) >= minValidValuesBefore2 &
    rowSums(!is.na(methodValues[, columnsAfter2])) >= minValidValuesAfter2
  methodValues$weightLongitudinalCategory[is] <- "longitudinal"
  
  return(methodValues)
  
}


#' Iterates through the cleaning steps until no value is changed.
#' 
#' @param values the values to process
#' 
#' @return returns the updated values
iterativeCleaning <- function(values) {
  
  if (any(is.na(values$sex))) {
    
    stop("Missing value for sex")
    
  }
  
  diff <- T
  iterationIndex <- 1
  inspectNegativeGrowth <- F
  
  while (diff) {
    
    previousValues <- values
    
    # Make categories based on the number of missing values
    
    print(paste(Sys.time(), iterationIndex, " Inspecting number of missing values"))
    
    values <- setLongitudinalCategories(values)
    
    print(paste("Iteration ", iterationIndex))
    print("Length")
    print(table(values$lengthLongitudinalCategory))
    print("Length Term Unrelated")
    print(table(values$lengthLongitudinalCategory[values$unrelated == 1 & values$pregnancy_duration_over_37w == 1]))
    print("Weight")
    print(table(values$weightLongitudinalCategory))
    print("Weight Term Unrelated")
    print(table(values$weightLongitudinalCategory[values$unrelated == 1 & values$pregnancy_duration_over_37w == 1]))
    
    
    # get growth curves
    
    print(paste(Sys.time(), iterationIndex, " Getting growth curves"))
    
    values <- getGrowth(values)
    refValues <- values[values$unrelated == 1 & values$pregnancy_duration_over_37w == 1, ]
    populationGrowth <- getPopulationGrowth(refValues)
    values <- getNormalizedGrowth(values, populationGrowth)
    
    
    # Inspect growth outliers
    
    print(paste(Sys.time(), iterationIndex, " Removing growth outliers"))
    
    values <- cleanGrowthOutliers(
      methodValues = values
    )
    
    
    # Inspect weight increase outliers
    
    print(paste(Sys.time(), iterationIndex, " Removing weight increase outliers"))
    
    values <- cleanWeightIncreaseOutliers(
      methodValues = values
    )
    
    
    # Update categories based on the number of missing values
    
    print(paste(Sys.time(), iterationIndex, " Inspecting number of missing values"))
    
    values <- setLongitudinalCategories(values)
    
    
    # Update growth curves
    
    print(paste(Sys.time(), iterationIndex, " Updating growth curves"))
    
    values <- getGrowth(values)
    
    refValues <- values[values$unrelated == 1 & values$pregnancy_duration_over_37w == 1, ]
    populationGrowth <- getPopulationGrowth(refValues)
    
    values <- getNormalizedGrowth(values, populationGrowth)
    
    
    # Impute missing values
    
    print(paste(Sys.time(), iterationIndex, " Imputation of missing values"))
    
    values <- imputeLengthMissingValues(
      methodValues = values, 
      growthCurves = populationGrowth
    )
    values <- imputeWeightMissingValues(
      methodValues = values, 
      growthCurves = populationGrowth
    )
    
    if (!inspectNegativeGrowth) {
      
      # Compare the new values to the previous iteration
      
      diff <- F
      
      for (ageI in 0:11) {
        
        for (pheno in c("length", "weight")) {
          
          colName <- paste0(pheno, ageI)
          
          if (sum(
            is.na(previousValues[[colName]]) & !is.na(values[[colName]]) |
            !is.na(previousValues[[colName]]) & is.na(values[[colName]]) |
            !is.na(previousValues[[colName]]) & !is.na(values[[colName]]) &
            previousValues[[colName]] != values[[colName]]) > 0) {
            
            diff <- T
            break
            
          }
        }
        if (diff) {
          break
        }
      }
      
      if (!diff) {
        
        inspectNegativeGrowth <- T
        
      }
    }
    
    if (inspectNegativeGrowth) {
      
      
      # Update categories based on the number of missing values
      
      print(paste(Sys.time(), iterationIndex, " Inspecting number of missing values"))
      
      values <- setLongitudinalCategories(values)
      
      
      # Update growth curves
      
      print(paste(Sys.time(), iterationIndex, " Updating growth curves"))
      
      values <- getGrowth(values)
      
      refValues <- values[values$unrelated == 1 & values$pregnancy_duration_over_37w == 1, ]
      populationGrowth <- getPopulationGrowth(refValues)
      
      values <- getNormalizedGrowth(values, populationGrowth)
      
      
      # Get population summary statistics
      
      print(paste(Sys.time(), iterationIndex, " Getting population summary statistics"))
      
      refValues <- values[values$unrelated == 1 & values$pregnancy_duration_over_37w == 1, ]
      populationLengthWeight <- getPopulationLengthWeight(refValues)
      
      
      # Correct for negative growth
      
      print(paste(Sys.time(), iterationIndex, " Correcting for negative growth"))
      
      values <- correctNegativeGrowth(
        methodValues = values, 
        growthCurves = populationGrowth, 
        populationAnchors = populationLengthWeight
      )
      
    }
    
    iterationIndex <- iterationIndex + 1
    
    if (iterationIndex > maxIt) {
      
      print("Maximal number of iterations reached.")
      break
      
    }
    
    # Compare the new values to the previous iteration
    
    diff <- F
    
    for (colName in c(weight_columns, length_columns)) {
        
        colName <- paste0(pheno, ageI)
        
        if (sum(
          is.na(previousValues[[colName]]) & !is.na(values[[colName]]) |
          !is.na(previousValues[[colName]]) & is.na(values[[colName]]) |
          !is.na(previousValues[[colName]]) & !is.na(values[[colName]]) &
          previousValues[[colName]] != values[[colName]]) > 0) {
          
          diff <- T
          break
          
        }
      }
  }
  
  return(values)
  
}


#' Iterates through the cleaning steps until no value is changed.
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


#' Iterates through the cleaning steps until no value is changed.
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
      original_length = originalLength,
      new_length = newLength
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
      geom_line(data = lengthDF, mapping = aes(x = x, y = original_length, group = 1), col = "grey40", linetype = "dotted") +
      geom_line(data = lengthDF, mapping = aes(x = x, y = new_length, group = 1), col = "darkblue", linetype = "solid") +
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
      original_weight = originalWeight,
      new_weight = newWeight
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
      geom_line(data = weightDF, mapping = aes(x = x, y = original_weight, group = 1), col = "grey40", linetype = "dotted") +
      geom_line(data = weightDF, mapping = aes(x = x, y = new_weight, group = 1), col = "blue3", linetype = "solid") +
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

