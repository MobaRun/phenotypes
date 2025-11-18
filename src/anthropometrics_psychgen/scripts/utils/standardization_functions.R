#' Imputes the values in x using the nearest neighbors in y.
#' 
#' @param x the values to impute.
#' @param y the values to use as nearest neighbors.
#' 
#' @return x with imputed values.
imputeXValue <- function(
    x,
    y,
    nNeighbors = 100
) {
    
    imputationDF <- data.frame(
        refX = x,
        refY = y,
        stringsAsFactors = F
    ) %>%
        filter(
            !is.na(refX) & !is.na(refY) & !is.infinite(refX) & !is.infinite(refY)
        ) %>%
        arrange(
            refY
        )
    
    for (i in 1:length(x)) {
        
        if (is.na(x[i]) || is.infinite(x[i])) {
            
            if (!is.na(y[i]) && !is.infinite(y[i])) {
                
                diff <- abs(imputationDF$refY - y[i])
                minDiff <- min(diff)
                minDiffI <- which(diff == minDiff)
                targetIMedian <- median(minDiffI)
                
                targetIDown <- targetIMedian - nNeighbors / 2
                targetIUp <- targetIMedian + nNeighbors / 2
                
                j <- 1:nrow(imputationDF)
                targetI <- j[j >= targetIDown & j <= targetIUp]
                
                x[i] <- median(imputationDF$refX[targetI])
                
            }
        }
        
    }
    
    return(x)
    
}


#' Standardizes the phenotypes separately for male and female children. Infinite values are replaced by NA.
#' Warning: gamlss uses data from the global environment, the following objects will be silently overwritten: 'maleTrainingDF' and 'femaleTrainingDF'.
#' 
#' @param trainingvalues the data frame containing the phenotypes to use as training set
#' @param values the association data frame containing the phenotypes to use for prediction
#' @param id the name of the column containing the identifiers
#' @param x the name of the variable to consider in x, set to constant if not in the data frame
#' @param imputeX a boolean indicating whether missing values in x should be imputed based on the nearest neighbors in y
#' @param y the name of the variable to consider in y
#' @param zY the name of the standardized y variable
#' @param formula the formula to use for the model
#' @param sigmaFormula the formula to use for sigma
#' @param formula2 the formula to use for the model if the first one produced NA for estimates
#' @param sigmaFormula2 the formula to use for sigma if the first one produced NA for estimates
#' @param family the family of distribution to fit
#' 
#' @return A data frame containing child sentrix id and standardized phenotypes.
standardizeBySex <- function(
    trainingvalues,
    values,
    id,
    x,
    imputeX = F,
    y,
    zY,
    formula,
    sigmaFormula,
    coefficients = "mu.coefficients",
    formula2 = formula,
    sigmaFormula2 = sigmaFormula,
    family
) {
    
    # Special case where there is no "x"
  
  x_found <- T
  
  if (! x %in% names(trainingvalues)) {
    
    trainingvalues[[x]] <- 0
    values[[x]] <- 0
    x_found <- F
    
  } else if (! x %in% names(values)) {
    
    stop(paste0("x variable (", x, ") found in training set but not in prediciton."))
    
  }
    
  
    # Create training data frame for males
    
    maleTrainingDF <- trainingvalues %>%
        filter(
            sex == 1
        ) %>%
        select(
            !!sym(id), !!sym(x), !!sym(y)
        ) %>%
        filter(
            !is.na(!!sym(id)) & !is.na(!!sym(x)) & !is.infinite(!!sym(x)) & !is.na(!!sym(y)) & !is.infinite(!!sym(y))
        )
    
    
    # Check that there are at least two different values
    
    nValues <- length(unique(maleTrainingDF[[y]][!is.na(maleTrainingDF[[y]])]))
    
    if (nValues < 2) {
        
        stop(glue("Less than two distinct values available for males for standardization of {y}."))
        
    }
    
    if (sum(is.infinite(maleTrainingDF[[y]])) > 0) {
        
        stop(glue("Infinite values for {y}."))
    }
    
    # Save in the global environment so that gamlss can access it
    
    assign("maleTrainingDF", maleTrainingDF, envir = .GlobalEnv)
    
    
    # Create training data frame for females and save it in the global environment
    
    femaleTrainingDF <- trainingvalues %>%
        filter(
            sex == 2
        ) %>%
        select(
            !!sym(id), !!sym(x), !!sym(y)
        ) %>%
        filter(
            !is.na(!!sym(id)) & !is.na(!!sym(x)) & !is.infinite(!!sym(x)) & !is.na(!!sym(y)) & !is.infinite(!!sym(y))
        )
    
    
    # Check that there are at least two different values
    
    nValues <- length(unique(femaleTrainingDF[[y]][!is.na(femaleTrainingDF[[y]])]))
    
    if (nValues < 2) {
        
        stop(glue("Less than two distinct values available for females for standardization of {y}."))
        
    }
    
    if (sum(is.infinite(femaleTrainingDF[[y]])) > 0) {
      
      stop(glue("Infinite values for {y}."))
    }
    
    
    # Save in the global environment so that gamlss can access it
    
    assign("femaleTrainingDF", femaleTrainingDF, envir = .GlobalEnv)
    
    
    # Create data frames for the prediction in males
    
    maleDF <- values %>%
        filter(
            sex == 1
        ) %>%
        select(
            !!sym(id), !!sym(x), !!sym(y)
        ) %>%
        filter(
            !is.na(!!sym(id)) & !is.na(!!sym(y)) & !is.infinite(!!sym(y))
        )
    
    
    # Impute missing x values if requested
    
    if (imputeX) {
    
    maleDF[[x]] <- imputeXValue(x = maleDF[[x]], y = maleDF[[y]])
    
    }
    
    maleDF <- maleDF %>%
        filter(
            !is.na(!!sym(x)) & !is.infinite(!!sym(x))
        )
    
    
    # Create data frames for the prediction in females
    
    femaleDF <- values %>%
        filter(
            sex == 2
        ) %>%
        select(
            !!sym(id), !!sym(x), !!sym(y)
        ) %>%
        filter(
            !is.na(!!sym(id)) & !is.na(!!sym(y)) & !is.infinite(!!sym(y))
        )
    
    
    # Impute missing x values if requested
    
    if (imputeX) {
    
    femaleDF[[x]] <- imputeXValue(x = femaleDF[[x]], y = femaleDF[[y]])
    
    }
    
    femaleDF <- femaleDF %>%
        filter(
            !is.na(!!sym(x)) & !is.infinite(!!sym(x))
        )
    
    
    # Train model for males
    
    maleModel <- gamlss(
        formula = formula,
        sigma.formula = sigmaFormula,
        family = family,
        data = maleTrainingDF
    )
    
    
    # Compute z-scores for males
    
    tryCatch(
      {
        
        maleDF[[zY]] <- centiles.pred(
          obj = maleModel, 
          xname = x, 
          xvalues = maleDF[[x]], 
          yval = maleDF[[y]], 
          type = "z-scores"
        )
        
      }, error = function(error) {
        
        print("Model 1 failed for males")
        print(error)
        print("Trying Model 2")
        
        maleModel <- gamlss(
          formula = formula2,
          sigma.formula = sigmaFormula2,
          family = family,
          data = maleTrainingDF
        )
        
        maleDF[[zY]] <- centiles.pred(
          obj = maleModel, 
          xname = x, 
          xvalues = maleDF[[x]], 
          yval = maleDF[[y]], 
          type = "z-scores"
        )
        
      }
    )
    
    
    # Train model for females
    
    femaleModel <- gamlss(
        formula = formula,
        sigma.formula = sigmaFormula,
        family = family,
        data = femaleTrainingDF
    )
    
    
    # Compute z-scores for females
    
    tryCatch(
      {
        
        femaleDF[[zY]] <- centiles.pred(
          obj = femaleModel, 
          xname = x, 
          xvalues = femaleDF[[x]], 
          yval = femaleDF[[y]], 
          type = "z-scores"
        )
        
      }, error = function(error) {
        
        print("Model 1 failed for females")
        print(error)
        print("Trying Model 2")
        
        femaleModel <- gamlss(
          formula = formula2,
          sigma.formula = sigmaFormula2,
          family = family,
          data = femaleTrainingDF
        )
        
        femaleDF[[zY]] <- centiles.pred(
          obj = femaleModel, 
          xname = x, 
          xvalues = femaleDF[[x]], 
          yval = femaleDF[[y]], 
          type = "z-scores"
        )
        
      }
    )
    
    
    # Make a data frame containing ids and z-scores for both males and females.
    
    zDF <- rbind(
      maleDF %>% 
        select(
          !!sym(id), !!sym(zY)
        ), 
      femaleDF %>% 
        select(
          !!sym(id), !!sym(zY)
        )
      )
    
    
    # Replace infinite with NA
    
    zDF[[zY]][is.infinite(zDF[[zY]])] <- NA
    
    
    # Make sure that at least two values are present in the output
    
    nValues <- length(unique(zDF[[zY]][!is.na(zDF[[zY]])]))
    
    if (nValues < 2) {
        
        stop(glue("Less than two distinct values remaining after standardization of {y}."))
        
    }
    
    
    # Return result
    
    return(zDF)
    
}


#' Standardizes the phenotypes. Infinite values are replaced by NA.
#' Warning: gamlss uses data from the global environment, the following objects will be silently overwritten: 'trainingDF'.
#' 
#' @param trainingvalues the data frame containing the phenotypes to use as training set
#' @param values the association data frame containing the phenotypes to use for prediction
#' @param id the name of the column containing the identifiers
#' @param x the name of the variable to consider in x, set to constant if not in the data frame
#' @param y the name of the variable to consider in y
#' @param zY the name of the standardized y variable
#' @param formula the formula to use for the model
#' @param sigmaFormula the formula to use for sigma
#' @param family the family of distribution to fit
#' 
#' @return A data frame containing child sentrix id and standardized phenotypes.
standardize <- function(
    trainingvalues,
    values,
    id,
    x,
    y,
    zY,
    formula,
    sigmaFormula,
    family
) {
    
    # Special case where there is no "x"
    
    if (! x %in% names(trainingvalues)) {
        
        trainingvalues[[x]] <- 0
        
    }
    
    if (! x %in% names(values)) {
        
        values[[x]] <- 0
        
    }
    
    # Make a training data frame
    
    trainingDF <- trainingvalues %>%
        select(
            !!sym(id), !!sym(x), !!sym(y)
        ) %>%
        filter(
            !is.na(!!sym(id)) & !is.na(!!sym(x)) & !is.na(!!sym(y)) & !is.infinite(!!sym(y))
        )
    
    
    # Check that there are at least two different values
    
    nValues <- length(unique(trainingDF[[y]][!is.na(trainingDF[[y]])]))
    
    if (nValues < 2) {
        
        stop(glue("Less than two distinct values available for standardization of {y}."))
        
    }
    
    
    # Save in the global environment so that gamlss can access it
    
    assign("trainingDF", trainingDF, envir = .GlobalEnv)
    
    
    # Model using gamlss
    
    model <- gamlss(
        formula = formula,
        sigma.formula = sigmaFormula,
        family = family,
        data = trainingDF
    )
    
    
    # Create a data frame for the results
    
    zDF <- values %>%
        select(
            !!sym(id), !!sym(x), !!sym(y)
        ) %>%
        filter(
            !is.na(!!sym(id)) & !is.na(!!sym(y)) & !is.infinite(!!sym(y))
        )
    
    
    # Impute missing x values using nearest neighbors
    
    zDF[[x]] <- imputeXValue(x = zDF[[x]], y = zDF[[y]])
    
    
    # Compute the z-scores
    
    zDF[[zY]] <- centiles.pred(
        obj = model, 
        xname = x, 
        xvalues = zDF[[x]], 
        yval = zDF[[y]], 
        type = "z-scores"
    )
    
    
    # Keep only id and z-scores for the output
    
    zDF %>% 
        select(
            !!sym(id), !!sym(zY)
        ) -> zDF
    
    
    # Replace infinite with NA
    
    zDF[[zY]][is.infinite(zDF[[zY]])] <- NA
    
    
    # Replace values higher/lower than 5 SD with NA
    
    zDF[[zY]][!is.na(zDF[[zY]]) & abs(zDF[[zY]]) > 5] <- NA
    
    
    # Make sure that at least two values are present in the output
    
    nValues <- length(unique(zDF[[zY]][!is.na(zDF[[zY]])]))
    
    if (nValues < 2) {
        
        stop(glue("Less than two distinct values remaining after standardization of {y}."))
        
    }
    
    
    # Return results
    
    return(zDF)
    
}


#' Centers and scales values distinguishing males and females.
#' 
#' @param values The association data frame containing the phenotypes to center and scale.
#' @param unrelatedValues The association data frame containing the phenotypes to use for the estimation of summary statistics. 
#' @param y the name of the variable to scale.
#' 
#' @return A data frame containing child sentrix id and standardized phenotypes.
centerAndScaleBySex <- function(
    values,
    unrelatedValues,
    y
) {
    
    
    # Check that there are at least two different values
    
    nValues <- length(unique(unrelatedValues[[y]][!is.na(unrelatedValues[[y]]) & unrelatedValues$sex == 1]))
    
    if (nValues < 2) {
        
        stop(glue("Less than two distinct values available for males for standardization of {y}."))
        
    }
    
    nValues <- length(unique(unrelatedValues[[y]][!is.na(unrelatedValues[[y]]) & unrelatedValues$sex == 2]))
    
    if (nValues < 2) {
        
        stop(glue("Less than two distinct values available for females for standardization of {y}."))
        
    }
    
    if (sum(is.infinite(unrelatedValues[[y]])) > 0) {
        
        stop(glue("Infinite values for {y}."))
    }
    
    # Get mean
    
    meanMales <- mean(unrelatedValues[[y]][unrelatedValues$sex == 1], na.rm = T)
    meanFemales <- mean(unrelatedValues[[y]][unrelatedValues$sex == 2], na.rm = T)
    
    
    # Center
    
    result <- ifelse(
        test = values$sex == 1, 
        yes = values[[y]] - meanMales, 
        no = values[[y]] - meanFemales
    )
    
    
    # Get sd
    
    sdMales <- sd(unrelatedValues[[y]][unrelatedValues$sex == 1], na.rm = T)
    sdFemales <- sd(unrelatedValues[[y]][unrelatedValues$sex == 2], na.rm = T)
    
    if (sdFemales == 0) {
        
        stop(glue("Null sd for females in phenotype {y}."))
        
    }
    if (sdMales == 0) {
        
        stop(glue("Null sd for males in phenotype {y}."))
        
    }
    
    
    # Scale
    
    result <- ifelse(
        test = values$sex == 1, 
        yes = result / sdMales, 
        no = result / sdFemales
    )
    
    # Remove infinite and outliers
    
    result[is.infinite(result) | abs(result) > 5] <- NA
    
    
    # Make sure that at least two values are present in the output
    
    nValues <- length(unique(result[!is.na(result)]))
    
    if (nValues < 2) {
        
        stop(glue("Less than two distinct values remaining after standardization of {y}."))
        
    }
    
    return(result)
    
}


#' Centers and scales values.
#' 
#' @param values The association data frame containing the phenotypes to center and scale.
#' @param unrelatedValues The association data frame containing the phenotypes to use for the estimation of summary statistics. 
#' @param y the name of the variable to scale.
#' 
#' @return A data frame containing child sentrix id and standardized phenotypes.
centerAndScale <- function(
    values,
    unrelatedValues,
    y
) {
    
    # Check that there are at least two different values
    
    nValues <- length(unique(unrelatedValues[[y]][!is.na(unrelatedValues[[y]])]))
    
    if (nValues < 2) {
        
        stop(glue("Less than two distinct values available for standardization of {y}."))
        
    }
    
    if (sum(is.infinite(unrelatedValues[[y]])) > 0) {
        
        stop(glue("Infinite values for {y}."))
    }
    
    # Get mean
    
    meanValue <- mean(unrelatedValues[[y]], na.rm = T)
    
    
    # Center
    
    result <- values[[y]] - meanValue
    
    
    # Get sd
    
    sdValue <- sd(unrelatedValues[[y]], na.rm = T)
    
    if (is.null(sdValue)) {
        
        stop(glue("sd null for phenotype {y}."))
        
    }
    
    if (sdValue == 0) {
        
        stop(glue("sd of zero for phenotype {y}."))
        
    }
    
    
    # Scale
    
    result <- result / sdValue
    
    
    # Remove infinite and outliers
    
    result[is.infinite(result) | abs(result) > 5] <- NA
    
    
    # Make sure that at least two values are present in the output
    
    nValues <- length(unique(result[!is.na(result)]))
    
    if (nValues < 2) {
        
        stop(glue("Less than two distinct values remaining after standardization of {y}."))
        
    }
    
    return(result)
    
}

