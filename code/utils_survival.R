suppressMessages(library(tidyverse))

## Function to compute survival quantiles at specified percentiles
compute_survival_quantiles <- function(fit, newX=NULL,
                                       percentiles = c(0.5, 0.75, 0.9)) {
    ## Predict survival curves
    predictions <- predict(fit, newdata=newX)
    survival_curves <- predictions$predictions
    time_points <- fit$failure.times

    ## Function to find the time at which the survival probability crosses a specific percentile
    find_quantile <- function(survival_probs, time_points, percentile) {
        index <- which.max(survival_probs <= percentile)
        if (index == 1) {
            return(NA)  ## No time found
        } else {
            return(time_points[index])
        }
    }

    ## Initialize a list to store quantiles for each individual
    quantiles_list <- list()

    ## Loop over each individual
    for (i in 1:nrow(survival_curves)) {
        quantiles <- sapply(percentiles, function(p) find_quantile(survival_curves[i, ], time_points, p))
        quantiles_list[[i]] <- quantiles
    }

    ## Convert the list to a data frame
    quantiles_df <- do.call(rbind, quantiles_list)
    colnames(quantiles_df) <- paste0("Q", percentiles * 100, "%")
    rownames(quantiles_df) <- paste0("Individual_", 1:nrow(quantiles_df))

    return(as.data.frame(quantiles_df))
}

predict_CQR <- function(X.test, fit, X.cal, Y.cal, alpha) {
    # Which quantile to predict? 1-alpha is a reasonable choice,
    # but in theory any other value can be used
    nominal_percentile <- c(1-alpha)

    ## Calibration
    # Predict the survival quantiles for the given nominal percentile
    pred.cal <- compute_survival_quantiles(fit, newX=X.cal,
                                           percentiles=nominal_percentile)
    pred.cal <- as.numeric(pred.cal[[1]])
    # Replace NAs with large numbers
    pred.cal[is.na(pred.cal)] <- max(Y.cal)
    # Compute the conformity scores for CQR
    scores <- pred.cal - Y.cal
    # Calculate the (1-alpha) quantile of the scores (with finite-sample correction)
    n <- length(scores)
    percentile_scores <- (1 - alpha) * (1+1/n)
    calibration <- as.numeric(quantile(scores, probs = percentile_scores, type = 1))

    ## Prediction for test data

    # Predict the survival quantiles for the given nominal percentile
    pred.test <- compute_survival_quantiles(fit, newX=X.test,
                                            percentiles=nominal_percentile)
    pred.test <- as.numeric(pred.test[[1]])
    # Compute the calibrated lower bounds (don't go below 0)
    lower.test <- pmax(0, pred.test - calibration)

    return(list(uncalibrated=pred.test, calibrated=lower.test))
}

predict_Candes <- function(X.test, fit, X.cal, Y.cal, C.cal, alpha, c0=NULL) {
    ## For special case with constant weights

    # Which quantile to predict? 1-alpha is a reasonable choice,
    # but in theory any other value can be used
    nominal_percentile <- c(1-alpha)

    ## Choose c0 if not supplied
    if(is.null(c0)) {
        c0 <- median(C.cal)
    }

    ## Remove calibration points with small censoring
    idx.keep <- which(C.cal >= c0)
    
    ## Calibration
    # Predict the survival quantiles for the given nominal percentile
    pred.cal <- compute_survival_quantiles(fit, newX=X.cal[idx.keep,],
                                           percentiles=nominal_percentile)
    pred.cal <- as.numeric(pred.cal[[1]])
    # Replace NAs with large numbers
    pred.cal[is.na(pred.cal)] <- max(Y.cal[idx.keep])
    # Compute the conformity scores for CQR
    scores <- pred.cal - pmin(Y.cal[idx.keep], c0)
    # Calculate the (1-alpha) quantile of the scores (with finite-sample correction)
    n <- length(scores)
    percentile_scores <- (1 - alpha) * (1+1/n)
    calibration <- as.numeric(quantile(scores, probs = percentile_scores, type = 1))

    ## Prediction for test data

    # Predict the survival quantiles for the given nominal percentile
    pred.test <- compute_survival_quantiles(fit, newX=X.test,
                                            percentiles=nominal_percentile)
    pred.test <- as.numeric(pred.test[[1]])
    # Compute the calibrated lower bounds (don't go below 0)
    lower.test <- pmin(c0, pmax(0, pred.test - calibration))

    return(list(uncalibrated=pred.test, calibrated=lower.test))
}


evaluate_bounds <- function(observed_time, lower_bound, event_time=NULL) {
    coverage_observed <- mean(lower_bound <= observed_time)
    if(is.null(event_time)) {
        coverage_event_time <- NA
    } else {
        coverage_event_time <- mean(lower_bound <= event_time)
    }
    mean_lower_bound <- mean(lower_bound)
    out <- tibble("Coverage (observed time)"=coverage_observed,
                  "Mean lower bound"=mean_lower_bound,
                  "Coverage (event time)"=coverage_event_time)
    return(out)
}
