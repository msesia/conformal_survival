library(tidyverse)

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

predict_CQR <- function(data.test, surv_model, data.cal, alpha) {
    # Which quantile to predict? 1-alpha is a reasonable choice,
    # but in theory any other value can be used
    probs <- c(1-alpha)

    ## Calibration
    # Predict the survival quantiles for the given nominal percentile
    pred.cal <- surv_model$predict_quantiles(data.cal, probs = probs)
    pred.cal <- as.numeric(pred.cal[[1]])

    # Replace NAs with large numbers
    Y.cal <- data.cal$time
    pred.cal[is.na(pred.cal)] <- max(Y.cal)
    # Compute the conformity scores for CQR
    scores <- pred.cal - Y.cal
    # Calculate the (1-alpha) quantile of the scores (with finite-sample correction)
    n <- length(scores)
    percentile_scores <- (1 - alpha) * (1+1/n)
    calibration <- as.numeric(quantile(scores, probs = percentile_scores, type = 1))

    ## Prediction for test data

    # Predict the survival quantiles for the given nominal percentile
    pred.test <- surv_model$predict_quantiles(data.test, probs = probs)
    pred.test <- as.numeric(pred.test[[1]])
    # Compute the calibrated lower bounds (don't go below 0)
    lower.test <- pmax(0, pred.test - calibration)

    return(list(uncalibrated=pred.test, calibrated=lower.test))
}

predict_Candes <- function(data.test, surv_model, data.cal, C.cal, alpha, c0=NULL) {
    ## For special case with constant weights

    # Which quantile to predict? 1-alpha is a reasonable choice,
    # but in theory any other value can be used
    probs <- c(1-alpha)

    ## Choose c0 if not supplied
    if(is.null(c0)) {
        c0 <- median(C.cal)
    }

    ## Remove calibration points with small censoring
    idx.keep <- which(C.cal >= c0)

    ## Calibration
    Y.cal <- data.cal$time

    # Predict the survival quantiles for the given nominal percentile
    pred.cal <- surv_model$predict_quantiles(data.cal[idx.keep,], probs = probs)
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
    pred.test <- surv_model$predict_quantiles(data.test, probs = probs)
    pred.test <- as.numeric(pred.test[[1]])
    # Compute the calibrated lower bounds (don't go below 0)
    lower.test <- pmin(c0, pmax(0, pred.test - calibration))

    return(list(uncalibrated=pred.test, calibrated=lower.test))
}

predict_prototype <- function(data.test, surv_model, imp_model, data.cal, alpha, c0=NULL) {
    ## For special case with constant weights

    # Initialize the censoring times equal to the observed times
    C.cal <- data.cal$time

    # Extract the covariate information (remove 'time' and 'status' columns, if present)
    X.cal <- data.cal |> select(-any_of(c("time", "status")))
    Y.cal <- data.cal$time

    # Find the indices of observations for which the censoring time needs to be imputed
    idx.event <- which(data.cal$status==TRUE)
    
    if(length(idx.event) > 0) {
        # Impute the missing censoring times using the oracle (data-generating) model
        C.cal[idx.event] <- imp_model$sample_censoring_times(X.cal[idx.event,], T=Y.cal[idx.event])
    }

    # Apply Candes' method using the imputed censoring times
    out <- predict_Candes(data.test, surv_model, data.cal, C.cal, alpha, c0=c0)
    return(out)
}
