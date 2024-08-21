library(tidyverse)

evaluate_bounds <- function(observed_time, lower_bound, event_time=NULL) {
    coverage_observed <- mean(lower_bound <= observed_time)
    if(is.null(event_time)) {
        coverage_event_time <- NA
    } else {
        coverage_event_time <- mean(lower_bound <= event_time)
    }
    mean_lower_bound <- mean(lower_bound)
    median_lower_bound <- median(lower_bound)
    out <- tibble("Coverage (observed time)"=coverage_observed,
                  "Mean lower bound"=mean_lower_bound,
                  "Median lower bound"=median_lower_bound,
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


weighted_calibration <- function(scores.cal, weights.cal, weight.new, alpha) {
  # Computes a weighted quantile for weighted conformal prediction.
  #
  # Args:
  #   scores.cal: A numeric vector of calibration scores (non-conformity scores).
  #   weights.cal: A numeric vector of weights corresponding to `scores.cal`.
  #               Represents the importance of each calibration score.
  #   weight.new: A numeric value or vector representing the weights for the new data points.
  #   alpha: A numeric value between 0 and 1 indicating the significance level (1 - alpha)
  #          for the conformal prediction interval.
  #
  # Returns:
  #   The calibration term (weighted quantile) that adjusts the prediction intervals.

  # Check if the input lengths are compatible
  if(length(scores.cal) != length(weights.cal)) {
    stop("The length of scores.cal is not compatible with the length of weights.cal!")
  }

  # Ensure alpha is a valid probability
  if(!is.numeric(alpha)) stop("alpha should be a real number between 0 and 1!")
  if(alpha > 1 | alpha < 0) stop("alpha should be a real number between 0 and 1!")

  # Combine the weights from calibration and new data
  weight <- c(weights.cal, weight.new)
  weight <- weight / sum(weight)  # Normalize weights to sum to 1

  # Combine scores and append an infinity value for the new data
  score_vec <- c(scores.cal, Inf)

  # Sort scores and reorder the weights accordingly
  sort_score <- sort(score_vec)
  order_score <- order(score_vec)
  sort_w <- weight[order_score]

  # Find the minimum index where the cumulative weight reaches or exceeds 1 - alpha
  idxw <- min(which(cumsum(sort_w) >= 1 - alpha))

  # Compute the calibration term as the score at the index determined by the cumulative weights
  calib_term <- sort_score[idxw]

  return(calib_term)
}



predict_Candes <- function(data.test, surv_model, cens_model, data.cal, C.cal, alpha, c0=NULL) {
    # Which quantile to predict? 1-alpha is a reasonable choice,
    # but in theory any other value can be used
    probs <- c(1-alpha)

    ## Choose c0 if not supplied
    if(is.null(c0)) {
        c0 <- median(C.cal)
    }

    # Extract the covariate information (remove 'time' and 'status' columns, if present)
    X.cal <- data.cal |> select(-any_of(c("time", "status")))
    X.test <- data.test |> select(-any_of(c("time", "status")))
    Y.cal <- data.cal$time

    ## Remove calibration points with small censoring
    idx.keep <- which(C.cal >= c0)

    ## Compute the conformity scores for the calibration data

    ## Predict the survival quantiles for the given nominal percentile
    pred.cal <- surv_model$predict_quantiles(data.cal[idx.keep,], probs = probs)
    pred.cal <- as.numeric(pred.cal[[1]])
    # Replace NAs with large numbers
    pred.cal[is.na(pred.cal)] <- max(Y.cal[idx.keep])
    # Compute the conformity scores for CQR
    scores.cal <- pred.cal - pmin(Y.cal[idx.keep], c0)

    ## Compute conformal weights
    weights.cal <- 1/pmax(1e-6, cens_model$predict_survival(X.cal[idx.keep,], failure.times=c0)$predictions)
    weights.test <- 1/pmax(1e-6, cens_model$predict_survival(X.test, failure.times=c0)$predictions)

    ## Prediction for test data
    n <- length(scores.cal)
    n.test <- nrow(data.test)
    lower.test <- rep(-1, n.test)

    pred.test.tmp <- surv_model$predict_quantiles(data.test, probs = probs)
    pred.test <- as.numeric(pred.test.tmp[[1]])

    lower.test <- sapply(1:n.test, function(i) {
        ## Compute weighted quantile of conformity scores
        calibration <- weighted_calibration(scores.cal, weights.cal, weights.test[i], alpha)

        ## Compute the calibrated lower bounds (don't go below 0)
        lower <- pmin(c0, pmax(0, pred.test[i] - calibration))
        return(lower)
    })

    return(list(uncalibrated=pred.test, calibrated=lower.test))
}


predict_prototype <- function(data.test, surv_model, cens_imputator, data.cal, alpha, c0=NULL) {
    # Initialize the censoring times equal to the observed times
    C.cal <- data.cal$time

    # Extract the covariate information (remove 'time' and 'status' columns, if present)
    X.cal <- data.cal |> select(-any_of(c("time", "status")))
    Y.cal <- data.cal$time

    # Find the indices of observations for which the censoring time needs to be imputed
    idx.event <- which(data.cal$status==TRUE)

    if(length(idx.event) > 0) {
        # Impute the missing censoring times using the oracle (data-generating) model
        C.cal[idx.event] <- cens_imputator$sample_censoring_times(X.cal[idx.event,], T=Y.cal[idx.event])
    }

    # Apply Candes' method using the imputed censoring times
    out <- predict_Candes(data.test, surv_model, cens_imputator$model, data.cal, C.cal, alpha, c0=c0)
    return(out)
}
