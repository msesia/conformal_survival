library(tidyverse)

evaluate_bounds <- function(observed_time, lower_bound, event_time=NULL, oracle=NULL) {
    coverage_observed <- mean(lower_bound <= observed_time)
    coverage_event_time <- NA
    mean_lower_bound_cover <- NA
    median_lower_bound_cover <- NA
    if(!is.null(event_time)) {
        coverage_event_time <- mean(lower_bound <= event_time)
        idx.cover <- which(lower_bound <= event_time)
        if(length(idx.cover)>0) {
            mean_lower_bound_cover <- mean(lower_bound[idx.cover])
            median_lower_bound_cover <- median(lower_bound[idx.cover])
        }
   }
    if(is.null(oracle)) {
        oracle_MSE <- NA
    } else {
        oracle_MSE <- mean((lower_bound-oracle)^2)
    }
    mean_lower_bound <- mean(lower_bound)
    median_lower_bound <- median(lower_bound)
    out <- tibble("Coverage (observed time)"=coverage_observed,
                  "Mean lower bound"=mean_lower_bound,
                  "Median lower bound"=median_lower_bound,
                  "Mean lower bound (cover)"=mean_lower_bound_cover,
                  "Median lower bound (cover)"=median_lower_bound_cover,
                  "Coverage (event time)"=coverage_event_time,
                  "Oracle MSE" = oracle_MSE)
    return(out)

}

predict_CQR <- function(data.test, surv_model, data.cal, alpha) {
    # Which quantile to predict? alpha is a reasonable choice,
    # but in theory any other value can be used
    probs <- c(alpha)

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

    return(lower.test)
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
  idx.above <- which(cumsum(sort_w) >= 1 - alpha)
  if(length(idx.above)>0) {
      idxw <- min(idx.above)
      ## Compute the calibration term as the score at the index determined by the cumulative weights
      calib_term <- sort_score[idxw]
  } else {
      calib_term <- Inf
  }


  return(calib_term)
}



predict_Candes <- function(data.test, surv_model, cens_model, data.cal, C.cal, alpha, tuning.package=NULL, c0=NULL) {
    # Which quantile to predict? alpha is a reasonable choice,
    # but in theory any other value can be used
    probs <- c(alpha)

    ## Function for tuning c0
    tune.c0 <- function(tuning.package) {
        ## Extract relevant stuff from tuning package
        data.train <- tuning.package$data.train
        surv_model_tune <- surv_model$clone(deep = TRUE)
        cens_model_tune <- cens_model$clone(deep = TRUE)
        C.train = tuning.package$C.train

        ## Simulate the lower bounds for different values of c0
        c0.candidates <- as.numeric(quantile(C.cal, c(seq(0.05,0.95,by=0.05))))

        greedy.calibration <- TRUE

        if(greedy.calibration) {
            lb.candidates <- sapply(c0.candidates, function(c0) {
                lower.tune <- predict_Candes(data.test, surv_model_tune, cens_model_tune, data.cal, C.cal, alpha, c0=c0)
                return(median(lower.tune))
            })
        } else {            
            ## Split the training data into two subsets
            train_indices.1 <- sample(1:nrow(data.train), size = 0.5 * nrow(data.train))
            data.train.1 <- data.train[train_indices.1, ]
            data.train.2 <- data.train[-train_indices.1, ]
            C.train.2 <- C.train[-train_indices.1]
            ## Fit the survival model on the training(1) data
            surv_model_tune$fit(survival::Surv(time, status) ~ ., data = data.train.1)
            ## Fit the censoring model on the training data
            cens_model_tune$fit(survival::Surv(time, 1-status) ~ ., data = data.train.1)
            if(FALSE) {
                ## Use bootstrap to make this more stable        
                B = 10
                lb.candidates.boot <- do.call(rbind, lapply(1:B, function(i) {
                    idx.test <- sample(nrow(data.test), 10, replace=TRUE)
                    idx.cal <- sample(nrow(data.train.2), nrow(data.train.2), replace=TRUE)
                    median.lb.candidates <- sapply(c0.candidates, function(c0) {
                        lower.tune <- predict_Candes(data.test[idx.test,], surv_model_tune, cens_model_tune, data.train.2[idx.cal,], C.train.2[idx.cal], alpha, c0=c0)
                        return(median(lower.tune))
                    })
                }))
                lb.candidates <- colMeans(lb.candidates.boot)
            } else {
                lb.candidates <- sapply(c0.candidates, function(c0) {
                    lower.tune <- predict_Candes(data.test, surv_model_tune, cens_model_tune, data.train.2, C.train.2, alpha, c0=c0)
                    return(median(lower.tune))
                })
            }
        }
        idx.selected <- which(diff(lb.candidates) < 0)[1]
        if(length(idx.selected)>0) {
            c0 <- c0.candidates[idx.selected]
        } else {
            c0 <- c0.candidates[1]
        }
        return(c0)
    }

    if((is.null(c0)) && (!is.null(tuning.package))) {
        c0 <- tune.c0(tuning.package)
        cat(sprintf("Tuned c0 = %.3f\n", c0))
    } else {
        ## Choose c0 if not supplied
        if(is.null(c0)) {
            c0 <- quantile(C.cal, 0.5)
        }
    }    

    # Extract the covariate information (remove 'time' and 'status' columns, if present)
    ##X.cal <- data.cal |> select(-any_of(c("time", "status")))
    ##X.test <- data.test |> select(-any_of(c("time", "status")))
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
    weights.cal <- 1/pmax(1e-6, cens_model$predict(data.cal[idx.keep,], time.points=c0)$predictions)
    weights.test <- 1/pmax(1e-6, cens_model$predict(data.test, time.points=c0)$predictions)
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
        lower <- pmin(c0, pmax(0, pmin(pred.test[i] - calibration, c0)))
        return(lower)
    })
    
    return(lower.test)
}


predict_prototype <- function(data.test, surv_model, cens_imputator, data.cal, alpha, c0=NULL, tuning.package=NULL, cutoffs="adaptive") {
    # Initialize the censoring times equal to the observed times
    C.cal <- data.cal$time

    ## Extract the covariate information (remove 'time' and 'status' columns, if present)
    ##X.cal <- data.cal |> select(-any_of(c("time", "status")))
    Y.cal <- data.cal$time

    ## Impute the missing censoring times for calibration data
    ## Find the indices of observations for which the censoring time needs to be imputed
    idx.event <- which(data.cal$status==TRUE)

    if(length(idx.event) > 0) {
        # Impute the missing censoring times using the oracle (data-generating) model
        C.cal[idx.event] <- cens_imputator$sample(data.cal[idx.event,], T=Y.cal[idx.event])
    }

    if(cutoffs=="adaptive") {
        ## Apply Gui's method using the imputed censoring times
        out <- predict_Gui(data.test, surv_model, cens_imputator$model, data.cal, C.cal, alpha)
    } else if(cutoffs=="adaptive-cqr") {
        ## Apply Gui's method using the imputed censoring times
        out <- predict_Gui(data.test, surv_model, cens_imputator$model, data.cal, C.cal, alpha, use_cqr=TRUE)
    } else if (cutoffs=="candes-fixed") {
        ## Apply Candes' method using the imputed censoring times
        out <- predict_Candes(data.test, surv_model, cens_imputator$model, data.cal, C.cal, alpha, c0=c0)
    } else if (cutoffs=="candes-tuning") {
        if(is.null(tuning.package)) {
            stop("Must provide 'tuning.package' input to allow automatic tuning of c0.")
        }
        ## Impute the missing censoring times for training data
        idx.event.train <- which(tuning.package$data.train$status==TRUE)
        C.train <- tuning.package$data.train$time
        if(length(idx.event.train) > 0) {
            Y.train <- tuning.package$data.train$event
            C.train[idx.event.train] <- cens_imputator$sample(tuning.package$data.train[idx.event.train,], T=Y.train[idx.event])
        }
        tuning.package$C.train <- C.train
        out <- predict_Candes(data.test, surv_model, cens_imputator$model, data.cal, C.cal, alpha, tuning.package=tuning.package, c0=c0)
    }
    else {
        stop(sprintf("Unknown input option: cutoffs = %s.", cutoffs))
    }
    return(out)
}


predict_Gui <- function(data.test, surv_model, cens_model, data.cal, C.cal, alpha, shift=0, use_cqr=FALSE, use_censoring_model=FALSE,
                        finite_sample_correction = TRUE) {

    ## Number of distinct values for tuning parameter alpha
    num_a <- 200

    Y.cal <- data.cal$time
    num_cal <- nrow(data.cal)
    num_test <- nrow(data.test)

    ## Construct sequences of predictions based on parameters a

    if(!use_cqr) {
        ## List of percentiles for estimated survival distribution
        probs <- seq(0.01, 0.99, length.out=num_a)

        ## Predict the survival quantiles for the given nominal percentiles
        if(use_censoring_model) {
            pred.T.cal <- pmax(surv_model$predict_quantiles(data.cal, probs = probs) - shift, 0)
            pred.T.test <- pmax(surv_model$predict_quantiles(data.test, probs = probs) - shift,0 )
            beta = 1/pmax(1,log(num_cal))
            pred.C.cal <- cens_model$predict_quantiles(data.cal, probs = 1-beta)[[1]]
            pred.C.test <- cens_model$predict_quantiles(data.test, probs = 1-beta)[[1]]
            pred.cal <- t(sapply(1:nrow(pred.T.cal), function(i) pmin(pred.T.cal[i, ], pred.C.cal[i])))
            pred.test <- t(sapply(1:nrow(pred.T.test), function(i) pmin(pred.T.test[i, ], pred.C.test[i])))
        } else {
            pred.cal <- pmax(surv_model$predict_quantiles(data.cal, probs = probs) - shift, 0)
            pred.test <- pmax(surv_model$predict_quantiles(data.test, probs = probs) - shift, 0)
        }
    } else {
        pred.cal.raw <- as.numeric(surv_model$predict_quantiles(data.cal, probs = alpha)[[1]])
        max.shift <- max(pred.cal.raw)
        #shifts <- rev(seq(0, max.shift, length.out=num_a))
        shifts <- seq(0, 1, length.out=num_a)
        shifts.mat <- matrix(shifts, num_cal, length(shifts), byrow = TRUE)
        pred.cal.raw <- matrix(pred.cal.raw, nrow = length(pred.cal.raw), ncol = length(shifts))
        pred.cal <- pmax(pred.cal.raw * shifts.mat, 0)
        pred.test.raw <- as.numeric(surv_model$predict_quantiles(data.test, probs = alpha)[[1]])
        pred.test.raw <- matrix(pred.test.raw, nrow = length(pred.test.raw), ncol = length(shifts))
        shifts.mat <- matrix(shifts, num_test, length(shifts), byrow = TRUE)
        pred.test <- pmax(pred.test.raw * shifts.mat, 0)
    }

    ## Compute the weights a function of a
    weights.cal.mat <- matrix(0, num_cal, num_a)
    for(i in 1:num_cal) {
        c0_seq = as.numeric(pred.cal[i,])
        weights.cal.mat[i,] <- 1/pmax(1e-6, cens_model$predict(data.cal[i,], time.points=c0_seq)$predictions)
    }

    ## Estimate alpha as a function of a
    compute_alpha_hat <- function(a) {
        weights.cal <- weights.cal.mat[,a]
        idx.keep.num <- which((Y.cal<pred.cal[,a])*(pred.cal[,a]<=C.cal)==TRUE)
        idx.keep.den <- which((pred.cal[,a]<=C.cal)==TRUE)
        if(length(idx.keep.den)==0) return(1)
        if(length(idx.keep.num)==0) return(0)
        out <- sum(weights.cal[idx.keep.num]) / sum(weights.cal[idx.keep.den])
        return(out)
    }
    alpha_hat_values <- sapply(1:num_a, function(a) compute_alpha_hat(a))
    alpha_hat_values <- cummax(alpha_hat_values)

    ## Find the smallest a such that alpha_hat(a) < alpha
    if(finite_sample_correction) {
        alpha_adjusted <- alpha - (1-alpha)/num_cal - (1) / sqrt(num_cal)
        alpha_adjusted <- pmax(0, alpha_adjusted)
    } else {
        alpha_adjusted <- alpha
    }
    idx.valid <- which(alpha_hat_values<=alpha_adjusted)
    if(length(idx.valid)>0) {
        a.star <- max(idx.valid)
        ## Compute the prediction corresponding to a.star
        lower <- as.numeric(pred.test[,a.star])
    } else {
        a.star <- num_a
        ## Be careful that in this case a.star may not lead to valid predictions
        lower <- as.numeric(pred.test[,a.star]) * 0
    }

    return(lower)
}
