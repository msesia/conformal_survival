library(tidyverse)

evaluate_bounds <- function(observed_time, lower_bound, status=NULL, event_time=NULL, oracle=NULL, method=NULL) {
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
    ## Evaluate lower and upper bounds on coverage
    coverage_lower <- mean(lower_bound <= observed_time)
    if(!is.null(status)) {
        idx.status <- which(status==TRUE)
        if(length(idx.status)>0) {
            coverage_upper <- 1-mean(lower_bound[idx.status] > observed_time[idx.status]) * mean(status==TRUE)
        } else {
            coverage_upper <- NA
        }
    } else {
        coverage_upper <- NA
    }
    if(is.null(oracle)) {
        oracle_MSE <- NA
    } else {
        oracle_MSE <- mean((lower_bound-oracle)^2)
    }
    mean_lower_bound <- mean(lower_bound)
    median_lower_bound <- median(lower_bound)
    out <- tibble(#"Coverage (observed time)"=coverage_observed,
                  "Mean lower bound"=mean_lower_bound,
#                  "Median lower bound"=median_lower_bound,
#                  "Mean lower bound (cover)"=mean_lower_bound_cover,
#                  "Median lower bound (cover)"=median_lower_bound_cover,
                  "Coverage (event time)"=coverage_event_time,
                  "Coverage (lower bound)"=coverage_lower,
                  "Coverage (upper bound)"=coverage_upper,
#                  "Oracle MSE" = oracle_MSE
    )
    if(!is.null(method)) {
        out <- out %>% mutate(Method = method) %>% select(Method, everything())
    }
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
    if(percentile_scores<1) {
        calibration <- as.numeric(quantile(scores, probs = percentile_scores, type = 1))
    } else {
        calibration <- Inf
    }

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


impute_drcosarc <- function(data.cal, cens_imputator) {
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

    return(C.cal)
}

predict_drcosarc <- function(data.test, surv_model, cens_imputator, data.cal, alpha, c0=NULL, tuning.package=NULL, cutoffs="adaptive",
                              finite_sample_correction = FALSE, doubly_robust = TRUE, alpha_hat_values = NULL) {

    ## Impute the missing censoring times for calibration data
    C.cal <- impute_drcosarc(data.cal, cens_imputator)

    if(cutoffs=="adaptive") {
        ## Apply Gui's method using the imputed censoring times
        pred.calibrated <- predict_Gui(data.test, surv_model, cens_imputator$model, data.cal, C.cal, alpha, finite_sample_correction=finite_sample_correction,
                                       alpha_hat_values = alpha_hat_values)
    } else if (cutoffs=="candes-fixed") {
        ## Apply Candes' method using the imputed censoring times
        pred.calibrated <- predict_Candes(data.test, surv_model, cens_imputator$model, data.cal, C.cal, alpha, c0=c0)
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
        pred.calibrated <- predict_Candes(data.test, surv_model, cens_imputator$model, data.cal, C.cal, alpha, tuning.package=tuning.package, c0=c0)
    }
    else {
        stop(sprintf("Unknown input option: cutoffs = %s.", cutoffs))
    }

    if(doubly_robust) {
        ## No calibration (trust the survival model's output)
        pred.nominal <- surv_model$predict_quantiles(data.test, probs=alpha)[[1]]
        out <- pmin(pred.calibrated, pred.nominal)
    } else {
        out <- pred.calibrated
    }

    return(out)
}

generate_probs_seq <- function(num_a) {
    probs <- seq(0.01, 0.99, length.out=num_a)
    return(probs)
}

compute_alpha_hat_Gui <- function(data.cal, C.cal, surv_model, cens_model, num_a, shift=0, use_censoring_model=TRUE, censoring_beta=NULL,
                                  use_cqr=FALSE, alpha_cqr=NULL) {
    Y.cal <- data.cal$time
    num_cal <- nrow(data.cal)

    if(use_censoring_model) stopifnot(!is.null(censoring_beta))

    ## List of percentiles for estimated survival distribution
    probs <- generate_probs_seq(num_a)

    ## Estimate miscoverage as a function of calibration parameter a
    estimate_alpha_hat <- function(num_a) {
        ## Construct sequences of predictions based on parameters a
        if(!use_cqr) {
            ## Predict the survival quantiles for the given nominal percentiles
            pred.cal <- as.matrix(pmax(surv_model$predict_quantiles(data.cal, probs = probs) - shift, 0))
            if(use_censoring_model) {
                pred.C.cal <- as.numeric(cens_model$predict_quantiles(data.cal, probs = 1-censoring_beta)[[1]])
                pred.cal <- pmin(pred.cal, pred.C.cal)                
            }
        } else {
            stopifnot(!is.null(alpha_cqr))
            shifts <- seq(0, 1, length.out=num_a)
            pred.cal.raw <- as.numeric(surv_model$predict_quantiles(data.cal, probs = alpha_cqr)[[1]])
            pred.cal.raw <- matrix(pred.cal.raw, nrow = num_cal, ncol = length(shifts))
            shifts.mat <- matrix(shifts, num_cal, length(shifts), byrow = TRUE)
            pred.cal <- pmax(pred.cal.raw * shifts.mat, 0)
        }

        ## Compute the weights a function of a
        weights.cal.mat <- matrix(0, num_cal, num_a)
        for(i in 1:num_cal) {
            c0_seq <- as.numeric(pred.cal[i,])
            cens_surv_probs <- cens_model$predict(data.cal[i,], time.points=c0_seq)$predictions
            weight_tol <- 1e-4
            idx.too.small <- which(cens_surv_probs<weight_tol)
            if(length(idx.too.small)>0) {
                ## RAISE WARNING IF INVERSE WEIGHTS ARE TOO SMALL!
                if(use_censoring_model) {
                    stop("Error: unstable weights! Try using option use_censoring_model=TRUE.\n")
                } else {
                    stop("Error: unstable weights! Try using option use_censoring_model=TRUE.\n")
                }
            }
            weights.cal.mat[i,] <- 1/cens_surv_probs
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

        return(alpha_hat_values)
    }

    alpha_hat_values <- estimate_alpha_hat(num_a)

    return(alpha_hat_values)
}

predict_Gui <- function(data.test, surv_model, cens_model, data.cal, C.cal, alpha, shift=0, use_censoring_model=TRUE, censoring_beta=NULL,
                        finite_sample_correction = FALSE, alpha_hat_values = NULL, use_cqr=FALSE) {

    Y.cal <- data.cal$time
    num_cal <- nrow(data.cal)
    num_test <- nrow(data.test)

    if(is.null(censoring_beta)) {
        # Default value for censoring beta constant (for numerical stability)
        censoring_beta = 1/pmax(1,log(num_cal)) # This is the value suggested by Gui et al., which we used in the paper        
    }

    if(is.null(alpha_hat_values)) {
        ## Number of distinct values for tuning parameter alpha
        num_a <- 200

        ## List of percentiles for estimated survival distribution
        probs <- generate_probs_seq(num_a)

        ## Compute alpha_hat values
        alpha_hat_values <- compute_alpha_hat_Gui(data.cal, C.cal, surv_model, cens_model, num_a, shift=shift, use_censoring_model=use_censoring_model,
                                                  censoring_beta = censoring_beta, use_cqr=use_cqr, alpha_cqr=alpha)
    } else {
        num_a <- length(alpha_hat_values)
        ## List of percentiles for estimated survival distribution
        probs <- generate_probs_seq(num_a)
    }

    ## Find the smallest a such that alpha_hat(a) < alpha
    if(finite_sample_correction) {
        alpha_adjusted <- alpha - (1-alpha)/num_cal - (1) / sqrt(num_cal)
        alpha_adjusted <- pmax(0, alpha_adjusted)
    } else {
        alpha_adjusted <- alpha
    }

    idx.valid <- which(alpha_hat_values<=alpha_adjusted)

    ## Make predictions for test data
    if(!use_cqr) {
        pred.test <- as.matrix(pmax(surv_model$predict_quantiles(data.test, probs = probs) - shift, 0))
        if(use_censoring_model) {
            pred.C.test <- as.numeric(cens_model$predict_quantiles(data.test, probs = 1-censoring_beta)[[1]])
            pred.test <- pmin(pred.test, pred.C.test)
        }
    } else {
        ## List of shifts to correct CQR predictions
        shifts <- seq(0, 1, length.out=num_a)
        shifts.mat <- matrix(shifts, num_test, length(shifts), byrow = TRUE)
        pred.test.raw <- as.numeric(surv_model$predict_quantiles(data.test, probs = alpha)[[1]])
        pred.test.raw <- matrix(pred.test.raw, nrow = num_test, ncol = length(shifts))
        shifts.mat <- matrix(shifts, num_test, length(shifts), byrow = TRUE)
        pred.test <- pmax(pred.test.raw * shifts.mat, 0)
    }

    if(length(idx.valid)>0) {
        a.star <- max(idx.valid)
        ## Compute the prediction corresponding to a.star
        lower <- as.numeric(pred.test[,a.star])
    } else {
        if(use_cqr) {
            ## Be careful that in this case a.star may not lead to valid predictions
            a.star <- num_a
            lower <- as.numeric(pred.test[,a.star]) * 0
        } else {
            print("Warning!")
            lower <- predict_Gui(data.test, surv_model, cens_model, data.cal, C.cal, alpha, shift=0, use_cqr=TRUE, finite_sample_correction = finite_sample_correction)
        }
    }

    return(lower)
}
