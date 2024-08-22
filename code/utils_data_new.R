# Base class for Survival Data Generation
SurvivalDataGenerator <- R6::R6Class("SurvivalDataGenerator",
  public = list(
    covariate_generator = NULL,
    survival = NULL,
    censoring = NULL,

    initialize = function(covariate_generator, survival_generator, censoring_generator) {
      self$covariate_generator <- covariate_generator
      self$survival <- survival_generator
      self$censoring <- censoring_generator
    },

    sample = function(num_samples) {
      # Generate covariates, true survival times, and censoring times
      X <- self$covariate_generator(num_samples)
      T <- self$survival$sample(X)
      C <- self$censoring$sample(X)

      # Determine observed time and event indicator
      time <- pmin(T, C)
      status <- as.integer(time == T)

      # Create the final dataset
      data <- data.frame(event_time = T, censoring_time = C, time = time,
                         status = status, X)
      return(data)
    }

  )
)


SurvivalDistribution <- R6::R6Class("SurvivalDistribution",
  public = list(
    # Abstract methods
    # Description: Initialize the survival distribution with necessary functions.
    # Inputs: Should be defined by derived classes based on their specific needs.
    initialize = function(...) {
      stop("This method should be implemented by the derived class.")
    },

    # Abstract method to sample survival times based on the distribution
    # Inputs:
    #   - X: A matrix of covariates
    #   - T: (Optional) A vector of true survival times. If provided, times C are sampled conditional on C > T.
    # Outputs:
    #   - A numeric vector of survival times for each individual.
    sample = function(X, T=NULL) {
      stop("This method should be implemented by the derived class.")
    },

    # Abstract method to predict survival probabilities at given time points
    # Inputs:
    #   - X: A matrix of covariates.
    #   - time.points: A numeric vector of time points to predict survival probabilities.
    # Outputs:
    #   - A list containing:
    #       - predictions: A matrix of survival probabilities for each individual at each time point.
    #       - time.points: The time points used for predictions.
    predict = function(X, time.points) {
      stop("This method should be implemented by the derived class.")
    },

    # Abstract method to predict quantiles of survival times
    # Inputs:
    #   - X: A matrix of covariates.
    #   - probs: A numeric vector of probabilities for which to compute the quantiles (default c(0.25, 0.5, 0.75)).
    # Outputs:
    #   - A data frame where each row corresponds to an individual and each column corresponds to a quantile.
    predict_quantiles = function(X, probs = c(0.25, 0.5, 0.75)) {
      stop("This method should be implemented by the derived class.")
    }
  )
)


LogNormalDistribution <- R6::R6Class("LogNormalDistribution",
  inherit = SurvivalDistribution,
  public = list(
      mu_fun = NULL,                 # Function to calculate mean based on covariates X
      sigma_fun = NULL,              # Function to calculate standard deviation based on covariates X

    initialize = function(mu_fun, sigma_fun) {
        self$mu_fun <- mu_fun
        self$sigma_fun <- sigma_fun
    },

    sample = function(X, T=NULL, max_reps=100) {
        ## Ensure that we are working with a data frame excluding time and status if present
        if (is.data.frame(X)) {
            X <- X |> select(-any_of(c("time", "status")))
        }
        num_samples <- nrow(X)
        mu_x <- self$mu_fun(X)
        sigma_x <- self$sigma_fun(X)
        if(is.null(T)) {
            log_out <- rnorm(num_samples, mean = mu_x, sd = sigma_x)
            out <- exp(log_out)
        } else {
            ## Initialize output vector for conditional samples
            out <- rep(NA, num_samples)

            ## Loop through each individual and sample conditionally
            for (i in 1:num_samples) {
                iter <- 0
                while (iter < max_reps) {
                    log_C <- rnorm(1, mean = mu_x[i], sd = sigma_x[i])
                    C <- exp(log_C)

                    ## Check if the sampled value satisfies the condition C > T[i]
                    if (is.null(T) || C > T[i]) {
                        out[i] <- C
                        break  # Exit the loop for this individual when C > T[i]
                    }
                    iter <- iter + 1
                }
                
                ## If maximum iterations exceeded, set the sample equal to T[i]
                if (is.na(out[i]) && !is.null(T)) {
                    out[i] <- T[i]
                    warning(paste("Max iterations exceeded for individual", i, "; sample set to T[i]."))
                }
            }            
        }
        return(out)
    },

    predict = function(X, time.points) {
        ## Ensure that we are working with a data frame excluding time and status if present
        if (is.data.frame(X)) {
            X <- X |> select(-any_of(c("time", "status")))
        }

        ## Compute the mean (mu) and standard deviation (sigma) for each individual
        mu_x <- self$mu_fun(X)
        sigma_x <- self$sigma_fun(X)

        ## Initialize a matrix to store survival probabilities
        log_t <- log(time.points)
        log_t[log_t <= 0] <- -Inf  ## Handle non-positive time points

        ## Vectorized calculation of survival probabilities
        log_t_matrix <- matrix(log_t, nrow = nrow(X), ncol = length(time.points), byrow = TRUE)
        mu_matrix <- matrix(mu_x, nrow = nrow(X), ncol = length(time.points))
        sigma_matrix <- matrix(sigma_x, nrow = nrow(X), ncol = length(time.points))

        ## Compute the survival probabilities for all individuals at all time points
        z_scores <- (log_t_matrix - mu_matrix) / sigma_matrix
        predictions <- 1 - pnorm(z_scores)
        predictions[log_t_matrix == -Inf] <- 1  ## Survival probability is 1 at time 0 or negative times

        return(list(predictions = predictions, time.points = time.points))
    },

    predict_quantiles = function(X, probs = c(0.25, 0.5, 0.75)) {
        ## Ensure that we are working with a data frame excluding time and status if present
        if (is.data.frame(X)) {
            X <- X |> select(-any_of(c("time", "status")))
        }

        ## Compute the mean (mu) and standard deviation (sigma) for each individual
        mu_x <- self$mu_fun(X)
        sigma_x <- self$sigma_fun(X)

        ## Vectorized calculation of quantiles
        qnorm_matrix <- matrix(qnorm(probs), nrow = nrow(X), ncol = length(probs), byrow = TRUE)
        mu_matrix <- matrix(mu_x, nrow = nrow(X), ncol = length(probs))
        sigma_matrix <- matrix(sigma_x, nrow = nrow(X), ncol = length(probs))

        ## Calculate quantiles using the inverse survival function formula
        quantiles_matrix <- exp(mu_matrix + sigma_matrix * qnorm_matrix)

        ## Convert the quantiles matrix to a data frame
        quantiles_df <- as.data.frame(quantiles_matrix)
        colnames(quantiles_df) <- paste0("Q", probs * 100, "%")
        rownames(quantiles_df) <- paste0("Individual_", 1:nrow(quantiles_df))

        return(quantiles_df)
    },

    # Abstract method that does not do anything. (For compatibility)
    fit = function(formula, data) {
    }

    )
  )



ExponentialDistribution <- R6::R6Class("ExponentialDistribution",
  inherit = SurvivalDistribution,
  public = list(
      rate_fun = NULL,                 # Function to calculate exponential rate based on covariates X

    initialize = function(rate_fun) {
        self$rate_fun <- rate_fun
    },

    sample = function(X, T = NULL) {
        ## Ensure that we are working with a data frame excluding time and status if present
        if (is.data.frame(X)) {
            X <- X |> select(-any_of(c("time", "status")))
        }
        num_samples <- nrow(X)
        ## Compute censoring rates
        censoring_rate <- self$rate_fun(X)

        ## Sample times from an exponential distribution
        out <- rexp(num_samples, rate = censoring_rate)

        ## If conditioning values T are provided, shift samples so that C > T
        if (!is.null(T)) {
            out <- out + T
        }

        return(out)
    },

    predict = function(X, time.points) {
        ## Ensure that X is a matrix (this should be the case for survival models)
        if (is.data.frame(X)) {
            X <- X |> select(-any_of(c("time", "status")))
        }

        ## Compute censoring rates using the rate function
        censoring_rate <- self$rate_fun(X)  # Vector of length n_individuals

        ## Compute survival probabilities in a vectorized way
                                        # pexp applied element-wise to time.points across censoring rates
                                        # Expand censoring_rate to match the dimensions of time.points
        censoring_matrix <- matrix(censoring_rate, nrow = length(censoring_rate), ncol = length(time.points), byrow = FALSE)
        time_matrix <- matrix(time.points, nrow = length(censoring_rate), ncol = length(time.points), byrow = TRUE)

        predictions <- pexp(time_matrix, rate = censoring_matrix, lower.tail = FALSE)

        ## Return predictions and time points
        return(list(predictions = predictions, time.points = time.points))
    },


    predict_quantiles = function(X, probs = c(0.25, 0.5, 0.75)) {
        ## Ensure that X is a matrix (this should be the case for survival models)
        if (is.data.frame(X)) {
            X <- X |> select(-any_of(c("time", "status")))
        }

        ## Compute the censoring rates (λ) for each individual based on covariates
        censoring_rate <- self$rate_fun(X)  # λ: a vector of length n_individuals

        ## Compute the exponential quantile function for each probability
        log_term <- -log(1 - probs)  # Vectorized for probs: a vector of length n_probs

        ## Calculate quantiles for all individuals and all probabilities in a vectorized manner
        quantiles_matrix <- outer(1 / censoring_rate, log_term, "*")  # n_individuals x n_probs matrix

        ## Convert the quantiles matrix to a data frame for easier handling
        quantiles_df <- as.data.frame(quantiles_matrix)

        ## Set appropriate column names for quantiles
        colnames(quantiles_df) <- paste0("Q", probs * 100, "%")

        ## Set appropriate row names for individuals
        rownames(quantiles_df) <- paste0("Individual_", 1:nrow(quantiles_df))

        return(quantiles_df)
    },


    ## Abstract method that does not do anything. (For compatibility)
    fit = function(formula, data) {
    }

    )
  )
