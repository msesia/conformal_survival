CensoringModel <- R6::R6Class("CensoringModel",
  public = list(
    model = NULL,  # Model for estimating the censoring distribution

    # Constructor to initialize the class with a survival model
    initialize = function(model) {
      self$model <- model
    },

    # Method to fit the censoring distribution using the provided survival model
    # Inputs:
    #   - model: The survival model used to estimate the censoring distribution.
    #   - data: The dataset containing survival information with event indicators.
    # Description: This method flips the event indicators and fits the censoring distribution.
    fit = function(data) {
      data.cens <- data

      # Flip the event indicators: 1 if censored, 0 if event
      data.cens$status <- 1 - data.cens$status

      # Fit the censoring model using the survival data
      self$model$fit(survival::Surv(time, status) ~ ., data = data.cens)
    },

    # Method to sample censoring times based on the fitted censoring distribution
    # Inputs:
    #   - new_data: A data.frame containing covariates for the individuals.
    #   - T: (Optional) A vector of true survival times. If provided, censoring times are sampled such that C > T.
    #   - max_reps: (Optional) The maximum number of sampling repetitions to enforce the C > T condition.
    # Outputs:
    #   - A vector of sampled censoring times (C) for each individual. If T is provided, C > T.
    sample = function(new_data, T = NULL, max_reps = 100) {
      # Predict survival curves for the covariates X using the fitted censoring distribution model
      predictions <- self$model$predict(new_data = new_data)

      # Extract the survival probabilities and failure times
      survival_probs <- predictions$predictions
      failure_times <- predictions$time.points

      # Function to sample missing censoring times based on the estimated censoring distribution
      sample_time <- function(surv_probs, failure_times, min_time = -Inf, max_reps = 100) {
        valid_time <- FALSE
        sampled_time <- NA
        reps <- 0  # Initialize the repetition counter

        # Keep sampling until a valid time greater than min_time is found, or the max number of repetitions is exceeded
        while (!valid_time && reps < max_reps) {
          u <- runif(1)  # Generate a uniform random number

          # Find the index where the survival probability is greater than or equal to u
          valid_indices <- which(surv_probs >= u)

          # Check if valid_indices is not empty before calling max()
          if (length(valid_indices) > 0) {
              idx <- max(valid_indices)
          } else {
              idx <- length(failure_times)  # Or handle the case when no valid index is found
          }

          # Get the corresponding failure time
          sampled_time <- failure_times[idx]

          # Check if the sampled time is valid (greater than min_time)
          if (sampled_time > min_time) {
            valid_time <- TRUE
          }

          reps <- reps + 1  # Increment repetition counter
        }

        # If the maximum number of repetitions is exceeded, return T (true survival time)
        if (!valid_time) {
            warning("Max repetitions exceeded; returning T as the censoring time.")
            return(min_time)  # Return the true survival time as the censoring time
        }

        return(sampled_time)
      }

      # Initialize the censoring times vector
      C <- numeric(nrow(new_data))

      # Loop through each individual to sample censoring times
      for (i in 1:nrow(new_data)) {
        # If true survival times are provided, ensure C > T
        if (!is.null(T)) {
          C[i] <- sample_time(survival_probs[i, ], failure_times, min_time = T[i], max_reps = max_reps)
        } else {
          C[i] <- sample_time(survival_probs[i, ], failure_times, max_reps = max_reps)
        }
      }

      return(C)
    }

  )
)
