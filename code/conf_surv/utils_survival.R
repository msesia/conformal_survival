# Function to ensure `new_data` is treated as a matrix
# If `new_data` is a vector, it will be converted to a matrix with one row.
# If `new_data` is already a matrix, it will be returned unchanged.
# If `new_data` is neither a vector nor a matrix, an error is raised.
ensure_matrix <- function(new_data) {
  # Check if the input is a vector
  if (is.vector(new_data)) {
    # Convert the vector to a matrix with one row
    new_data <- matrix(new_data, nrow = 1)

  # Check if the input is already a matrix
  } else if (is.matrix(new_data)) {
    # If it's already a matrix, return it unchanged
    return(new_data)

  # Raise an error if the input is neither a vector nor a matrix
  } else {
    stop("Input must be either a vector or a matrix")
  }

  # Return the resulting matrix
  return(new_data)
}

# Function to interpolate survival probability at specific time points using monotone interpolation
get_survival_prob_at_time <- function(pred, time_points) {
  # Extract the survival times and survival probabilities from the prediction object
  times <- pred$time.interest
  survival_probs <- pred$survival
  
  # Initialize a matrix to hold the interpolated probabilities
  num_individuals <- nrow(survival_probs)
  num_time_points <- length(time_points)
  interpolated_probs <- matrix(NA, nrow = num_individuals, ncol = num_time_points)
  
  # Perform monotone interpolation using PCHIP for each individual and each time point
  for (i in 1:num_individuals) {
    # Create the monotone interpolation function for the current individual
    interp_function <- splinefun(x = times, y = survival_probs[i, ], method = "monoH.FC")
    
    # Interpolate at each time point and ensure the probabilities are within [0, 1]
    interpolated_probs[i, ] <- pmin(pmax(interp_function(time_points), 0), 1)
  }
  
  return(interpolated_probs)
}


SurvivalModelWrapper <- R6::R6Class("SurvivalModelWrapper",
  public = list(
    model = NULL,                # Holds the trained model object.
    formula = NULL,              # Stores the formula used to fit the model.
    time.points = NULL,          # A sequence of time points for which survival probabilities are calculated.
    use_covariates = NULL,       # List of relevant covariates (e.g., c("X1", "X3")) to use when fitting the censoring model

    # Constructor
    # Description: Default constructor
    initialize = function(use_covariates = NULL) {
      self$use_covariates <- use_covariates
    },

    ## Abstract method to fit a survival model
    # Description:
    #   This method is abstract and must be implemented in the subclass. It is responsible for fitting
    #   a survival model using the specified formula and data. The subclass should define the specific
    #   survival model to be fitted (e.g., parametric survival model using `survreg`, Cox proportional hazards,
    #   random survival forest). The fitted model should be stored in the `model` field of the class for later use
    #   in prediction methods. After fitting the model, relevant time points (e.g., failure times) should also be
    #   stored in 'self$time.points' for use in survival predictions.
    #
    # Inputs:
    #   - formula: A survival formula of the form `Surv(time, status) ~ predictors`. The formula should
    #              specify the time-to-event outcome (`Surv(time, status)`) and the predictor variables.
    #   - data: A data.frame containing the predictor variables, the time variable, and the status variable.
    #           The data frame must have a column for the time-to-event (e.g., `time`), a binary status
    #           variable indicating the event occurrence (e.g., `status`), and the predictor variables.
    #
    # Outputs:
    #   - None. The method does not return a value, but the fitted model should be stored in the `model` field
    #     of the class for subsequent use in prediction methods like `predict` or `predict_quantiles`.
    #
    # Notes for Subclass Implementation:
    #   - The subclass should implement this method to define the specific survival model being used.
    #   - Once the model is fitted, it should be stored in the `self$model` field, making it accessible for
    #     other methods (e.g., `predict`, `predict_quantiles`).
    #   - If the model involves failure times (e.g., parametric models like `survreg`), these times can also
    #     be stored in a field such as `self$time.points` to be used later for survival predictions.
    #
    # Example of Subclass Implementation Using `survreg`:
    #   - In a subclass using a parametric survival model via `survreg`:
    #     fit = function(formula, data, dist = "weibull") {
    #       # Fit a parametric survival model using the `survreg` function from the `survival` package
    #       self$model <- survival::survreg(formula, data = data, dist = dist)
    #
    #       # Store the range of observed times to help generate predictions later
    #       self$time.points <- seq(min(data$time), max(data$time), length.out = 100)
    #     }
    #
    fit = function(formula, data) {
      stop("This method should be implemented in the subclass.")
    },


    # Default method to predict survival curves
    # Description:
    #   This method predicts survival curves for a given set of new data. It is intended to be used
    #   by subclasses of `SurvivalModelWrapper` and relies on the `predict_quantiles()` method,
    #   which should be defined in the subclass. If the `predict_quantiles()` method is not defined,
    #   this method should be implemented in the subclass. The method supports interpolating survival
    #   probabilities at specified custom failure times. If no custom failure times are provided, the default
    #   failure times from the model are used.
    #
    #   Note: At least one of `predict` or `predict_quantiles` needs to be implemented in the subclass.
    #
    #
    # Inputs:
    #   - new_data: A data.frame containing new predictor variables for which to predict survival curves.
    #               The data frame must have the same structure as the data used to train the model.
    #   - time.points: (Optional) A numeric vector of custom failure times at which survival probabilities
    #                    should be interpolated. If NULL, the default failure times from the model are used.
    #
    # Outputs:
    #   - A list with two components:
    #     - `predictions`: A matrix of survival probabilities. Each row corresponds to an individual from
    #                      `new_data`, and each column corresponds to a failure time. The value in each
    #                      cell represents the survival probability at the corresponding failure time.
    #     - `time.points`: A numeric vector of the failure times at which survival probabilities were
    #                        calculated (either the provided custom times or the model's default times).
    #
    # Steps:
    #   1. The function begins by checking if `time.points` is provided. If not, it defaults to using
    #      `self$time.points`, which should have been set during model fitting. If `time.points` is
    #      still NULL, an error is raised.
    #
    #   2. It generates survival quantiles for the individuals in `new_data` by calling `self$predict_quantiles()`.
    #      This method should return the predicted survival times (quantiles) for each probability in the
    #      specified `probs` sequence (default: `seq(0.01, 0.99, by = 0.01)`).
    #
    #   3. A matrix `survival_probs` is initialized to store the interpolated survival probabilities for
    #      each individual at the provided `time.points`.
    #
    #   4. For each individual, the function constructs an interpolation function using `approxfun()`
    #      to perform linear monotone decreasing interpolation between the predicted quantiles and their
    #      corresponding probabilities. The survival probabilities are then interpolated at the specified
    #      `time.points`.
    #
    #   5. Finally, the method returns a list containing the matrix of predicted survival probabilities
    #      (`predictions`) and the vector of failure times (`time.points`).
    #
    # Example Usage:
    #   survival_predictions <- model$predict(new_data = test_data)
    #   survival_predictions_custom <- model$predict(new_data = test_data, time.points = c(100, 200, 300))
    #
    predict = function(new_data, time.points = NULL) {
      # If time.points is not provided, use the default values from the model
      if (is.null(time.points)) {
        time.points <- self$time.points
      }

      # Ensure time.points is provided either by the user or from the model
      if (is.null(time.points)) {
        stop("Error: time.points are not provided, and there is no default value set.")
      }

      # Step 1: Predict quantiles for each probability in `probs`
      probs <- seq(0.01, 0.99, by = 0.01)
      survival_times <- self$predict_quantiles(new_data, probs = probs)

      # Step 2: Initialize a matrix to hold survival probabilities
      survival_probs <- matrix(NA, nrow = nrow(survival_times), ncol = length(time.points))
        
      # Step 3: Loop through each individual and interpolate the survival probabilities
      for (i in 1:nrow(survival_times)) {       
        if(!all(is.na(survival_times[i, ]))) {
          # Create an interpolation function for each individual
          interp_fun <- approxfun(rev(as.numeric(survival_times[i, ])), rev(probs), rule = 2, method = "linear", ties = "ordered")
          # Apply the interpolation function to the provided time.points
          survival_probs[i, ] <- interp_fun(time.points)
        } else {
            survival_probs[i, ] <- 1
        }
      }

      # Step 4: Return the predicted survival curves and failure times
      list(predictions = survival_probs, time.points = time.points)
    },


    # Default method to predict survival quantiles
    # Description:
    #   This method predicts survival quantiles for a given set of new data. It uses the predicted
    #   survival curves generated by the `predict` method, so the `predict` method
    #   must be implemented in the subclass. If `predict_quantiles` is not overridden in the subclass,
    #   this default implementation will be used.
    #
    #   Note: At least one of `predict` or `predict_quantiles` needs to be implemented in the subclass.
    #
    # Inputs:
    #   - new_data: A data.frame containing new predictor variables for which to predict survival quantiles.
    #               The data frame must have the same structure as the data used to train the model.
    #   - probs: A numeric vector of quantile probabilities (e.g., `c(0.25, 0.5, 0.75)` for the 25th, 50th,
    #            and 75th percentiles). These represent the survival probabilities at which the quantile
    #            survival times will be calculated.
    #
    # Outputs:
    #   - A data.frame where each row corresponds to an individual in `new_data`, and each column corresponds
    #     to one of the requested quantiles. The values in the data frame represent the survival times (in terms
    #     of time points) at which the specified survival probabilities (`probs`) are reached.
    #
    # Method Overview:
    #   1. The method begins by calling `self$predict(new_data)` to obtain the survival curves
    #      for each individual in `new_data`. The survival curves represent the survival probabilities
    #      at various time points.
    #
    #   2. The method defines an internal function `find_quantile` that, for each individual, locates the
    #      first time point where the survival probability drops below the specified percentile (1 - quantile).
    #      The corresponding time point is recorded as the quantile survival time.
    #
    #   3. The method loops over each individual, applying `find_quantile` to compute the survival times
    #      corresponding to the specified quantiles (`probs`).
    #
    #   4. Finally, the method returns a data frame where each row contains the predicted quantiles for
    #      each individual, and each column corresponds to a specific quantile probability (e.g., Q25%, Q50%, Q75%).
    #
    # Example Usage:
    #   quantiles_df <- model$predict_quantiles(new_data = test_data, probs = c(0.25, 0.5, 0.75))
    #   quantiles_df <- model$predict_quantiles(new_data = test_data)  # Uses default probs: c(0.25, 0.5, 0.75)
    #
    predict_quantiles = function(new_data, probs = c(0.25, 0.5, 0.75)) {
        ## Predict survival curves
        predictions <- self$predict(new_data)
        survival_curves <- predictions$predictions
        time_points <- predictions$time.points  # Time points associated with the survival curves

        # Add padding to ensure interpolation will work
        survival_curves <- cbind(1,survival_curves,0)
        time_points <- c(0,time_points,max(time_points)+1)
        
        ## Function to find the survival time corresponding to a given survival percentile using built-in interpolation
        find_quantile <- function(survival_probs, time_points, percentile) {
            target_prob <- 1 - percentile  # Convert percentile to survival probability threshold
            ## Use linear interpolation
            interpolated_time <- approx(x = rev(survival_probs), y = rev(time_points), xout = target_prob, rule = 2, ties="ordered")$y
            return(interpolated_time)
        }
       
        ## Initialize a list to store quantiles for each individual
        quantiles_list <- list()

        ## Loop over each individual
        for (i in 1:nrow(survival_curves)) {
            # For each individual, find the survival times at the specified percentiles
            quantiles <- sapply(probs, function(p) find_quantile(survival_curves[i, ], time_points, p))
            quantiles_list[[i]] <- quantiles
        }

        ## Convert the list of quantiles to a data frame
        quantiles_df <- do.call(rbind, quantiles_list)
        colnames(quantiles_df) <- paste0("Q", probs * 100, "%")
        rownames(quantiles_df) <- paste0("Individual_", 1:nrow(quantiles_df))

        return(as.data.frame(quantiles_df))
    },


    # Predict survival with interpolation of survival probabilities at custom failure times (utility)
    predict_interpolate = function(survival_probs, original_failure_times, time.points) {
      # Initialize a matrix to store interpolated survival probabilities
      survival_probs_interp <- matrix(NA, nrow = nrow(survival_probs), ncol = length(time.points))

      # Loop through each individual and interpolate the survival probabilities at the custom times
      for (i in 1:nrow(survival_probs)) {
        # Create an interpolation function for each individual
        interp_fun <- approxfun(original_failure_times, survival_probs[i, ], rule = 2, ties = "ordered")

        # Apply the interpolation function to the custom failure times
        survival_probs_interp[i, ] <- interp_fun(time.points)
      }

      # Return the interpolated survival probabilities
      return(survival_probs_interp)
    },

    # Parse formula method (utility)
    # Description: Parses a survival formula and extracts the response variables (time and status) and covariates (predictors).
    # Input:
    #   - formula: A survival formula of the form `Surv(time, status) ~ predictors`.
    #   - data: A data.frame containing the time, status, and predictor variables.
    # Output: A list with three components:
    #   - time: A vector of survival times.
    #   - status: A vector indicating whether the event (death) occurred (1) or was censored (0).
    #   - covariates: A matrix of covariate data (predictor variables), excluding the intercept.
    parse_formula = function(formula, data) {
      self$formula <- formula
      response <- model.response(model.frame(formula, data))
      time <- response[, 1]
      status <- response[, 2]
      covariates <- model.matrix(formula, data)[, -1, drop = FALSE]  # Remove intercept
      list(time = time, status = status, covariates = covariates)
    },

  ## Select relevant columns
  select_columns = function(new_data) {
      if(!is.null(self$use_covariates)) {
          new_data_sel <- new_data %>% select(time, status, self$use_covariates)
      } else {
          new_data_sel <- new_data
      }
      return(new_data_sel)
  }
 
  ),
 
)


GRF_SurvivalForestWrapper <- R6::R6Class("GRF_SurvivalForestWrapper",
  inherit = SurvivalModelWrapper,
  public = list(

    # Fit method
    fit = function(formula, data) {
      data <- self$select_columns(data)
      parsed_data <- self$parse_formula(formula, data)

      # Fit the survival forest model
      self$model <- grf::survival_forest(parsed_data$covariates,
                                         Y = parsed_data$time,
                                         D = parsed_data$status,
                                         num.trees = 100)
      self$time.points <- self$model$failure.times  # Extract the default failure times
      if(length(self$time.points)==0) {
          self$time.points = seq(min(data$time), max(data$time), length.out=100)
      }
    },

    # Predict survival curves with optional custom failure times
    predict = function(new_data, time.points = NULL) {
      new_data <- self$select_columns(new_data)
      
      # Generate the design matrix from the new data using the stored formula
      covariates_new <- model.matrix(self$formula, new_data)[, -1, drop = FALSE]

      # Predict survival curves
      predictions <- predict(self$model, newdata = covariates_new)
      
      # Use default failure times if custom ones are not provided
      original_failure_times <- self$time.points
      survival_probs <- predictions$predictions

      # If predictions fail, replace survival_probs with ones
      if(length(survival_probs)==0){
          survival_probs = matrix(1, nrow = nrow(new_data), ncol = length(self$time.points))
      }

      # If custom failure times are provided, use the default interpolation method
      if (!is.null(time.points)) {
        survival_probs_interp <- self$predict_interpolate(survival_probs, original_failure_times, time.points)
        return(list(predictions = survival_probs_interp, time.points = time.points))
      }
      
      # If no custom failure times are provided, return the original predictions
      return(list(predictions = survival_probs, time.points = original_failure_times))
    }
  )
)


randomForestSRC_SurvivalWrapper <- R6::R6Class("randomForestSRC_SurvivalWrapper",
  inherit = SurvivalModelWrapper,
  public = list(

    # Constructor
    # Description: Default constructor
    initialize = function(use_covariates = NULL) {
      self$use_covariates <- use_covariates
    },

    # Fit method
    fit = function(formula, data, ntree = 100, ...) {
      data <- self$select_columns(data)

      ## Extract the model frame based on the formula and data
      mf <- model.frame(formula, data)

      # Ensure the response is a survival object (Surv)
      response <- model.response(mf)
      if (!inherits(response, "Surv")) {
        stop("The left-hand side of the formula must be a survival object (Surv(time, status)).")
      }
      # Fit the survival forest model using randomForestSRC
      self$model <- randomForestSRC::rfsrc(formula, data = data, ntree = ntree, save.memory=TRUE, ...)

      # Store the failure times from the model
      self$time.points <- self$model$time.interest

      if(length(self$time.points)==0) {
          self$time.points = seq(min(data$time), max(data$time), length.out=100)
      }
    },

    # Predict survival curves
    predict = function(new_data, time.points = NULL) { 
      new_data <- self$select_columns(new_data)
      # Ensure that new_data is correctly formatted
      if (!is.data.frame(new_data)) {
        stop("new_data must be a data frame.")
      }
      
      # Use the model's failure times if custom failure times are not provided
      if (is.null(time.points)) {
        time.points <- self$time.points
      }

      ## Predict survival curves using the trained random forest model
      predictions <- predict(self$model, newdata = new_data)
      ## Extract survival probabilities for each individual at each failure time
      survival_probs <- get_survival_prob_at_time(predictions, time.points)

      # If predictions fail, replace survival_probs with ones
      if(length(survival_probs)==0){
          survival_probs = matrix(1, nrow = nrow(new_data), ncol = length(self$time.points))
      }
      
        ## Return predictions and failure times
        list(predictions = survival_probs, time.points = time.points)
    }
  )
)


SurvregModelWrapper <- R6::R6Class("SurvregModelWrapper",
  inherit = SurvivalModelWrapper,
  public = list(
    dist = NULL,  # Distribution parameter

    # Constructor
    # Description: Initializes the SurvregModelWrapper object with the specified distribution.
    # Inputs:
    #   - dist: The distribution to be used in the survreg model (e.g., "weibull", "lognormal", etc.).
    initialize = function(use_covariates = NULL, dist = "weibull") {
      self$use_covariates <- use_covariates
      self$dist <- dist
    },

    # Fit method
    fit = function(formula, data) {
      data <- self$select_columns(data)
      parsed_data <- self$parse_formula(formula, data)
      self$model <- survival::survreg(formula, data = data, dist = self$dist, control = survival::survreg.control(maxiter = 1000))
      self$time.points <- seq(min(parsed_data$time), max(parsed_data$time), length.out = 100)
    },

    # Predict quantiles
    predict_quantiles = function(new_data, probs = c(0.25, 0.5, 0.75)) {

      # Predict quantiles for each probability in probs
      quantiles_matrix <- sapply(probs, function(p) {
        predict(self$model, newdata = new_data, type = "quantile", p = p)
      })      

      # Ensure quantiles_matrix is a matrix
      quantiles_matrix <- ensure_matrix(quantiles_matrix)

      # Set column names to represent the quantiles
      colnames(quantiles_matrix) <- paste0("Q", probs * 100, "%")
      rownames(quantiles_matrix) <- paste0("Individual_", 1:nrow(quantiles_matrix))

      # Return the quantile estimates as a data frame
      as.data.frame(quantiles_matrix)
    }
  )
)

CoxphModelWrapper <- R6::R6Class("CoxphModelWrapper",
  inherit = SurvivalModelWrapper,
  public = list(

    model = NULL,             # To store the fitted coxph model
    formula_env = NULL,        # To store the formula environment

    initialize = function(use_covariates = NULL) {
      self$use_covariates <- use_covariates
    },

    fit = function(formula, data) {
        data <- self$select_columns(data)
        # Fit the coxph model and store it
        parsed_data <- self$parse_formula(formula, data)
        self$model <- survival::coxph(formula, data = data, x = TRUE, y = TRUE)
        self$time.points <- seq(min(parsed_data$time), max(parsed_data$time), length.out = 100)

        # Capture and store the formula's environment
        #self$formula_env <- environment(formula)

        # Ensure the environment is set correctly for terms evaluation
        #environment(self$model$terms) <- self$formula_env
    },

    predict = function(new_data, time.points = NULL) {
        new_data <- self$select_columns(new_data)

        ## Use the model's failure times if custom failure times are not provided
        if (is.null(time.points)) {
            time.points <- self$time.points
        }

       
        ## Step 1: Calculate linear predictors for the new data
        linear_predictors <- predict(self$model, newdata = new_data, type = "lp")

        ## Step 2: Extract the baseline cumulative hazard
        baseline_hazard <- survival::basehaz(self$model, centered = TRUE)

       
        ## Step 3: Interpolate baseline hazard to match the requested time.points
        extended_hazard_times <- c(0, baseline_hazard$time, max(baseline_hazard$time) * 10)
        extended_hazard_probs <- c(0, baseline_hazard$hazard, 1)
        baseline_hazard_at_times <- approx(extended_hazard_times, extended_hazard_probs, xout = time.points)$y

        ## Step 4: Compute survival probabilities for all individuals
        ## Survival probability is calculated as: exp(-baseline hazard * exp(linear predictor))
        survival_probs <- t(exp(-outer(baseline_hazard_at_times, exp(linear_predictors), "*")))

        return(list(predictions = survival_probs, time.points = time.points))
    }


  )
)
