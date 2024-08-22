# Load required libraries
suppressMessages(library(tidyverse))

## Source utility functions for data generation and analysis
source("../conf_surv/utils_data_new.R")
source("../conf_surv/utils_survival_new.R")
source("../conf_surv/utils_censoring.R")
source("../conf_surv/utils_conformal.R")
source("../conf_surv/utils_decensoring.R")

######################
## Input parameters ##
######################

## Flag to determine if input should be parsed from command line
parse_input <- TRUE

if(parse_input) {
    ## Reading command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    ## Checking if the correct number of arguments is provided
    if (length(args) < 5) {
        stop("Insufficient arguments provided. Expected 4 arguments.")
    }
    ## Assigning command line arguments to variables
    setup <- as.integer(args[1])
    setting <- as.integer(args[2])
    num_samples_train <- as.integer(args[3])
    num_samples_cal <- as.integer(args[4])
    batch <- as.integer(args[5])
} else {
    setup <- 1
    setting <- 3
    num_samples_train <- 200
    num_samples_cal <- 500
    batch <- 1
}

######################
## Fixed parameters ##
######################

# Sample sizes
num_samples_test <- 1000

## Nominal level
alpha <- 0.1

## Number of repetitions
batch_size <- 10

####################
## Prepare output ##
####################

## Store important parameters
header <- tibble(setup = setup, setting = setting, n_train = num_samples_train, n_cal = num_samples_cal, alpha = alpha, batch = batch)

## Generate a unique and interpretable file name based on the input parameters
output_file <- paste0("results/setup_", setup, "/", "setting", setting, "_train", num_samples_train, "_cal", num_samples_cal, "_batch", batch, ".txt")

## Print the output file name to verify
cat("Output file name:", output_file, "\n")

############################
# Define data distribution #
############################

# Initialize the covariate model
covariate_generator <- function(num_samples) {
    num_features = 10
    matrix(runif(num_samples * num_features, 0, 4), nrow = num_samples)
}

if(setting==3) {
    # Initialize the survival time distribution
    surv_mu_fun <- function(X) 0*X[,1] + 2 * (X[,1]>2) + X[,1] * (X[,1]<2)
    surv_sigma_fun <- function(X) 0*X[,1] + 0.5
    survival_generator <- LogNormalDistribution$new(mu_fun = surv_mu_fun, sigma_fun = surv_sigma_fun)
    # Initialize the censoring time distribution
    cens_rate_fun <- function(X) 0*X[,1] + 0.25 + (6+X[,1])/100
    censoring_generator <- ExponentialDistribution$new(rate_fun = cens_rate_fun)
} else if(setting==4) {
    # Initialize the survival time distribution
    surv_mu_fun <- function(X) 0*X[,1] + 3 * (X[,1]>2) + 1.5 * X[,1] * (X[,1]<2)
    surv_sigma_fun <- function(X) 0*X[,1] + 0.5
    survival_generator <- LogNormalDistribution$new(mu_fun = surv_mu_fun, sigma_fun = surv_sigma_fun)
    # Initialize the censoring time distribution
    cens_mu_fun <- function(X) 0*X[,1] + 2 + (2-X[,1])/50
    cens_sigma_fun <- function(X) 0*X[,1] + 0.5
    censoring_generator <- LogNormalDistribution$new(mu_fun = cens_mu_fun, sigma_fun = cens_sigma_fun)
}

# Initialize the data generator
generator <- SurvivalDataGenerator$new(covariate_generator, survival_generator, censoring_generator)

##############################################
## Define the survival and censoring models ##
##############################################

surv_model <- CoxphModelWrapper$new()
cens_base_model <- CoxphModelWrapper$new()

# Create an instance of the CensoringModel class with the model
cens_model <- CensoringModel$new(model = cens_base_model)

#######################################
# Define function to analyze the data #
#######################################

analyze_data <- function(data.train, data.cal, data.test, surv_model, cens_model, generator=NULL, C.cal.oracle=NULL) {

    # Fit the survival model on the training data
    surv_model$fit(survival::Surv(time, status) ~ ., data = data.train)

    # Fit the censoring model on the training data
    cens_model$fit(data = data.train)

    ## Fit the Kaplan-Meier survival model for the de-censoring method
    surv_object <- survival::Surv(time = data.train$time, event = data.train$status)
    km_fit <- survival::survfit(surv_object ~ 1)

    # Apply all methods
    predictions <- c()

    # Construct oracle lower bound
    if(!is.null(generator)) {
        predictions$oracle <- generator$survival$predict_quantiles(select(data.test, -time, -status), probs=c(alpha))[[1]]
    }

    # Construct nominal lower bound
    predictions$nominal <- surv_model$predict_quantiles(data.test, probs=alpha)[[1]]

    # Apply naive CQR
    predictions$cqr <- predict_CQR(data.test, surv_model, data.cal, alpha)

    # Apply CQR with de-censoring
    predictions$cqr.decensor <- predict_decensoring(data.test, surv_model, km_fit, data.cal, alpha, R=10)

    if(!is.null(C.cal.oracle)) {
        ## Apply Candes' method with "oracle" censoring model
        predictions$candes.oracle <- predict_Candes(data.test, surv_model, generator$censoring, data.cal, C.cal.oracle, alpha)

        ## Apply Gui's method with "oracle" censoring model
        predictions$gui.oracle <- predict_Gui(data.test, surv_model, generator$censoring, data.cal, C.cal.oracle, alpha)
    }

    # Apply prototype (Candes)
    predictions$prototype.candes <- predict_prototype(data.test, surv_model, cens_model, data.cal, alpha)

    # Apply prototype (Gui)
    predictions$prototype.gui <- predict_prototype(data.test, surv_model, cens_model, data.cal, alpha, cutoffs="adaptive")


    return(predictions)
}

#######################################
## Define function to run experiment ##
#######################################

run_experiment <- function(random.state) {
    set.seed(random.state)

    ## Generate training, calibration, and test data
    ## (including true event and censoring times)
    data.train.oracle <- generator$sample(num_samples_train)
    data.cal.oracle <- generator$sample(num_samples_cal)
    data.test.oracle <- generator$sample(num_samples_test)
    C.cal.oracle <- data.cal.oracle$censoring_time

    ## Remove true event and censoring times from the data (right-censoring)
    data.train <- data.train.oracle |> select(-event_time, -censoring_time)
    data.cal <- data.cal.oracle |> select(-event_time, -censoring_time)
    data.test <- data.test.oracle |> select(-event_time, -censoring_time)

    ## Compute prediction lower bounds
    predictions <- analyze_data(data.train, data.cal, data.test, surv_model, cens_model, generator=generator, C.cal.oracle=C.cal.oracle)

    ## Evaluate results
    results <- do.call(rbind, lapply(names(predictions), function(name) {
        res = evaluate_bounds(data.test$time,
                              predictions[[name]],
                              event_time=data.test.oracle$event_time)

        cbind(Method = name, res)
    }))

    return(results)
}


## Function to run multiple experiments and gather results
## Args:
##   batch_size: Number of repetitions for each experimental setting
## Returns:
##   A tibble containing the combined results of all experiments
run_multiple_experiments <- function(batch_size) {
    results_df <- data.frame()  # Initialize an empty data frame to store cumulative results

    # Print a progress bar header
    cat("Running experiments\n")
    pb <- txtProgressBar(min = 0, max = batch_size, style = 3)  # Initialize progress bar

    # Loop over each repetition
    for (i in 1:batch_size) {
        random.state <- batch*1000 + i
        res <- run_experiment(random.state)  # Run experiment and get the result

        ## Combine the results with experiment metadata
        result_df <- tibble(Seed = random.state) |> cbind(header) |> cbind(res)

        # Add the result to the cumulative data frame
        results_df <- rbind(results_df, result_df)

        # Write the cumulative results to the CSV file
        write.csv(results_df, output_file, row.names = FALSE)

        setTxtProgressBar(pb, i)  # Update progress bar
    }

    close(pb)  # Close the progress bar

    return(results_df)  # Return the cumulative results data frame
}

#####################
## Run experiments ##
#####################

## Run the experiments with specified parameters
results <- run_multiple_experiments(batch_size)
