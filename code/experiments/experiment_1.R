# Load required libraries
suppressMessages(library(tidyverse))
library(survival)

## Source utility functions for data generation and analysis
source("../conf_surv/utils_data.R")
source("../conf_surv/utils_survival.R")
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
    if (length(args) < 9) {
        stop("Insufficient arguments provided. Expected 7 arguments.")
    }
    
    ## Assigning command line arguments to variables
    setup <- as.integer(args[1])
    setting <- as.integer(args[2])
    surv_model_type <- args[3]
    cens_model_type <- args[4]
    num_samples_train <- as.integer(args[5])
    num_samples_train_cens <- as.integer(args[6])
    num_samples_cal <- as.integer(args[7])
    num_feat_censor <- as.integer(args[8])
    batch <- as.integer(args[9])

} else {
    setup <- 1
    setting <- 3
    surv_model_type <- "grf"
    cens_model_type <- "cox"
    num_samples_train <- 200
    num_samples_cal <- 500
    num_feat_censor <- 10
    batch <- 1
}

######################
## Fixed parameters ##
######################

## Test sample size
num_samples_test <- 1000

## Nominal level
alpha <- 0.1

## Number of repetitions
batch_size <- 10

## Whether to use (approximate) finite sample correction for Gui's method
fsc <- FALSE

####################
## Prepare output ##
####################

## Store important parameters including model types
header <- tibble(setup = setup, 
                 setting = setting, 
                 surv_model_type = surv_model_type, 
                 cens_model_type = cens_model_type,
                 n_train = num_samples_train, 
                 n_train_cens = num_samples_train_cens, 
                 n_cal = num_samples_cal,
                 num_feat_censor = num_feat_censor,                 
                 alpha = alpha, 
                 batch = batch)

## Generate a unique and interpretable file name based on the input parameters
output_file <- paste0("results/setup_", setup, "/", 
                      "setting", setting, 
                      "_surv_", surv_model_type, 
                      "_cens_", cens_model_type,
                      "_train", num_samples_train, 
                      "_trainc", num_samples_train_cens, 
                      "_cal", num_samples_cal, 
                      "_nfc", num_feat_censor,
                      "_batch", batch, ".txt")

## Print the output file name to verify
cat("Output file name:", output_file, "\n")


############################
# Define data distribution #
############################

if(setting==1) {
    ## Setting 1 from Gui et al.
    ## Initialize the covariate model
    num_features <- 1
    covariate_generator <- function(num_samples) {
        matrix(runif(num_samples * num_features, 0, 4), nrow = num_samples)
    }
    ## Initialize the survival time distribution
    surv_mu_fun <- function(X) 0*X[,1] + 0.632 * X[,1]
    surv_sigma_fun <- function(X) 0*X[,1] + 2
    survival_generator <- LogNormalDistribution$new(mu_fun = surv_mu_fun, sigma_fun = surv_sigma_fun)
    ## Initialize the censoring time distribution
    cens_rate_fun <- function(X) 0*X[,1] + 0.1
    censoring_generator <- ExponentialDistribution$new(rate_fun = cens_rate_fun)

} else if(setting==2) {
    ## Setting 2 from Gui et al.
    ## Initialize the covariate model
    num_features <- 1
    covariate_generator <- function(num_samples) {
        matrix(runif(num_samples * num_features, 0, 4), nrow = num_samples)
    }
    ## Initialize the survival time distribution
    surv_mu_fun <- function(X) 0*X[,1] + 3 * (X[,1]>2) + X[,1] * (X[,1]<2)
    surv_sigma_fun <- function(X) 0*X[,1] + 0.5
    survival_generator <- LogNormalDistribution$new(mu_fun = surv_mu_fun, sigma_fun = surv_sigma_fun)
    ## Initialize the censoring time distribution
    cens_rate_fun <- function(X) 0*X[,1] + 0.1
    censoring_generator <- ExponentialDistribution$new(rate_fun = cens_rate_fun)

} else if(setting==3) {
    ## Setting 3 from Gui et al.
    ## Initialize the covariate model
    num_features <- 1
    covariate_generator <- function(num_samples) {
        matrix(runif(num_samples * num_features, 0, 4), nrow = num_samples)
    }
    ## Initialize the survival time distribution
    surv_mu_fun <- function(X) 0*X[,1] + 2 * (X[,1]>2) + X[,1] * (X[,1]<2)
    surv_sigma_fun <- function(X) 0*X[,1] + 0.5
    survival_generator <- LogNormalDistribution$new(mu_fun = surv_mu_fun, sigma_fun = surv_sigma_fun)
    ## Initialize the censoring time distribution
    cens_rate_fun <- function(X) 0*X[,1] + 0.25 + (X[,1]+6)/100
    censoring_generator <- ExponentialDistribution$new(rate_fun = cens_rate_fun)

    
} else if(setting==4) {
    ## Setting 4 from Gui et al.
    ## Initialize the covariate model
    num_features <- 1
    covariate_generator <- function(num_samples) {
        matrix(runif(num_samples * num_features, 0, 4), nrow = num_samples)
    }
    ## Initialize the survival time distribution
    surv_mu_fun <- function(X) 0*X[,1] + 3 * (X[,1]>2) + 1.5 * X[,1] * (X[,1]<2)
    surv_sigma_fun <- function(X) 0*X[,1] + 0.5
    survival_generator <- LogNormalDistribution$new(mu_fun = surv_mu_fun, sigma_fun = surv_sigma_fun)
    ## Initialize the censoring time distribution
    cens_mu_fun <- function(X) 0*X[,1] + 2 + (2-X[,1])/50
    cens_sigma_fun <- function(X) 0*X[,1] + 0.5
    censoring_generator <- LogNormalDistribution$new(mu_fun = cens_mu_fun, sigma_fun = cens_sigma_fun)
    
} else if(setting==5) {
    ## Setting 5 from Gui et al.
    ## Initialize the covariate model
    num_features <- 10
    covariate_generator <- function(num_samples) {
        matrix(runif(num_samples * num_features, 0, 4), nrow = num_samples)
    }
    ## Initialize the survival time distribution
    surv_mu_fun <- function(X) 0*X[,1] + 0.126 * (X[,1] + sqrt(X[,3]*X[,5])) + 1
    surv_sigma_fun <- function(X) 0*X[,1] + 0.5
    survival_generator <- LogNormalDistribution$new(mu_fun = surv_mu_fun, sigma_fun = surv_sigma_fun)
    ## Initialize the censoring time distribution
    cens_rate_fun <- function(X) 0*X[,1] + X[,10]/10 + 1/20
    censoring_generator <- ExponentialDistribution$new(rate_fun = cens_rate_fun)

} else if(setting==6) {
    ## Setting 6 from Gui et al.
    ## Initialize the covariate model
    num_features <- 10
    covariate_generator <- function(num_samples) {
        matrix(runif(num_samples * num_features, 0, 4), nrow = num_samples)
    }
    ## Initialize the survival time distribution
    surv_mu_fun <- function(X) 0*X[,1] + 0.126 * (X[,1] + sqrt(X[,3]*X[,5])) + 1
    surv_sigma_fun <- function(X) 0*X[,1] + (X[,2]+2)/4
    survival_generator <- LogNormalDistribution$new(mu_fun = surv_mu_fun, sigma_fun = surv_sigma_fun)
    ## Initialize the censoring time distribution
    cens_rate_fun <- function(X) 0*X[,1] + X[,10]/10 + 1/20
    censoring_generator <- ExponentialDistribution$new(rate_fun = cens_rate_fun)
    
} else if(setting==7) {
    ## New setting
    ## Initialize the covariate model
    num_features <- 100
    covariate_generator <- function(num_samples) {
        matrix(runif(num_samples * num_features, 0, 1), nrow = num_samples)
    }
    ## Initialize the survival time distribution
    surv_mu_fun <- function(X) 0*X[,1] + X[,1]^0.25
    surv_sigma_fun <- function(X) 0*X[,1] + 0.1
    survival_generator <- LogNormalDistribution$new(mu_fun = surv_mu_fun, sigma_fun = surv_sigma_fun)
    ## Initialize the censoring time distribution
    cens_mu_fun <- function(X) 0*X[,1] + X[,1]^4 + 0.4
    cens_sigma_fun <- function(X) 0*X[,1] + 0.1
    censoring_generator <- LogNormalDistribution$new(mu_fun = cens_mu_fun, sigma_fun = cens_sigma_fun)

} else if(setting==8) {
    ## New setting
    ## Initialize the covariate model
    num_features <- 100
    covariate_generator <- function(num_samples) {
        matrix(runif(num_samples * num_features, 0, 1), nrow = num_samples)
    }
    ## Initialize the survival time distribution
    surv_mu_fun <- function(X) 0*X[,1] + (X[,2]>0.5) + (X[,3]<0.5) + (1-X[,1])^0.25
    surv_sigma_fun <- function(X) 0*X[,1] + 0.1 * (1-X[,1])
    survival_generator <- LogNormalDistribution$new(mu_fun = surv_mu_fun, sigma_fun = surv_sigma_fun)
    ## Initialize the censoring time distribution
    cens_mu_fun <- function(X) 0*X[,1] + (X[,2]>0.5) + (X[,3]<0.5) + (1-X[,1])^4 + 0.4
    cens_sigma_fun <- function(X) 0*X[,1] + 0.1 * X[,2]
    censoring_generator <- LogNormalDistribution$new(mu_fun = cens_mu_fun, sigma_fun = cens_sigma_fun)

} else if(setting==9) {
    ## Setting 3 from Candes et al.
    ## Initialize the covariate model
    num_features <- 100
    covariate_generator <- function(num_samples) {
        matrix(runif(num_samples * num_features, -1, 1), nrow = num_samples)
    }
    ## Initialize the survival time distribution
    surv_mu_fun <- function(X) 0*X[,1] + log(2) + 1 + 0.55*(X[,1]^2-X[,3]*X[,5])
    surv_sigma_fun <- function(X) 0*X[,1] + 1
    survival_generator <- LogNormalDistribution$new(mu_fun = surv_mu_fun, sigma_fun = surv_sigma_fun)
    ## Initialize the censoring time distribution
    cens_rate_fun <- function(X) 0*X[,1] + 0.4
    censoring_generator <- ExponentialDistribution$new(rate_fun = cens_rate_fun)

} else if(setting==10) {
    ## Setting 4 from Candes et al.
    ## Initialize the covariate model
    num_features <- 100
    covariate_generator <- function(num_samples) {
        matrix(runif(num_samples * num_features, -1, 1), nrow = num_samples)
    }
    ## Initialize the survival time distribution
    surv_mu_fun <- function(X) 0*X[,1] + log(2) + 1 + 0.55*(X[,1]^2-X[,3]*X[,5])
    surv_sigma_fun <- function(X) 0*X[,1] + abs(X[,10]) + 1
    survival_generator <- LogNormalDistribution$new(mu_fun = surv_mu_fun, sigma_fun = surv_sigma_fun)
    ## Initialize the censoring time distribution
    cens_rate_fun <- function(X) 0*X[,1] + 0.4
    censoring_generator <- ExponentialDistribution$new(rate_fun = cens_rate_fun)
}


# Initialize the data generator
generator <- SurvivalDataGenerator$new(covariate_generator, survival_generator, censoring_generator)

###################################################
## Instantiate the survival and censoring models ##
###################################################

# Define a mapping between model types and constructors for survival models
model_constructors_surv <- list(
    grf = GRF_SurvivalForestWrapper$new,
    survreg = function() SurvregModelWrapper$new(dist="lognormal"),
    rf = randomForestSRC_SurvivalWrapper$new,
    cox = CoxphModelWrapper$new
)

# Instantiate survival model based on the specified type
if (!is.null(model_constructors_surv[[surv_model_type]])) {
    surv_model <- model_constructors_surv[[surv_model_type]]()
} else {
    stop("Unknown survival model type!")
}


# List of covariates to use for censoring model
use.covariates <- paste("X", 1:min(num_features, num_feat_censor), sep="")

# Define a mapping between model types and constructors for censoring models
model_constructors_cens <- list(
    grf = function() GRF_SurvivalForestWrapper$new(use_covariates=use.covariates),
    survreg = function() SurvregModelWrapper$new(use_covariates=use.covariates, dist="lognormal"),
    rf = function() randomForestSRC_SurvivalWrapper$new(use_covariates=use.covariates),
    cox = function() CoxphModelWrapper$new(use_covariates=use.covariates)
)

# Instantiate censoring model based on the specified type
if (!is.null(model_constructors_cens[[cens_model_type]])) {
    cens_base_model <- model_constructors_cens[[cens_model_type]]()
} else {
    stop("Unknown censoring model type!")
}

# Create an instance of the CensoringModel class with the model
cens_model <- CensoringModel$new(model = cens_base_model)

#######################################
# Define function to analyze the data #
#######################################

analyze_data <- function(data.train, data.cal, data.test, surv_model, cens_model, generator=NULL, C.train.oracle=NULL, C.cal.oracle=NULL) {

    ## Fit the survival model on the training data
    surv_model$fit(Surv(time, status) ~ ., data = data.train)
    
    ## Fit the censoring model on a subset of the training data
    if(num_samples_train_cens < num_samples_train) {
        idx.train.cens <- sort(sample(1:nrow(data.train), num_samples_train_cens))
    } else {
        idx.train.cens <- 1:nrow(data.train)
    }
    cens_model$fit(data = data.train[idx.train.cens,])

    ## Fit the Kaplan-Meier survival model for the de-censoring method
    surv_object <- Surv(time = data.train$time, event = data.train$status)
    km_fit <- survival::survfit(surv_object ~ 1)

    ## Apply all methods
    predictions <- c()

    ## Construct oracle lower bound
    if(!is.null(generator)) {
        predictions$oracle <- generator$survival$predict_quantiles(select(data.test, -time, -status), probs=c(alpha))[[1]]
    }

    ## Construct nominal lower bound
    predictions$nominal <- surv_model$predict_quantiles(data.test, probs=alpha)[[1]]

    ## Apply naive CQR
    predictions$cqr <- predict_CQR(data.test, surv_model, data.cal, alpha)

    ## Apply CQR with de-censoring
    predictions$cqr.decensor <- predict_decensoring(data.test, surv_model, km_fit, data.cal, alpha, R=10)

    if(!is.null(C.cal.oracle)) {
        ## Apply Candes' method with "oracle" censoring model (with fixed c0)
        predictions$candes.oracle <- predict_Candes(data.test, surv_model, generator$censoring, data.cal, C.cal.oracle, alpha)

        ## ## Apply Candes' method with "oracle" censoring model (with fixed c0)
        ## if(!is.null(C.train.oracle)) {
        ##     tuning.package.oracle <- list(data.train = data.train, C.train = C.train.oracle)
        ##     predictions$candes.oracle.tuned <- predict_Candes(data.test, surv_model, generator$censoring, data.cal, C.cal.oracle, alpha,
        ##                                                       tuning.package=tuning.package.oracle)
        ## }

        ## Apply Gui's method with "oracle" censoring model
        predictions$gui.oracle <- predict_Gui(data.test, surv_model, generator$censoring, data.cal, C.cal.oracle, alpha, use_cqr=FALSE, use_censoring_model=TRUE,
                                              finite_sample_correction=fsc)
        ## Apply Gui's method with "oracle" censoring model, with CQR approach
        predictions$gui.oracle.cqr <- predict_Gui(data.test, surv_model, generator$censoring, data.cal, C.cal.oracle, alpha, use_cqr=TRUE, finite_sample_correction=fsc)
    }

    ## Apply prototype (Candes, with fixed c0)
    predictions$prototype.candes <- predict_prototype(data.test, surv_model, cens_model, data.cal, alpha, cutoffs="candes-fixed")

    ## ## ## Apply prototype (Candes, with tuned c0)
    ## tuning.package <- list(data.train = data.train)
    ## predictions$prototype.candes.tuned <- predict_prototype(data.test, surv_model, cens_model, data.cal, alpha, tuning.package=tuning.package, cutoffs="candes-tuning")

    ## Apply prototype (Gui)
    predictions$prototype.gui <- predict_prototype(data.test, surv_model, cens_model, data.cal, alpha, cutoffs="adaptive", finite_sample_correction=fsc)

    ## Apply prototype (Gui, CQR)
    predictions$prototype.gui.cqr <- predict_prototype(data.test, surv_model, cens_model, data.cal, alpha, cutoffs="adaptive-cqr")   

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
    C.train.oracle <- data.train.oracle$censoring_time
    C.cal.oracle <- data.cal.oracle$censoring_time

    ## Remove true event and censoring times from the data (right-censoring)
    data.train <- data.train.oracle |> select(-event_time, -censoring_time)
    data.cal <- data.cal.oracle |> select(-event_time, -censoring_time)
    data.test <- data.test.oracle |> select(-event_time, -censoring_time)

    ## Compute prediction lower bounds
    predictions <- analyze_data(data.train, data.cal, data.test, surv_model, cens_model, generator=generator,
                                C.train.oracle=C.train.oracle, C.cal.oracle=C.cal.oracle)

    ## Evaluate results
    results <- do.call(rbind, lapply(names(predictions), function(name) {
        res = evaluate_bounds(data.test$time, predictions[[name]],
                              event_time=data.test.oracle$event_time, oracle=predictions$oracle)

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
