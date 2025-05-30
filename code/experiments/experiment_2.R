# Load required libraries
suppressMessages(library(tidyverse))
suppressMessages(library(survival))

## Source utility functions for data generation and analysis
source("utils_data.R")
source("../conf_surv/utils_survival.R")
source("../conf_surv/utils_censoring.R")
source("../conf_surv/utils_conformal.R")
source("../conf_surv/utils_decensoring.R")

######################
## Input parameters ##
######################

# Available datasets:
# VALCT: Lung cancer trial data (137 obs, 6 vars, from survival::veteran).
# PBC: Liver disease data (137 obs, 6 vars, from survival::pbc, NAs imputed).
# GBSG: Breast cancer data (2232 obs, 7 vars, combined train/test).
# METABRIC: Breast cancer clinical data (1981 obs, 25 vars).

## Flag to determine if input should be parsed from command line
parse_input <- TRUE

if(parse_input) {
    ## Reading command line arguments
    args <- commandArgs(trailingOnly = TRUE)

    ## Checking if the correct number of arguments is provided
    if (length(args) < 6) {
        stop("Insufficient arguments provided. Expected 4 arguments.")
    }

    ## Assigning command line arguments to variables
    dataset <- args[1]
    surv_model_type <- args[2]
    cens_model_type <- args[3]
    train_prop_sub <- as.double(args[4])
    alpha <- as.double(args[5])
    batch <- as.integer(args[6])

} else {
    dataset <- "METABRIC"
    surv_model_type <- "grf"
    cens_model_type <- "grf"
    train_prop_sub = 1
    alpha <- 0.1
    batch <- 1
}

######################
## Fixed parameters ##
######################

## Relative calibration sample size
cal_prop = 0.2

## Relative test sample size
test_prop = 0.2

## Relative training sample size
train_prop = 1 - (cal_prop + test_prop)

## Number of repetitions
batch_size <- 10

## Whether to use (approximate) finite sample correction for Gui's method
fsc <- FALSE

####################
## Prepare output ##
####################

## Store important parameters including model types
header <- tibble(dataset = dataset,
                 surv_model_type = surv_model_type,
                 cens_model_type = cens_model_type,
                 train_prop_sub = train_prop_sub,
                 alpha = alpha,
                 batch = batch)

## Generate a unique and interpretable file name based on the input parameters
output_file <- paste0("results/data/", dataset,
                      "_surv_", surv_model_type,
                      "_cens_", cens_model_type,
                      "_train_", train_prop_sub,
                      "_alpha_", alpha,
                      "_batch", batch, ".txt")

## Print the output file name to verify
cat("Output file name:", output_file, "\n")


############################
# Define data distribution #
############################

data <- load_csv(dataset)

## Data features
num_features <- ncol(data) - 2

# Rename the columns
col.names <- c("time", "status", paste("X", 1:num_features, sep = ""))
colnames(data) <- col.names

## Use all features
num_feat_censor <- num_features

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
    num_samples_train <- nrow(data.train)
    num_samples_train_cens <- num_samples_train
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

    ## Apply drcosarc (Candes, with fixed c0)
    predictions$drcosarc.candes <- predict_drcosarc(data.test, surv_model, cens_model, data.cal, alpha, cutoffs="candes-fixed")

    ## ## ## Apply drcosarc (Candes, with tuned c0)
    ## tuning.package <- list(data.train = data.train)
    ## predictions$drcosarc.candes.tuned <- predict_drcosarc(data.test, surv_model, cens_model, data.cal, alpha, tuning.package=tuning.package, cutoffs="candes-tuning")

    ## Apply drcosarc (Gui)
    predictions$drcosarc.gui <- predict_drcosarc(data.test, surv_model, cens_model, data.cal, alpha, cutoffs="adaptive", finite_sample_correction=fsc)

    return(predictions)
}

#######################################
## Define function to run experiment ##
#######################################

run_experiment <- function(random.state) {
    set.seed(random.state)

    ## Split the dataset
    split_result <- split_data(data, train_prop = train_prop, cal_prop = cal_prop, test_prop = test_prop, seed = random.state)

    ## Access subsets
    data.train <- split_result$train
    data.cal <- split_result$cal
    data.test <- split_result$test

    ## Subsample the training set
    data.train = data.train[sample(1:nrow(data.train), ceiling(train_prop_sub*nrow(data.train))),]
    
    ## Compute prediction lower bounds
    predictions <- analyze_data(data.train, data.cal, data.test, surv_model, cens_model)

    ## Evaluate results
    results <- do.call(rbind, lapply(names(predictions), function(name) {
        res = evaluate_bounds(data.test$time, predictions[[name]],
                              event_time=data.test$event_time, status=data.test$status)

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
