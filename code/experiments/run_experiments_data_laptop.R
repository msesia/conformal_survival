suppressMessages(library(tidyverse))  # Load the tidyverse package without displaying messages

source("experiment_data_laptop.R")  # Source external script containing additional experiment-related functions

##################
## Data loading ##
##################

### Function to simulate a synthetic survival dataset ###
simulate_data <- function(num_samples_train,num_samples_cal, num_samples_test) {
    ## New setting (setting 8 in code, or setting 1 in paper, hard setting)
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

    ## Initialize the data generator
    generator <- SurvivalDataGenerator$new(covariate_generator, survival_generator, censoring_generator)


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

    return(list(train=data.train, cal=data.cal, test=data.test, test.oracle=data.test.oracle))
}

### Generate synthetic dataset with specified sample sizes ###
num_samples_train <- 200
num_samples_cal <- 200
num_samples_test <- 200

data <- simulate_data(num_samples_train, num_samples_cal, num_samples_test)

## Extract datasets for training, calibration, and testing
data.train <- data$train  # Training data
data.cal <- data$cal      # Calibration data
data.test <- data$test    # Test data

## Test data with true event times (used for evaluation, i.e., computing coverage)
data.test.oracle <- data$test.oracle

#######################
## Survival analysis ##
#######################

### Model configuration ###
surv_model_type <- "grf"  # Survival model type (e.g., "grf" for Generalized Random Forests)
cens_model_type <- "cox"  # Censoring model type (e.g., "cox" for Cox proportional hazards)

## Create directories for results
output_dir <- "results_data"  # Directory to store analysis results
dir.create("results", showWarnings=FALSE)
dir.create(sprintf("results/%s", output_dir), showWarnings=FALSE)

### Run survival analysis for different confidence levels (alpha) ###
for(alpha in c(0.1,0.2)) {
    ## Generate a unique file name based on model types and alpha value
    output_file <- paste0("results/", output_dir, "/",
                          "surv_", surv_model_type,
                          "_cens_", cens_model_type,
                          "_alpha", alpha,
                          ".txt")

    ## Store important parameters including model types
    header <- tibble(surv_model_type = surv_model_type,
                     cens_model_type = cens_model_type,
                     n_train = num_samples_train,
                     n_cal = num_samples_cal,
                     n_test = num_samples_test,
                     alpha = alpha)

    ## Run survival analysis and save results to file    
    predictions <- run_analysis(data.train, data.cal, data.test, surv_model_type, cens_model_type, alpha)

    ## Evaluate results
    ## Note: coverage can be evaluated only if 'oracle' test data with true event times is available
    results <- do.call(rbind, lapply(names(predictions), function(name) {
        if(is.null(data.test.oracle)) {
            res = evaluate_bounds(data.test$time, predictions[[name]])
        } else {
            res = evaluate_bounds(data.test$time, predictions[[name]], event_time=data.test.oracle$event_time)
        }

        cbind(Method = name, res)
    })) |> cbind(header)

    ## Write the results to the CSV file
    write.csv(results, output_file, row.names = FALSE)

    ## Print the output file name to verify
    cat("Results saved ub file file:", output_file, "\n")
}

### Function to load analysis results from the specified output directory ###
load_data <- function(output_dir) {
  idir <- sprintf("results/%s", output_dir)
  ifile.list <- list.files(idir)
  results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
  }))
  return(results)
}

## Load results from saved files
results <- load_data(output_dir)
