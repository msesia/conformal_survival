# Load required libraries
suppressMessages(library(tidyverse))

# Load useful scripts
source("../conf_surv/utils_data_new.R")
source("../conf_surv/utils_survival_new.R")
source("../conf_surv/utils_censoring.R")
source("../conf_surv/utils_conformal.R")
source("../conf_surv/utils_decensoring.R")

######################
## Input parameters ##
######################

num_samples_train <- 200
num_samples_cal <- 500

######################
## Fixed parameters ##
######################

# Sample sizes
num_samples_test <- 1000

## Nominal level
alpha <- 0.1


############################
# Define data distribution #
############################

# Initialize the covariate model
covariate_generator <- function(num_samples) {
    num_features = 10
    matrix(runif(num_samples * num_features, 0, 4), nrow = num_samples)
}

setting = 3

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

run.experiment <- function(random.state) {
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

random.state <- 2024
results <- run.experiment(random.state)
