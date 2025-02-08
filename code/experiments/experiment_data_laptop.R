## Load required libraries
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

run_analysis <- function(data.train, data.cal, data.test, surv_model_type, cens_model_type, alpha) {

    ######################
    ## Fixed parameters ##
    ######################

    ## Whether to use (approximate) finite sample correction for Gui's method
    fsc <- FALSE

    ###################################################
    ## Instantiate the survival and censoring models ##
    ###################################################

    ## Define a mapping between model types and constructors for survival models
    model_constructors_surv <- list(
        grf = GRF_SurvivalForestWrapper$new,
        survreg = function() SurvregModelWrapper$new(dist="lognormal"),
        rf = randomForestSRC_SurvivalWrapper$new,
        cox = CoxphModelWrapper$new
    )

    ## Instantiate survival model based on the specified type
    if (!is.null(model_constructors_surv[[surv_model_type]])) {
        surv_model <- model_constructors_surv[[surv_model_type]]()
    } else {
        stop("Unknown survival model type!")
    }

    ## Define a mapping between model types and constructors for censoring models
    model_constructors_cens <- list(
        grf = function() GRF_SurvivalForestWrapper$new(),
        survreg = function() SurvregModelWrapper$new(, dist="lognormal"),
        rf = function() randomForestSRC_SurvivalWrapper$new(),
        cox = function() CoxphModelWrapper$new()
    )

    ## Instantiate censoring model based on the specified type
    if (!is.null(model_constructors_cens[[cens_model_type]])) {
        cens_base_model <- model_constructors_cens[[cens_model_type]]()
    } else {
        stop("Unknown censoring model type!")
    }

    ## Create an instance of the CensoringModel class with the model
    cens_model <- CensoringModel$new(model = cens_base_model)

    #######################################
    ## Define function to analyze the data #
    #######################################

    analyze_data <- function(data.train, data.cal, data.test, surv_model, cens_model) {
        ## Fit the survival model on the training data
        surv_model$fit(Surv(time, status) ~ ., data = data.train)

        ## Fit the censoring model on the training data
        cens_model$fit(data = data.train)

        ## Fit the Kaplan-Meier survival model for the de-censoring method
        surv_object <- Surv(time = data.train$time, event = data.train$status)
        km_fit <- survival::survfit(surv_object ~ 1)

        ## Apply all methods
        predictions <- c()

        ## Construct uncalibrated lower bound
        predictions$uncalibrated <- surv_model$predict_quantiles(data.test, probs=alpha)[[1]]

        ## Apply naive CQR
        predictions$cqr <- predict_CQR(data.test, surv_model, data.cal, alpha)

        ## Apply CQR with de-censoring
        predictions$cqr.decensor <- predict_decensoring(data.test, surv_model, km_fit, data.cal, alpha, R=10)

        ## Apply drcosarc (Candes, with fixed c0)
        predictions$drcosarc.candes <- predict_drcosarc(data.test, surv_model, cens_model, data.cal, alpha, cutoffs="candes-fixed")

        ## Apply drcosarc (Gui)
        predictions$drcosarc.gui <- predict_drcosarc(data.test, surv_model, cens_model, data.cal, alpha, cutoffs="adaptive", finite_sample_correction=fsc)

        return(predictions)
    }

    #######################################
    ## Define function to run experiment ##
    #######################################

    ## Compute prediction lower bounds
    predictions <- analyze_data(data.train, data.cal, data.test, surv_model, cens_model)
    
    return(predictions)
}
