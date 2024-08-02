suppressMessages(library(tidyverse))

## Function to generate covariates X
generate_covariates <- function(num_samples, num_features) {
    ## Generate i.i.d. uniform random variables on [-1,1]
    X <- matrix(runif(num_samples * num_features, -1, 1), nrow = num_samples)
    return(X)
}

## Function to generate true survival times T|X
generate_true_survival <- function(X) {
    num_samples <- nrow(X)
    ## Log-normal model (Mvt + Heterosc.) from Candes et al. (Table 1)
    mu_x <- log(2) + 1 + 0.55*(X[, 1]^2 - X[,3]*X[,5])
    sigma_x <- 0.1 * (abs(X[,10]) + 1)
    log_T <- rnorm(num_samples, mean=mu_x, sd=sigma_x)
    T <- exp(log_T)
    return(T)
}

## Function to generate censoring times C|X
generate_censoring_times <- function(X) {
    num_samples <- nrow(X)
    censor_rate = 0.2
    C <- rexp(num_samples, rate = censor_rate)
    return(C)
}

## Function to combine and create observed data
generate_survival_data <- function(num_samples, num_features) {
    X <- generate_covariates(num_samples, num_features)
    T <- generate_true_survival(X)
    C <- generate_censoring_times(X)
    observed_time <- pmin(T, C)
    event <- observed_time == T
    data <- data.frame(event_time=T, censoring_time=C, observed_time, event, X)
    return(data)
}
