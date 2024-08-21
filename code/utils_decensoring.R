predict_decensoring <- function(data.test, surv_model, km_fit, data.cal, alpha, R=1) {
    ## From Qi et al.

    impute_kaplan_meier <- function(data.cal, km_fit, R) {
      # Initialize an empty list to store the results
      imputed_data_list <- vector("list", R)

      # Repeat the imputation process R times
      for (r in 1:R) {
        # Copy the original data
        data.cal.imputed <- data.cal

        # Extract the original survival times
        Y.cal <- data.cal$time

        # Find the indices of censored observations
        idx.censored <- which(data.cal$status == FALSE)

        # Impute the missing survival times using the Kaplan-Meier fit
        if(length(idx.censored) > 0) {
          Y.cal[idx.censored] <- sapply(idx.censored, function(i) sample_km_conditional(km_fit, Y.cal[i], 1))
        }

        # Update the imputed data frame with the new survival times
        data.cal.imputed$time <- Y.cal

        # Store the imputed dataset in the list
        imputed_data_list[[r]] <- data.cal.imputed
      }

      # Combine the results into one data frame with R many rows
      combined_data <- do.call(rbind, imputed_data_list)

      return(combined_data)
    }

    ## Call the function to obtain the de-censored data with R independent imputations
    data.cal.decensored <- impute_kaplan_meier(data.cal, km_fit, R)

    ## Apply CQR
    out <- predict_CQR(data.test, surv_model, data.cal.decensored, alpha)
    return(out)
}


sample_km_conditional <- function(km_fit, c, n) {
  # Extract the times and survival probabilities from the Kaplan-Meier fit
  km_times <- km_fit$time
  km_survival <- km_fit$surv

  # Condition on T > c (only keep times and probabilities where times > c)
  km_times_conditional <- km_times[km_times > c]
  km_survival_conditional <- km_survival[km_times > c]

  # Check if there are any times greater than c
  if(length(km_times_conditional) == 0) {
    stop("No survival times greater than c.")
  }

  # Normalize the conditional survival probabilities so they sum to 1
  km_survival_conditional <- km_survival_conditional / max(km_survival_conditional)

  # Inverse transform sampling function
  sample_from_km <- function(n, km_times, km_survival) {
    u <- runif(n)  # Uniform random numbers
    sapply(u, function(x) {
      # Find the time corresponding to the sampled probability
      idx <- max(which(km_survival >= x))
      return(km_times[idx])
    })
  }

  # Sample n random survival times from the conditioned Kaplan-Meier curve
  random_survival_times <- sample_from_km(n, km_times_conditional, km_survival_conditional)

  return(random_survival_times)
}
