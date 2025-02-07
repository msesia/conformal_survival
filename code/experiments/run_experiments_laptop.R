suppressMessages(library(tidyverse))

source("experiment_1_laptop.R")

setup <- 1 # The only thing this does is determine the sub-directory where the results are saved


setting <- 7
surv_model_type <- "grf"
cens_model_type <- "cox"
num_samples_train <- 100
num_samples_cal <- 100
batch <- 1
batch_size <- 1
num_samples_test <- 10


# Run all the experiments
for(alpha in c(0.1,0.2)) {
  run_all_experiments(setting, surv_model_type, cens_model_type, num_samples_train, num_samples_cal, num_samples_test, batch, alpha, batch_size=batch_size)
}

# Load the results from sub-directory corresponding to given 'setup' number
load_data <- function(setup) {
  idir <- sprintf("results/setup_%d", setup)
  ifile.list <- list.files(idir)
  results <- do.call("rbind", lapply(ifile.list, function(ifile) {
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
  }))
  return(results)
}

results <- load_data(setup)
