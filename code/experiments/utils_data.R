suppressMessages(library(tidyverse))

load_csv <- function(data.name) {
    data <- read_csv(sprintf("../../data/data_csv/%s.csv", data.name), show_col_types=FALSE)
    return(data)
}

split_data <- function(data, train_prop = 0.6, cal_prop = 0.2, test_prop = 0.2, seed = NULL) {
  # Ensure proportions sum to 1
  if (abs(train_prop + cal_prop + test_prop - 1) > .Machine$double.eps) {
    stop("Proportions must sum to 1")
  }

  # Set seed for reproducibility, if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Determine the total number of rows
  n <- nrow(data)

  # Calculate sizes for each split
  train_size <- round(train_prop * n)
  cal_size <- round(cal_prop * n)

  # Shuffle data
  shuffled_data <- data %>% slice_sample(n = n)

  # Split the data
  train_data <- shuffled_data[1:train_size, ]
  cal_data <- shuffled_data[(train_size + 1):(train_size + cal_size), ]
  test_data <- shuffled_data[(train_size + cal_size + 1):n, ]

  # Return a list of data subsets
  list(train = train_data, cal = cal_data, test = test_data)
}
