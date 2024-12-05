library(tidyverse)
#library(hdf5r)

# Helper function to calculate mode
get_mode <- function(x) {
    uniq_vals <- na.omit(x) # Remove NA values
    uniq_vals[which.max(tabulate(match(x, uniq_vals)))]
}

load_data <- function(dataset.name) {
    if (dataset.name == "VALCT") {
        # Load data from survival::veteran
        # 137 observations, 6 variables
        df <- survival::veteran |> select(time, status, everything())
        col.names <- c("time", "status", paste("X", 1:(ncol(df)-2), sep = ""))
        colnames(df) <- col.names
    } else if (dataset.name == "PBC") {
        # Load data from survival::PBC
        # 137 observations, 6 variables
        df <- survival::pbc |> select(time, status, everything()) %>% select(-id) %>%
            mutate(status = ifelse(status==0, 0, 1))
        # Impute NAs with median for numeric and mode for non-numeric
        df <- df %>%
          mutate(across(
            everything(),
            ~ if_else(is.na(.),
                      if (is.numeric(.)) median(., na.rm = TRUE) else get_mode(.),
                      .),
            .names = "{col}" # Optional: create imputed columns
          ))
        col.names <- c("time", "status", paste("X", 1:(ncol(df)-2), sep = ""))
        colnames(df) <- col.names
    } else if (dataset.name == "GBSG") {
        ifile <- "../../data/gbsg_cancer_train_test.h5"
        h5_file <- hdf5r::H5File$new(ifile, mode = "r")
        train_group <- h5_file[["train"]]
        test_group <- h5_file[["test"]]
        dataset.1.x <- train_group[["x"]]
        dims.1.x <- dataset.1.x$dims
        data.1.x <- dataset.1.x[rep(TRUE, dims.1.x[1]), rep(TRUE, dims.1.x[2])]
        dataset.2.x <- test_group[["x"]]
        dims.2.x <- dataset.2.x$dims
        data.2.x <- dataset.2.x[rep(TRUE, dims.2.x[1]), rep(TRUE, dims.2.x[2])]
        data.x <- rbind(t(data.1.x), t(data.2.x))
        data.e <- c(train_group[["e"]][], test_group[["e"]][])
        data.t <- c(train_group[["t"]][], test_group[["t"]][])
        df <- tibble(time=data.t, status=data.e) %>% cbind(data.x)
        colnames(df) <- c("time", "status", paste("X", 1:7, sep = ""))
    } else if (dataset.name=="METABRIC") {
        ifile <- "../../data/brca_metabric_clinical_data.tsv"
        df <- read_tsv(ifile)
        df <- df %>% mutate(time=`Overall Survival (Months)`, status=`Overall Survival Status`) %>%
            mutate(status = ifelse(status=="0:LIVING",0,1)) %>%
            select(time, status, everything()) %>%
            filter(!is.na(time), !is.na(status)) %>%
            select(-c("Study ID", "Patient ID", "Sample ID", "Patient's Vital Status", 'Overall Survival (Months)', 'Overall Survival Status', 'Patient\'s Vital Status',
                     'Relapse Free Status (Months)', 'Relapse Free Status', 'Cohort', 'Number of Samples Per Patient', 'Sample Type', 'Oncotree Code', 'Cancer Type'))
        col.names <- c("time", "status", paste("X", 1:(ncol(df)-2), sep = ""))
        colnames(df) <- col.names
    } else {
        df <- NULL
    }
    return(df)
}

load_csv <- function(data.name) {
    if(data.name=="VALCT") {
        data <- read_csv(sprintf("../../data/data_csv/%s.csv", data.name),
                         col_types = cols(
                             time = col_double(),
                             status = col_double(),
                             trt = col_double(),
                             celltype = col_factor(),
                             karno = col_double(),
                             diagtime = col_double(),
                             age = col_double(),
                             prior = col_double()
                         )
                         )
    } else if(data.name=="PBC") {
        data <- read_csv(sprintf("../../data/data_csv/%s.csv", data.name),
                         col_types = cols(
                             time = col_double(),
                             status = col_double(),
                             trt = col_double(),
                             age = col_double(),
                             sex = col_factor(),
                             ascites = col_double(),
                             hepato = col_double(),
                             spiders = col_double(),
                             edema = col_double(),
                             bili = col_double(),
                             chol = col_double(),
                             albumin = col_double(),
                             copper = col_double(),
                             alk.phos = col_double(),
                             ast = col_double(),
                             trig = col_double(),
                             platelet = col_double(),
                             protime = col_double(),
                             stage = col_double()
                         )
                         )
    } else if(data.name=="GBSG") {
        data <- read_csv(sprintf("../../data/data_csv/%s.csv", data.name),
                         col_types = cols(
                             time = col_double(),
                             status = col_double(),
                             `1` = col_double(),
                             `2` = col_double(),
                             `3` = col_double(),
                             `4` = col_double(),
                             `5` = col_double(),
                             `6` = col_double(),
                             `7` = col_double()
                         )
                         )
    } else if(data.name=="METABRIC") {
        data <- read_csv(sprintf("../../data/data_csv/%s.csv", data.name),
                         col_types = cols(
                             time = col_double(),
                             status = col_double(),
                             `Age at Diagnosis` = col_double(),
                             `Type of Breast Surgery` = col_factor(),
                             `Cancer Type Detailed` = col_factor(),
                             Cellularity = col_factor(),
                             Chemotherapy = col_factor(),
                             `Pam50 + Claudin-low subtype` = col_factor(),
                             `ER status measured by IHC` = col_factor(),
                             `ER Status` = col_factor(),
                             `Neoplasm Histologic Grade` = col_double(),
                             `HER2 status measured by SNP6` = col_factor(),
                             `HER2 Status` = col_factor(),
                             `Tumor Other Histologic Subtype` = col_factor(),
                             `Hormone Therapy` = col_factor(),
                             `Inferred Menopausal State` = col_factor(),
                             `Integrative Cluster` = col_factor(),
                             `Primary Tumor Laterality` = col_factor(),
                             `Lymph nodes examined positive` = col_double(),
                             `Mutation Count` = col_double(),
                             `Nottingham prognostic index` = col_double(),
                             `PR Status` = col_factor(),
                             `Radio Therapy` = col_factor(),
                             `3-Gene classifier subtype` = col_factor(),
                             `TMB (nonsynonymous)` = col_double(),
                             `Tumor Size` = col_double(),
                             `Tumor Stage` = col_double()
                         )
                         ) %>%
            mutate(time = pmax(0.001, time))            
    } else {
        data <- NULL
    }

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
