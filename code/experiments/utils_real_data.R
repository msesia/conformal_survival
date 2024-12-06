suppressMessages(library(tidyverse))
suppressMessages(library(survival))
suppressMessages(library(caret))
suppressMessages(library(car))
suppressMessages(library(KMsurv))

##library(hdf5r)

# Helper function to calculate mode
get_mode <- function(x) {
    uniq_vals <- na.omit(x) # Remove NA values
    uniq_vals[which.max(tabulate(match(x, uniq_vals)))]
}


preprocess_data <- function(df) {
    if(is.null(df)) return(NULL)

    ## Replace zero survival times with half of smallest observed value
    min.time.nnz <- min(df$time[df$time>0])
    df <- df %>%
        mutate(time = pmax(0.5*min.time.nnz, time))

    ## Impute NAs with median for numeric and mode for non-numeric
    df <- df %>%
        mutate(across(
            everything(),
            ~ if_else(is.na(.),
                      if (is.numeric(.)) median(., na.rm = TRUE) else get_mode(.),
                      .),
            .names = "{col}" # Optional: create imputed columns
        ))
    ## Remove columns with only one unique value (for factor or character columns)
    df <- df[, !(sapply(df, function(col) {
        is.character(col) || is.factor(col)
    }) & sapply(df, function(col) {
        length(unique(col)) == 1
    }))]

    ## Transform factors in to dummies
    df <- as_tibble(model.matrix(~ . - 1, data = df))
    ## Remove redundant features
    cox_model <- coxph(Surv(time, status) ~ ., data = df)
    if(!is.null(alias(cox_model)$Complete)) {
        alias_cols <- which(alias(cox_model)$Complete)
        if(length(alias_cols)>0) {
            df <- df[,-alias_cols]
        }
    }

    remove_rare_binary <- function(df, threshold = 0.02) {
        df[, !sapply(df, function(col) {
            ## Ensure the column is binary
            if (length(unique(col)) == 2) {
                ## Calculate the frequency of each value
                freq_table <- table(col)
                min_freq <- min(freq_table) / sum(freq_table)
                ## Return TRUE if the minor value frequency is below the threshold
                return(min_freq < threshold)
            }
            return(FALSE)  # Not binary
        })]
    }
    df <- remove_rare_binary(df, threshold = 0.02)

    ## Remove highly correlated features
    cor_matrix <- cor(df, use = "pairwise.complete.obs")
    high_cor <- caret::findCorrelation(cor_matrix, cutoff = 0.75)  # Correlation threshold
    if(length(high_cor)>0) {
        df <- df[, -high_cor]
    }

    ## Replace col names
    col.names <- c("time", "status", paste("X", 1:(ncol(df)-2), sep = ""))
    colnames(df) <- col.names

    return(df)
}

load_data <- function(dataset.name) {
    if (dataset.name == "VALCT") {
        ## Load data from survival::veteran
        ## 137 observations, 6 variables
        df <- survival::veteran |> select(time, status, everything())
    } else if (dataset.name == "PBC") {
        ## Load data from survival::PBC
        ## 137 observations, 6 variables
        df <- survival::pbc |> select(time, status, everything()) %>% select(-id) %>%
                                   mutate(status = ifelse(status==0, 0, 1))
    } else if (dataset.name=="OVARIAN") {
        df <- survival::ovarian %>%
            mutate(time=futime, status=fustat) %>%
            select(time, status, age, resid.ds, rx, ecog.ps)
    } else if (dataset.name=="BMT") {
        data(bmt)
        df <- bmt %>% mutate(time=t1, status=d1) %>%
            select(time, status, group, starts_with("z"))
    } else if (dataset.name=="GBSG2") {
        df <- TH.data::GBSG2 %>%
            mutate(status=cens) %>%
            select(-cens) %>%
            select(time, status, everything())        
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
    } else if (dataset.name=="METABRIC") {
        ifile <- "../../data/brca_metabric_clinical_data.tsv"
        df <- read_tsv(ifile)
        df <- df %>% mutate(time=`Overall Survival (Months)`, status=`Overall Survival Status`) %>%
            mutate(status = ifelse(status=="0:LIVING",0,1)) %>%
            select(time, status, everything()) %>%
            filter(!is.na(time), !is.na(status)) %>%
            select(-c("Study ID", "Patient ID", "Sample ID", "Patient's Vital Status", 'Overall Survival (Months)', 'Overall Survival Status', 'Patient\'s Vital Status',
                      'Relapse Free Status (Months)', 'Relapse Free Status', 'Cohort', 'Number of Samples Per Patient', 'Sample Type', 'Oncotree Code', 'Cancer Type'))
    } else {
        df <- NULL
        return(NULL)
    }

    return(df)
}

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
