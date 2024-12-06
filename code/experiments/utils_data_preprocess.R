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
    if (is.null(df)) return(NULL)

    ## Replace zero survival times with half of the smallest observed value
    min.time.nnz <- min(df$time[df$time > 0], na.rm = TRUE)
    df <- df %>%
        mutate(time = pmax(0.5 * min.time.nnz, time))

    ## Impute NAs with median for numeric and mode for non-numeric
    df <- df %>%
        mutate(across(
            everything(),
            ~ ifelse(is.na(.),
                     ifelse(is.numeric(.), median(., na.rm = TRUE), get_mode(.)),
                     .)
        ))

    ## Function to process factor/character columns
    process_factors <- function(df, threshold = 0.02) {
        df[] <- lapply(df, function(col) {
            if (is.character(col)) {
                col <- factor(col)  ## Convert characters to factors
            }
            if (is.factor(col)) {
                freq_table <- prop.table(table(col))  ## Calculate level frequencies
                rare_levels <- names(freq_table[freq_table < threshold])

                ## Handle binary factors
                if (length(levels(col)) == 2) {
                    if (length(rare_levels) == 1) {
                        ## Return NULL for binary factors with rare levels
                        return(NULL)
                    }
                }
                ## Merge rare levels into "Other" for factors with more than two levels
                if (length(rare_levels) > 0 && length(levels(col)) > 2) {
                    levels(col) <- c(levels(col), "Other")
                    col[col %in% rare_levels] <- "Other"
                    col <- factor(col)  ## Drop unused levels
                }
            }
            return(col)
        })

        ## Remove NULL columns (binary factors with rare levels removed)
        df <- Filter(Negate(is.null), df)

        ## Remove factors with fewer than 2 levels
        df <- df[, sapply(df, function(col) !(is.factor(col) && length(levels(col)) < 2))]

        return(as.data.frame(df))  ## Return as a data frame
    }

    ## Apply the function to clean the data frame
    threshold <- 0.02  ## Minimum frequency for rare levels
    df <- process_factors(df, threshold = threshold)

    ## Transform factors into dummy variables
    if (any(sapply(df, is.factor))) {  # Check if any factors exist before applying model.matrix
        df <- as_tibble(model.matrix(~ . - 1, data = df))
    }

    ## Remove redundant features using alias check
    cox_model <- coxph(Surv(time, status) ~ ., data = df)
    if (!is.null(alias(cox_model)$Complete)) {
        alias_cols <- which(alias(cox_model)$Complete)
        if (length(alias_cols) > 0) {
            df <- df[, -alias_cols]
        }
    }

    ## Function to remove highly correlated features
    remove_highly_correlated <- function(df, excluded_columns = c("time", "status"), cutoff = 0.75) {
        ## Identify columns to include in the correlation analysis
        cor_columns <- setdiff(names(df), excluded_columns)

        ## Iteratively remove highly correlated features
        while (TRUE) {
            ## Calculate the correlation matrix for the selected columns
            cor_matrix <- cor(df[, cor_columns, drop = FALSE], use = "pairwise.complete.obs")

            ## Find highly correlated features
            high_cor <- caret::findCorrelation(cor_matrix, cutoff = cutoff, verbose = FALSE)

            ## Break loop if no highly correlated features are found
            if (length(high_cor) == 0) break

            ## Remove the first correlated variable
            cor_columns <- cor_columns[-high_cor[1]]
        }

        ## Return the updated data frame, retaining excluded columns
        return(df[, c(excluded_columns, cor_columns)])
    }

    ## Apply the correlation filter
    df <- remove_highly_correlated(df, excluded_columns = c("time", "status"), cutoff = 0.75)

    ## Replace column names with consistent naming
    col.names <- c("time", "status", paste0("X", seq_len(ncol(df) - 2)))
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
    } else if (dataset.name=="COLON") {
        df <- survival::colon %>%
            select(-id, -study, -etype) %>%
            select(time, status, everything())
    } else if (dataset.name=="HEART") {
        df <- survival::heart %>%
            mutate(time=stop-start, status=event) %>%
            select(-id, -start, -stop, -event) %>%
            select(time, status, everything())
    } else if (dataset.name=="RETINOPATHY") {
        df <- survival::retinopathy %>%
            mutate(time=futime) %>%
            select(-id, -futime) %>%
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
