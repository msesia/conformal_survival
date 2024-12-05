options(width = 300)

library(tidyverse)
library(kableExtra)

load_data <- function() {
    idir <- sprintf("results_hpc/data")
    ifile.list <- list.files(idir)
    results <- do.call("rbind", lapply(ifile.list, function(ifile) {
        df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
    }))
    return(results)
}

method.values <- c("oracle", "nominal", "cqr", "cqr.decensor",
#                   "candes.oracle",
#                   "gui.oracle",
#                   "gui.oracle.cqr",
#                   "prototype.candes",
                   "prototype.gui"
#                   "prototype.gui.cqr",                   
                   )

method.labels <- c("Oracle", "Uncalibrated", "Naive CQR", "KM Decensoring",
#                   "Prototype (Candes, oracle)",
#                   "Prototype (Candes)",
#                   "Prototype (Gui et al., oracle)",
#                   "DR-COSAR (fixed)",
                   "DR-COSAR"                   
                   )


results <- load_data()

df <- results %>%    
    filter(Method %in% method.values) %>%
    mutate(Method = factor(Method, levels = method.values, labels = method.labels)) %>%
    pivot_longer(c("Coverage (lower bound)", "Coverage (upper bound)", "Mean lower bound"), names_to="metric", values_to="value") %>%
    group_by(dataset, surv_model_type, cens_model_type, train_prop_sub, alpha, Method, metric) %>%
    summarise(value.mean = mean(value), value.se = 2*sd(value)/sqrt(n())) %>%
    ##    mutate(value = sprintf("%.2f (%.2f)", value.mean, value.se)) %>%
    mutate(value = value.mean) %>%
    select(-value.mean, -value.se) %>%
    pivot_wider(names_from = metric, values_from = value) %>%
    mutate(Coverage.point = round(0.5*(as.numeric(`Coverage (lower bound)`) + as.numeric(`Coverage (upper bound)`)),2)) %>%    
    mutate(
        Coverage = ifelse(Coverage.point < 0.9,
                          sprintf("\\textcolor{red}{%.2f (%.2f--%.2f})", Coverage.point, `Coverage (lower bound)`, `Coverage (upper bound)`),
                          sprintf("%.2f (%.2f--%.2f)", Coverage.point, `Coverage (lower bound)`, `Coverage (upper bound)`)),
        `Mean LPB` = sprintf("%.2f", `Mean lower bound`)
    ) %>%
    select(-`Coverage (lower bound)`, -`Coverage (upper bound)`)

df.tab <- df %>% filter(train_prop_sub==1, surv_model_type=="rf") %>%
    ungroup() %>%
    select(-surv_model_type, -cens_model_type, -train_prop_sub, -alpha) %>%
    select(dataset, Method, Coverage, `Mean LPB`) %>%
    arrange(dataset, Method)

df.tab %>%
    select(-dataset) %>%
    kable("latex", booktabs = TRUE, longtable = FALSE, align = "c", escape = FALSE) %>%
    kable_styling(latex_options = c("hold_position", "repeat_header")) %>%
    group_rows(index = table(df.tab$dataset))    


