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
#                   "DR-COSARC (fixed)",
                   "DR-COSARC"
                   )


results <- load_data()


make_table_small <- function(surv_model_plot, alpha_plot=0.1) {

    df <- results %>%
        filter(alpha==alpha_plot) %>%
        filter(Method %in% method.values) %>%
        mutate(Method = factor(Method, levels = method.values, labels = method.labels)) %>%
        mutate(Lower = `Coverage (lower bound)`, Upper=`Coverage (upper bound)`, LPB=`Mean lower bound`) %>%
        pivot_longer(c(Lower, Upper, LPB), names_to="metric", values_to="value") %>%
        group_by(dataset, surv_model_type, cens_model_type, train_prop_sub, alpha, Method, metric) %>%
        summarise(mean = mean(value, na.rm=T), se = 2*sd(value, na.rm=T)/sqrt(n())) %>%
        pivot_wider(names_from = metric, values_from = c(mean,se)) %>%
        mutate(mean_Point=0.5*(mean_Lower+mean_Upper), se_Point=sqrt(se_Lower^2+se_Upper^2)/2) %>%
        mutate(
            `Estimated Coverage` = ifelse(mean_Point+se_Point < 1-alpha_plot,
                              sprintf("\\textcolor{red}{%.2f (%.2f--%.2f})", mean_Point, mean_Lower, mean_Upper),
                              sprintf("%.2f (%.2f--%.2f)", mean_Point, mean_Lower, mean_Upper)),
            `Mean LPB` = sprintf("%.2f", mean_LPB)
        ) %>%
        select(-mean_Point, -mean_Lower, -mean_Upper, -se_Point, -se_LPB, -se_Lower, -se_Upper)

    df.tab <- df %>% filter(train_prop_sub==1, surv_model_type==surv_model_plot) %>%
        ungroup() %>%
        select(-surv_model_type, -cens_model_type, -train_prop_sub, -alpha) %>%
        select(dataset, Method, `Estimated Coverage`, `Mean LPB`) %>%
        arrange(dataset, Method)

    latex_table <- df.tab %>%
        select(-dataset) %>%
        kable("latex", booktabs = TRUE, longtable = FALSE, align = "c", escape = FALSE) %>%
        kable_styling(latex_options = c("hold_position", "repeat_header")) %>%
        group_rows(index = table(df.tab$dataset))

    output_file <- sprintf("tables/table_data_%s_alpha_%s_small.tex", surv_model_plot, alpha_plot)
    writeLines(latex_table, output_file)

}

make_table_large <- function(surv_model_plot, alpha_plot=0.1) {

    df <- results %>%
        filter(alpha==alpha_plot) %>%
        filter(Method %in% method.values) %>%
        mutate(Method = factor(Method, levels = method.values, labels = method.labels)) %>%
        mutate(LPB=`Mean lower bound`)  %>%
        pivot_longer(c(`Coverage (lower bound)`,`Coverage (upper bound)`, LPB), names_to="metric", values_to="value") %>%
        group_by(dataset, surv_model_type, cens_model_type, train_prop_sub, alpha, Method, metric) %>%
        summarise(mean = mean(value, na.rm=T), se = 2*sd(value, na.rm=T)/sqrt(n())) %>%
        pivot_wider(names_from = metric, values_from = c(mean,se)) %>%
        mutate(`Coverage (lower bound)` = sprintf("%.2f (%.2f)", `mean_Coverage (lower bound)`, `se_Coverage (lower bound)`),
               `Coverage (upper bound)` = sprintf("%.2f (%.2f)", `mean_Coverage (upper bound)`, `se_Coverage (upper bound)`),
               `LPB` = sprintf("%.2f (%.2f)", `mean_LPB`, `se_LPB`)) %>%
        select(`Coverage (lower bound)`, `Coverage (upper bound)`, LPB)

    df.tab <- df %>% filter(train_prop_sub==1, surv_model_type==surv_model_plot) %>%
        ungroup() %>%
        select(-surv_model_type, -cens_model_type, -train_prop_sub, -alpha) %>%
        select(dataset, Method, `Coverage (lower bound)`, `Coverage (upper bound)`, LPB) %>%
        arrange(dataset, Method)

    latex_table <- df.tab %>%
        select(-dataset) %>%
        kable("latex", booktabs = TRUE, longtable = FALSE, align = "c", escape = FALSE) %>%
        kable_styling(latex_options = c("hold_position", "repeat_header")) %>%
        group_rows(index = table(df.tab$dataset))

    output_file <- sprintf("tables/table_data_%s_alpha_%s_large.tex", surv_model_plot, alpha_plot)
    writeLines(latex_table, output_file)

}


for(surv_model_plot in c("grf", "cox", "survreg", "rf")) {
    make_table_small(surv_model_plot, alpha=0.1)
    make_table_large(surv_model_plot, alpha=0.1)
    make_table_small(surv_model_plot, alpha=0.2)
    make_table_large(surv_model_plot, alpha=0.2)
}


results <- load_data()

## Make small summary figure
pp <- results %>%
    filter(alpha==0.1) %>%
    filter(Method %in% method.values) %>%
    mutate(Method = factor(Method, levels = method.values, labels = method.labels)) %>%
    mutate(Lower = `Coverage (lower bound)`, Upper=`Coverage (upper bound)`, LPB=`Mean lower bound`) %>%
    pivot_longer(c(Lower, Upper, LPB), names_to="metric", values_to="value") %>%
    group_by(dataset, surv_model_type, cens_model_type, train_prop_sub, alpha, Method, metric) %>%
    summarise(mean = mean(value, na.rm=T), se = 2*sd(value, na.rm=T)/sqrt(n())) %>%
    pivot_wider(names_from = metric, values_from = c(mean,se)) %>%
    mutate(mean_Point=0.5*(mean_Lower+mean_Upper), se_Point=sqrt(se_Lower^2+se_Upper^2)/2) %>%
    ggplot(aes(x=Method, y=mean_Point)) +
    geom_boxplot() +
    geom_hline(aes(yintercept=1-alpha), color="red", linetype=2) +
    ylab("Estimated coverage") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))    
out.file <- sprintf("figures/data_summary_small.pdf")
ggsave(out.file, plot = pp, width = 4, height = 2.5, units = "in")


## Make large summary figure
pp <- results %>%
    filter(alpha %in% c(0.05, 0.1, 0.2)) %>%
    filter(Method %in% method.values) %>%
    mutate(Method = fct_rev(factor(Method, levels = method.values, labels = method.labels))) %>%
    mutate(Lower = `Coverage (lower bound)`, Upper=`Coverage (upper bound)`, LPB=`Mean lower bound`) %>%
    pivot_longer(c(Lower, Upper, LPB), names_to="metric", values_to="value") %>%
    group_by(dataset, surv_model_type, cens_model_type, train_prop_sub, alpha, Method, metric) %>%
    summarise(mean = mean(value, na.rm=T), se = 2*sd(value, na.rm=T)/sqrt(n())) %>%
    pivot_wider(names_from = metric, values_from = c(mean,se)) %>%
    mutate(mean_Point=0.5*(mean_Lower+mean_Upper), se_Point=sqrt(se_Lower^2+se_Upper^2)/2) %>%
    mutate(Target = sprintf("Target: %.2f", 1-alpha)) %>%
    ggplot(aes(y=Method, x=mean_Point)) +
    geom_boxplot() +
    geom_vline(aes(xintercept=1-alpha), color="red", linetype=2) +
    facet_grid(.~Target, scales="free")+
    xlab("Estimated coverage") +
    ylab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))    
out.file <- sprintf("figures/data_summary_full.pdf")
ggsave(out.file, plot = pp, width = 5, height = 2, units = "in")
