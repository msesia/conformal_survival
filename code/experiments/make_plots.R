options(width = 300)

library(tidyverse)

load_data <- function(setup) {
    idir <- sprintf("results_hpc/setup_%d", setup)
    ifile.list <- list.files(idir)
    results <- do.call("rbind", lapply(ifile.list, function(ifile) {
        df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
    }))
    return(results)
}

# Custom color palette using named colors
method_colors <- c(
  "Oracle" = "black",
  "None" = "grey",
  "Naive CQR" = "orchid1",
  "Qi et al." = "tan2",
  "Prototype (Candes, oracle censoring model)" = "blue",
  "Prototype (Candes)" = "blue4",
  "Prototype (Gui et al., oracle censoring model)" = "green",
  "Prototype" = "green4",
  "Prototype (Gui et al., CQR)" = "green3"
)

method.values <- c("oracle", "nominal", "cqr", "cqr.decensor",
#                   "candes.oracle",
##                   "prototype.candes",
#                   "gui.oracle",
#                   "gui.oracle.cqr",
                   "prototype.gui"
#                   "prototype.gui.cqr",                   
                   )

method.labels <- c("Oracle", "None", "Naive CQR", "Qi et al.",
#                   "Prototype (Candes, oracle)",
#                   "Prototype (Candes)",
#                   "Prototype (Gui et al., oracle)",
                   "Prototype"                   
                   )

setting.values <- c(8,7,10,9,1:6)
setting.labels <- c(1:10)

make_figure_0 <- function(plot_surv_model=NULL, plot_cens_model="grf", plot_n_train=500, plot_n_cal=500) {
    setting.list <- c(9,1:6)
    df <- results %>% filter(setting %in% setting.list, surv_model_type==plot_surv_model)
    out.file <- sprintf("figures/exp1_fig0_bp_surv_%s_cens_%s_n%d_n%d.pdf",
                        plot_surv_model, plot_cens_model, plot_n_train, plot_n_cal)
    plot.height <- 8
    summary <- df %>%
        mutate(`Surv. model` = surv_model_type) %>%
        filter(setting %in% setting.list, cens_model_type==plot_cens_model, n_train==plot_n_train, n_cal==plot_n_cal) %>%
        group_by(setup, setting, n_train, n_train_cens, n_cal, num_feat_censor, alpha) %>%
        mutate(`Mean lower bound` = `Mean lower bound` / mean(`Mean lower bound`[Method=="oracle"])) %>%
        ungroup() %>%
        pivot_longer(c("Coverage (event time)", "Mean lower bound"), names_to="metric", values_to="value") %>%
        mutate(metric = factor(metric, c("Coverage (event time)", "Mean lower bound"), c("Coverage", "Lower Bound"))) %>%
        filter(Method %in% method.values) %>%
        mutate(Method = fct_rev(factor(Method, levels = method.values, labels = method.labels))) %>%
        mutate(Setting = factor(setting, setting.values, setting.labels))
    if(nrow(summary)==0) return(0)
    df.limits <- tibble(metric=c("Coverage","Coverage", "Lower Bound", "Lower Bound"), value=c(0,1,0,1), Method="Oracle")
    ## Plotting the data
    if(nrow(summary)==0) return()
    pp <- summary %>%
        ggplot(aes(y = Method, x = value, color=Method)) +
        geom_boxplot() +
        scale_color_manual(values = method_colors) + # Apply custom color scale mapped to Method labels
        geom_vline(data = summary %>% filter(metric=="Coverage"), aes(xintercept = 0.9), linetype = "dashed", color="red") +
        geom_point(data=df.limits, alpha=0) +
        xlab("") +
        scale_x_continuous(breaks = seq(0, 1, by = 0.25)) +
        theme_bw() +
        theme(legend.position = "none") +
        facet_grid(Setting ~ metric, scales = "free", labeller = labeller(Setting = label_both, metric = label_value))
    ggsave(out.file, plot = pp, width = 6, height = plot.height, units = "in")
}
results <- load_data(0)
make_figure_0(plot_cens_model="grf", plot_surv_model="grf", plot_n_train=1000, plot_n_cal=1000)


make_figure_1 <- function(plot_surv_model=NULL, plot_cens_model="grf", plot_n_train=500, plot_n_cal=500) {
    setting.list <- c(7,8,10)
    df <- results %>% filter(setting %in% setting.list, surv_model_type==plot_surv_model)
    out.file <- sprintf("figures/exp1_fig1_bp_surv_%s_cens_%s_n%d_n%d.pdf",
                        plot_surv_model, plot_cens_model, plot_n_train, plot_n_cal)
    plot.height <- 3.5
    summary <- df %>%
        mutate(`Surv. model` = surv_model_type) %>%
        filter(setting %in% setting.list, cens_model_type==plot_cens_model, n_train==plot_n_train, n_cal==plot_n_cal) %>%
        group_by(setup, setting, n_train, n_train_cens, n_cal, num_feat_censor, alpha) %>%
        mutate(`Mean lower bound` = `Mean lower bound` / mean(`Mean lower bound`[Method=="oracle"])) %>%
        ungroup() %>%
        pivot_longer(c("Coverage (event time)", "Mean lower bound"), names_to="metric", values_to="value") %>%
        mutate(metric = factor(metric, c("Coverage (event time)", "Mean lower bound"), c("Coverage", "Lower Bound"))) %>%
        filter(Method %in% method.values) %>%
        mutate(Method = fct_rev(factor(Method, levels = method.values, labels = method.labels))) %>%
        mutate(Setting = factor(setting, setting.values, setting.labels))
    if(nrow(summary)==0) return(0)
    df.limits <- tibble(metric=c("Coverage","Coverage", "Lower Bound", "Lower Bound"), value=c(0,1,0,1), Method="Oracle")
    ## Plotting the data
    if(nrow(summary)==0) return()
    pp <- summary %>%
        ggplot(aes(y = Method, x = value, color=Method)) +
        geom_boxplot() +
        scale_color_manual(values = method_colors) + # Apply custom color scale mapped to Method labels
        geom_vline(data = summary %>% filter(metric=="Coverage"), aes(xintercept = 0.9), linetype = "dashed", color="red") +
        geom_point(data=df.limits, alpha=0) +
        xlab("") +
        scale_x_continuous(breaks = seq(0, 1, by = 0.25)) +
        theme_bw() +
        theme(legend.position = "none") +
        facet_grid(Setting ~ metric, scales = "free", labeller = labeller(Setting = label_both, metric = label_value))
    ggsave(out.file, plot = pp, width = 6, height = plot.height, units = "in")
}
results <- load_data(1)
make_figure_1(plot_cens_model="grf", plot_surv_model="grf", plot_n_train=1000, plot_n_cal=1000)



make_boxplot <- function(plot_setting=1, plot_surv_model=NULL, plot_cens_model="grf", plot_n_train=500, plot_n_cal=500) {
    if(!is.null(plot_surv_model)) {
        df <- results %>% filter(setting==plot_setting, surv_model_type==plot_surv_model)
        out.file <- sprintf("figures/exp1_%d_bp_surv_%s_cens_%s_n%d_n%d.pdf",
                            plot_setting, plot_surv_model, plot_cens_model, plot_n_train, plot_n_cal)
        plot.height <- 2
    } else {
        df <- results
        out.file <- sprintf("figures/exp1_%d_bp_surv_all_cens_%s_n%d_n%d.pdf", plot_setting, plot_cens_model, plot_n_train, plot_n_cal)
        plot.height <- 6.5
    }
    summary <- df %>%
        mutate(`Surv. model` = surv_model_type) %>%
        filter(setting==plot_setting, cens_model_type==plot_cens_model, n_train==plot_n_train, n_cal==plot_n_cal) %>%
        pivot_longer(c("Coverage (event time)", "Mean lower bound"), names_to="metric", values_to="value") %>%
        mutate(metric = factor(metric, c("Coverage (event time)", "Mean lower bound"), c("Coverage", "Lower Bound"))) %>%
        filter(Method %in% method.values) %>%
        mutate(Method = fct_rev(factor(Method, levels = method.values, labels = method.labels)))
    if(nrow(summary)==0) return(0)
    df.limits <- tibble(metric=c("Coverage","Coverage", "Lower Bound"), value=c(0,1,0), Method="Oracle")
    ## Plotting the data
    if(nrow(summary)==0) return()
    pp <- summary %>%
        ggplot(aes(y = Method, x = value, color=Method)) +
        geom_boxplot() +
        scale_color_manual(values = method_colors) + # Apply custom color scale mapped to Method labels
        geom_vline(data = summary %>% filter(metric=="Coverage"), aes(xintercept = 0.9), linetype = "dashed", color="red") +
        geom_point(data=df.limits, alpha=0) +
        xlab("") +
        theme_bw() +
        theme(legend.position = "none")
    if(!is.null(plot_surv_model)) {
        pp <- pp +
        facet_grid(. ~ metric, scales = "free")
    } else {
        pp <- pp +
        facet_grid(`Surv. model` ~ metric, scales = "free",
                   labeller = labeller(`Surv. model` = label_both, metric = label_value))
    }
    ggsave(out.file, plot = pp, width = 6, height = plot.height, units = "in")
}

results <- load_data(1)

for(pps in 1:10) {
    make_boxplot(plot_setting=pps, plot_cens_model="grf", plot_n_train=1000, plot_n_cal=1000)
    make_boxplot(plot_setting=pps, plot_surv_model="grf", plot_cens_model="grf", plot_n_train=1000, plot_n_cal=1000)
    make_boxplot(plot_setting=pps, plot_cens_model="cox", plot_n_train=1000, plot_n_cal=1000)
}



make_plot_2 <- function(plot_setting=1) {
    df <- results %>% filter(setting %in% plot_setting)
    summary <- df %>%
        pivot_longer(c("Coverage (event time)", "Mean lower bound"), names_to="metric", values_to="value") %>%
        mutate(metric = factor(metric, c("Coverage (event time)", "Mean lower bound"), c("Coverage", "Lower Bound"))) %>%
        filter(Method %in% method.values) %>%
        mutate(Method = factor(Method, levels = method.values, labels = method.labels)) %>%
        group_by(setup, setting, n_train, n_train_cens, n_cal, num_feat_censor, alpha, Method, metric) %>%
        summarise(SE = sd(value)/sqrt(n()), value=mean(value), N=n()) %>%
        mutate(Setting = factor(setting, setting.values, setting.labels))        
    df.limits <- tibble(metric=c("Coverage","Coverage", "Lower Bound"), value=c(0.3,1,0), n_train_cens=100, Method="Oracle")
    ## Plotting the data
    if(nrow(summary)==0) return()
    pp <- summary %>%
        filter(n_train_cens >= 50) %>%
        ggplot(aes(y = value, x = n_train_cens, color=Method, shape=Method)) +
        geom_point() +
        geom_line() +
#        geom_errorbar(aes(ymin = value - 2*SE, ymax = value + 2*SE), width = 0.2) +
        scale_color_manual(values = method_colors) + # Apply custom color scale mapped to Method labels
##        geom_hline(data = summary %>% filter(metric=="Coverage"), aes(yintercept = 0.9), linetype = "dashed", color="red") +
        geom_point(data=df.limits, alpha=0) +
        scale_x_log10() +
        xlab("Number of training samples for censoring model") +
        ylab("") +
        theme_bw()
    out.file <- sprintf("figures/exp2_%s.pdf", paste(plot_setting,collapse="_"))
    if(length(plot_setting)==1) {
        plot.height <- 2
        pp <- pp + facet_wrap(metric ~ ., scales = "free", nrow=1)

    } else {
        plot.height <- 4
        pp <- pp + facet_wrap(Setting ~ metric, scales = "free", labeller = labeller(Setting = label_both, metric = label_value))
    }
    ggsave(out.file, plot = pp, width = 6, height = plot.height, units = "in")
}

results <- load_data(2)

make_plot_2(plot_setting=8)
make_plot_2(plot_setting=c(7,10))


make_plot_3 <- function(plot_setting=1) {
    df <- results %>% filter(setting %in% plot_setting)
    summary <- df %>%
        pivot_longer(c("Coverage (event time)", "Mean lower bound"), names_to="metric", values_to="value") %>%
        mutate(metric = factor(metric, c("Coverage (event time)", "Mean lower bound"), c("Coverage", "Lower Bound"))) %>%
        filter(Method %in% method.values, surv_model_type=="grf") %>%
        mutate(Method = factor(Method, levels = method.values, labels = method.labels)) %>%
        group_by(setup, setting, n_train, n_train_cens, n_cal, num_feat_censor, alpha, Method, metric, surv_model_type) %>%
        summarise(SE = sd(value)/sqrt(n()), value=mean(value), N=n()) %>%
        filter(n_train>=50) %>%
        mutate(Setting = factor(setting, setting.values, setting.labels))                
    df.limits <- tibble(metric=c("Coverage","Coverage", "Lower Bound"), value=c(0.3,1,0), n_train=100, Method="Oracle")
    ## Plotting the data
    if(nrow(summary)==0) return()
    pp <- summary %>%
        ggplot(aes(y = value, x = n_train, color=Method, shape=Method)) +
        geom_point() +
        geom_line() +
#        geom_errorbar(aes(ymin = value - 2*SE, ymax = value + 2*SE), width = 0.2) +
        scale_color_manual(values = method_colors) + # Apply custom color scale mapped to Method labels
##        geom_hline(data = summary %>% filter(metric=="Coverage"), aes(yintercept = 0.9), linetype = "dashed", color="red") +
        geom_point(data=df.limits, alpha=0) +
        scale_x_log10() +
        xlab("Number of training samples") +
        theme_bw()
    out.file <- sprintf("figures/exp3_%s.pdf", paste(plot_setting,collapse="_"))
    if(length(plot_setting)==1) {
        plot.height <- 2
        pp <- pp + facet_wrap(metric ~ ., scales = "free", nrow=1)

    } else {
        plot.height <- 4
        pp <- pp + facet_wrap(Setting ~ metric, scales = "free", labeller = labeller(Setting = label_both, metric = label_value))
    }
    ggsave(out.file, plot = pp, width = 6, height = plot.height, units = "in")
}

results <- load_data(3)

make_plot_3(plot_setting=8)
make_plot_3(plot_setting=c(7,10))




make_plot_4 <- function(plot_setting=1) {
    df <- results %>% filter(setting %in% plot_setting)
    summary <- df %>%
        pivot_longer(c("Coverage (event time)", "Mean lower bound"), names_to="metric", values_to="value") %>%
        mutate(metric = factor(metric, c("Coverage (event time)", "Mean lower bound"), c("Coverage", "Lower Bound"))) %>%
        filter(Method %in% method.values) %>%
        mutate(Method = factor(Method, levels = method.values, labels = method.labels)) %>%
        group_by(setup, setting, n_train, n_train_cens, n_cal, num_feat_censor, alpha, Method, metric) %>%
        summarise(SE = sd(value)/sqrt(n()), value=mean(value), N=n()) %>%
        mutate(Setting = factor(setting, setting.values, setting.labels))                        
    df.limits <- tibble(metric=c("Coverage","Coverage","Lower Bound"), value=c(0.3,1,0), n_cal=100, Method="Oracle")
    ## Plotting the data
    if(nrow(summary)==0) return()
    pp <- summary %>%
        filter(n_train>=100) %>%
        ggplot(aes(y = value, x = n_cal, color=Method, shape=Method)) +
        geom_point(alpha=0.75) +
        geom_line(alpha=0.75) +
##        geom_errorbar(aes(ymin = value - 2*SE, ymax = value + 2*SE), width = 0.2) +
        scale_color_manual(values = method_colors) + # Apply custom color scale mapped to Method labels
##        geom_hline(data = summary %>% filter(metric=="Coverage"), aes(yintercept = 0.9), linetype = "dashed", color="red") +
        geom_point(data=df.limits, alpha=0) +
        facet_wrap(metric ~ ., scales = "free", nrow=1) +
        scale_x_log10() +
        xlab("Number of calibration samples") +
        ylab("") +
        theme_bw()
    out.file <- sprintf("figures/exp4_%s.pdf", paste(plot_setting,collapse="_"))
    if(length(plot_setting)==1) {
        plot.height <- 2
        pp <- pp + facet_wrap(metric ~ ., scales = "free", nrow=1)

    } else {
        plot.height <- 4
        pp <- pp + facet_wrap(Setting ~ metric, scales = "free", labeller = labeller(Setting = label_both, metric = label_value))
    }
    ggsave(out.file, plot = pp, width = 6, height = plot.height, units = "in")
}

results <- load_data(4)

make_plot_4(plot_setting=8)
make_plot_4(plot_setting=c(7,10))



make_plot_5 <- function(plot_setting=1) {
    df <- results %>% filter(setting %in% plot_setting)
    summary <- df %>%
        pivot_longer(c("Coverage (event time)", "Mean lower bound"), names_to="metric", values_to="value") %>%
        mutate(metric = factor(metric, c("Coverage (event time)", "Mean lower bound"), c("Coverage", "Lower Bound"))) %>%
        filter(Method %in% method.values) %>%
        mutate(Method = factor(Method, levels = method.values, labels = method.labels)) %>%
        group_by(setup, setting, n_train, n_train_cens, n_cal, num_feat_censor, alpha, Method, metric) %>%
        summarise(SE = sd(value)/sqrt(n()), value=mean(value), N=n()) %>%
        mutate(Setting = factor(setting, setting.values, setting.labels))                                
    df.limits <- tibble(metric=c("Coverage","Coverage", "Lower Bound"), value=c(0.3,1,0), num_feat_censor=100, Method="Oracle")
    ## Plotting the data
    if(nrow(summary)==0) return()
    pp <- summary %>%
        ggplot(aes(y = value, x = num_feat_censor, color=Method, shape=Method)) +
        geom_point(alpha=0.75) +
        geom_line(alpha=0.75) +
##        geom_errorbar(aes(ymin = value - 2*SE, ymax = value + 2*SE), width = 0.2) +
        scale_color_manual(values = method_colors) + # Apply custom color scale mapped to Method labels
##        geom_hline(data = summary %>% filter(metric=="Coverage"), aes(yintercept = 0.9), linetype = "dashed", color="red") +
        geom_point(data=df.limits, alpha=0) +
        facet_wrap(metric ~ ., scales = "free", nrow=1) +
        scale_x_log10() +
        xlab("Number of covarariates in the censoring model") +
        ylab("") +
        theme_bw()
    out.file <- sprintf("figures/exp5_%s.pdf", paste(plot_setting,collapse="_"))
    if(length(plot_setting)==1) {
        plot.height <- 2
        pp <- pp + facet_wrap(metric ~ ., scales = "free", nrow=1)

    } else {
        plot.height <- 4
        pp <- pp + facet_wrap(Setting ~ metric, scales = "free", labeller = labeller(Setting = label_both, metric = label_value))
    }
    ggsave(out.file, plot = pp, width = 6, height = plot.height, units = "in")
}

results <- load_data(5)

make_plot_5(plot_setting=8)
make_plot_5(plot_setting=c(7,10))

for(pps in 1:10) {
    make_plot_5(plot_setting=pps)
}




# Custom shape mapping (optional)
method_shapes <- c(
  "Oracle" = 16,  # Circle
  "None" = 17, # Triangle
  "CQR" = 15,     # Square
  "Qi et al." = 15, # Square (same as CQR)
  "Candes (oracle)" = 16,  # Circle
  "Candes (oracle, tuned)" = 16,  # Circle
  "Candes (prototype)" = 18, # Diamond
  "Candes (prototype, tuned)" = 18, # Diamond
  "Gui et al. (oracle)" = 19, # Filled circle
  "Gui et al. (oracle, cqr)" = 15, # Square
  "Gui et al. (prototype)" = 18, # Diamond (same as Candes prototype)
  "Gui et al. (prototype, cqr)" = 15 # Square
)


make_plot_vs_n <- function(plot_n_cal=500) {
    summary <- results %>%
        pivot_longer(c("Coverage (observed time)", "Mean lower bound", "Median lower bound", "Mean lower bound (cover)", "Median lower bound (cover)",
                       "Coverage (event time)", "Oracle MSE"),
                     names_to="metric", values_to="value") %>%
        group_by(setup, setting, n_train, n_cal, alpha, Method, metric) %>%
        summarise(SE = sd(value)/sqrt(n()), value=mean(value))
    ## Plotting the data
    summary %>%
        filter(setting==1, n_cal==plot_n_cal, metric != "Coverage (observed time)") %>%
        filter(Method %in% c("nominal", "cqr", "cqr.decensor", "candes.oracle", "gui.oracle")) %>%
        mutate(Method = factor(Method, levels = method.values, labels = method.labels)) %>%
        ggplot(aes(x = n_train, y = value, color = Method, shape = Method)) +
        geom_point(alpha=0.5) +
        geom_line(alpha=0.5) +
        geom_errorbar(aes(ymin = value - 2*SE, ymax = value + 2*SE), width = 0.2) +
        facet_wrap(metric ~ ., scales = "free", nrow=1) +
        scale_color_manual(values = method_colors) + # Apply custom color scale mapped to Method labels
        scale_shape_manual(values = method_shapes) + # Apply custom shape scale (optional)
        geom_hline(data = summary %>% filter(metric=="Coverage (event time)"), aes(yintercept = 0.9), linetype = "dashed") +
        scale_x_log10() +
        theme_bw()

}

make_plot_vs_n(plot_n_cal=200)


## init_settings <- function(idx.exclude=NULL, names_ACODE=FALSE) {
##     cbPalette <<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#6e57d2", "red")
##     method.values <<- c("lb_simes", "lb_storey_simes", "lb_fisher", "lb_wmw_k2", "lb_wmw_k3", "lb_lmp", "lb_auto") # "lb_wmw_k3", "lb_wmw_k4",
##     if(names_ACODE) {
##         method.labels <<- c("ACODE (Simes)", "ACODE (Storey-Simes)", "ACODE (Fisher)", "ACODE (WMW)", "ACODE (LMP, Lehmann k=3)", "ACODE (LMP,G-hat)", "ACODE (Adaptive)")
##     } else {
##         method.labels <<- c("Simes", "Storey-Simes", "Fisher", "WMW", "LMPI (Lehmann, k=3)", "LMPI (G-hat)", "Adaptive") #"WMW (k=3)",
##     }
##     classifier.values <<- c("occ-auto", "bc-auto", "auto")
##     classifier.labels <<- c("One-Class", "Binary", "Automatic")
##     color.scale <<- cbPalette[c(1,1,4,3,6,7,8)]
##     shape.scale <<- c(2,6,3,1,0,9,8)
##     alpha.scale <<- c(0.75,0.75,0.75,0.75,0.75,0.75,1)
##     data.values <<- c("pendigits", "creditcard", "cover", "shuttle", "mammography", "aloi")
##     data.labels <<- c("Pendigits", "Creditcard", "Covertype", "Shuttle", "Mammography", "ALOI")
##     if(length(idx.exclude)>0) {
##         ## Exclude Storey-simes
##         method.values <<- method.values[-idx.exclude]
##         method.labels <<- method.labels[-idx.exclude]
##         color.scale <<- color.scale[-idx.exclude]
##         shape.scale <<- shape.scale[-idx.exclude]
##         alpha.scale <<- alpha.scale[-idx.exclude]
##     }
## }
