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
  "Nominal" = "grey",
  "Naive CQR" = "orchid1",
  "Qi et al." = "tan2",
  "Prototype (Candes, oracle censoring model)" = "blue",
  "Prototype (Candes)" = "blue4",
  "Prototype (Gui et al., oracle censoring model)" = "green",
  "Prototype (Gui et al.)" = "green4"
)

method.values <- c("nominal", "cqr.decensor", "cqr", 
#                   "candes.oracle",
                   "prototype.candes",
#                   "gui.oracle",
                   "prototype.gui",
                   "oracle"
                   )

method.labels <- c("Nominal", "Qi et al.", "Naive CQR", 
#                   "Prototype (Candes, oracle)",
                   "Prototype (Candes)",
#                   "Prototype (Gui et al., oracle)",
                   "Prototype (Gui et al.)",
                   "Oracle"
                   )


make_boxplot <- function(plot_setting=1, plot_model="grf", plot_n_train=500, plot_n_cal=500) {
    summary <- results %>%
        filter(setting==plot_setting, surv_model_type==plot_model, cens_model_type==plot_model, n_train==plot_n_train, n_cal==plot_n_cal) %>%
        pivot_longer(c("Coverage (event time)", "Mean lower bound"), names_to="metric", values_to="value") %>%
        mutate(metric = factor(metric, c("Coverage (event time)", "Mean lower bound"), c("Coverage", "Lower Bound"))) %>%
        filter(Method %in% method.values) %>%
        mutate(Method = fct_rev(factor(Method, levels = method.values, labels = method.labels)))
    df.limits <- tibble(metric=c("Coverage","Coverage", "Lower Bound"), value=c(0,1,0), Method="Oracle")
    ## Plotting the data
    pp <- summary %>%
        ggplot(aes(y = Method, x = value, color=Method)) +
        geom_boxplot() +
        facet_grid(.~metric, scales = "free") +
        scale_color_manual(values = method_colors) + # Apply custom color scale mapped to Method labels
        geom_vline(data = summary %>% filter(metric=="Coverage"), aes(xintercept = 0.9), linetype = "dashed", color="red") +
        geom_point(data=df.limits, alpha=0) +
        xlab("") +
        theme_bw() +
        theme(legend.position = "none")
    ggsave(sprintf("figures/experiments_boxplots_setting%d_%s_n%d_n%d.png",
                   plot_setting, plot_model, plot_n_train, plot_n_cal), plot = pp, width = 6, height = 2.5, units = "in")
}

results <- load_data(1)

for(pps in 1:4) {
##    make_boxplot(plot_setting=pps, plot_n_train=500, plot_n_cal=500)
    make_boxplot(plot_setting=pps, plot_model="cox", plot_n_train=1500, plot_n_cal=1500)
    make_boxplot(plot_setting=pps, plot_model="grf", plot_n_train=1500, plot_n_cal=1500)
}



# Custom shape mapping (optional)
method_shapes <- c(
  "Oracle" = 16,  # Circle
  "Nominal" = 17, # Triangle
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
        geom_errorbar(aes(ymin = value - SE, ymax = value + SE), width = 0.2) +
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
