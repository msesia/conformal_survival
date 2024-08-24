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
  "CQR" = "orange",
  "Qi et al." = "red",
  "Candes (oracle)" = "darkblue",
  "Candes (prototype)" = "lightblue",
  "Gui et al. (oracle)" = "darkgreen",
  "Gui et al. (prototype)" = "green"
)

# Custom shape mapping (optional)
method_shapes <- c(
  "Oracle" = 16,  # Circle
  "Nominal" = 17, # Triangle
  "CQR" = 15,     # Square
  "Qi et al." = 15, # Square (same as CQR)
  "Candes (oracle)" = 16,  # Circle
  "Candes (prototype)" = 18, # Diamond
  "Gui et al. (oracle)" = 19, # Filled circle
  "Gui et al. (prototype)" = 18 # Diamond (same as Candes prototype)
)

results <- load_data(1)

summary <- results %>%
    pivot_longer(c("Coverage (observed time)", "Mean lower bound", "Median lower bound", "Mean lower bound (cover)", "Median lower bound (cover)",
                   "Coverage (event time)", "Oracle MSE"),
                 names_to="metric", values_to="value") %>%
    group_by(setup, setting, n_train, n_cal, alpha, Method, metric) %>%
    summarise(SE = sd(value)/sqrt(n()), value=mean(value))

method.values <- c("oracle", "nominal", "cqr", "cqr.decensor", "candes.oracle", "prototype.candes", "gui.oracle", "prototype.gui")
method.labels <- c("Oracle", "Nominal", "CQR", "Qi et al.", "Candes (oracle)", "Candes (prototype)", "Gui et al. (oracle)", "Gui et al. (prototype)")


methods.show <- c(1, 2, 3, 4, 5, 6, 7, 8)

#methods.show <- c(1, 2, 3, 4, 6, 8)

## Plotting the data
summary %>%
    filter(setting==1,
           metric != "Coverage (observed time)",
           Method %in% method.values[methods.show]) %>%
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
