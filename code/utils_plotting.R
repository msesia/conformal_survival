suppressMessages(library(tidyverse))

## Function to plot overlayed survival curves
plot_survival_curves <- function(fit, newX=NULL, individuals = NULL) {
    ## Predict survival curves
    predictions <- predict(fit, newdata=newX)
    survival_curves <- as.data.frame(predictions$predictions)
    time_points <- fit$failure.times

    df <- survival_curves
    colnames(df) <- time_points
    df$individual <- seq(1, nrow(df))
    df <- df %>%
        pivot_longer(cols=-individual,
                     names_to = "time",
                     values_to = "survival_probability") |>
    mutate(time=as.numeric(time))

    ## Filter for specific individuals if provided
    if (!is.null(individuals)) {
        df <- df %>% filter(individual %in% individuals)
    }

    ## Plot using ggplot
    df |>
    ggplot(aes(x = time,
               y = survival_probability, color = factor(individual))) +
        geom_line(linewidth = 1) +
        labs(title = "Predicted Survival Curves",
             x = "Time",
             y = "Survival Probability",
             color = "Individual") +
        theme_minimal()
}
