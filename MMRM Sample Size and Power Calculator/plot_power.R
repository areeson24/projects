# Function to plot sample size / power results for trial scenarios
plot_power <- function(data, xvar, yvar, pow, n,
                       xlabel, ylabel, title) {
  
  # Filter data for specified power / sample size
  if (tolower(yvar) == "ntot") {
    plot_tmp <- data |> filter(power == pow)
  }
  else if (tolower(yvar) == "power") {
    plot_tmp <- data |> filter(ntot == n)
  }
  
  # Create bar plot
  plot <- plot_tmp |>
    ggplot(aes(x = as.factor(!!sym(xvar)), y = !!sym(yvar), 
               fill = as.factor(dropout), group = as.factor(dropout))) +
    geom_bar(stat = "identity", color = "black",
             position = "dodge") +
    geom_text(aes(label = !!sym(yvar)), 
              position = position_dodge(width = 0.9), 
              vjust = -0.3, size = 3, fontface = "bold") +
    labs(x = xlabel, y = ylabel, title = title,
         fill = "Exponential dropout rate") +
    scale_y_continuous(expand = c(0, 0), 
                       limits = c(0, 1.1*max(data[[yvar]]))) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  # Print plot
  print(plot)

}