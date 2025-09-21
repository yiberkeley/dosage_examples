library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# Fixed version that adapts to actual methods in data
plot_tmle_jittered <- function(sim_results,
                               title = "TMLE Results by Alpha Level",
                               subtitle = "Estimates (solid) vs Truth (dashed) with 95% CI",
                               y_label = "Estimate",
                               include_static = TRUE,
                               static_name = "static",
                               alpha_offset = 0.0001,
                               dodge_width = 0.08,
                               estimate_size = 4,
                               truth_size = 3,
                               error_bar_width = 0.008) {
  
  if(length(sim_results) == 0) {
    cat("No results to plot\n")
    return(NULL)
  }
  
  # Get method names
  method_names <- names(sim_results[[1]])
  
  # Extract estimates and truths
  estimates <- sapply(method_names, function(m) {
    sapply(sim_results, function(s) s[[m]]$estimate)
  })
  
  truths <- sapply(method_names, function(m) {
    sapply(sim_results, function(s) s[[m]]$truth)
  })
  
  # Calculate summary statistics
  summary_data <- data.frame(
    Method = method_names,
    Mean_Estimate = colMeans(estimates, na.rm = TRUE),
    Mean_Truth = colMeans(truths, na.rm = TRUE),
    Empirical_SE = apply(estimates, 2, sd, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  # Calculate 95% CI
  summary_data$CI_Lower <- summary_data$Mean_Estimate - 1.96 * summary_data$Empirical_SE
  summary_data$CI_Upper <- summary_data$Mean_Estimate + 1.96 * summary_data$Empirical_SE
  
  # Parse method names
  summary_data <- summary_data %>%
    mutate(
      Method_Type = case_when(
        grepl("^tmle_ltmle_superlearner_alpha_", Method) ~ "TMLE SL",
        grepl("^mm_alpha_", Method) ~ "MM",
        grepl("^simple_alpha_", Method) ~ "Simple",
        Method == static_name ~ "Static",
        TRUE ~ "Other"
      ),
      Alpha = case_when(
        Method == static_name ~ 0,
        grepl("_alpha_", Method) ~ {
          alpha_str <- str_extract(Method, "(?<=_alpha_)[0-9.]+|(?<=_alpha_)0$")
          as.numeric(alpha_str)
        },
        TRUE ~ NA_real_
      )
    )
  
  # Handle alpha_0
  summary_data$Alpha[grepl("_alpha_0$", summary_data$Method)] <- 0
  
  # Add offset for log scale
  summary_data <- summary_data %>%
    mutate(Alpha_plot = ifelse(Alpha == 0, alpha_offset, Alpha))
  
  # Filter and prepare data
  plot_data <- summary_data %>%
    filter(!is.na(Alpha)) %>%
    filter(Method_Type != "Other")
  
  # Get unique method types actually in the data
  unique_methods <- sort(unique(plot_data$Method_Type))
  n_methods <- length(unique_methods)
  
  # Create dodge positions based on actual methods
  dodge_positions <- seq(-dodge_width/2, dodge_width/2, length.out = n_methods)
  names(dodge_positions) <- unique_methods
  
  # Add dodge positions to data
  plot_data <- plot_data %>%
    mutate(Alpha_dodged = Alpha_plot * exp(dodge_positions[as.character(Method_Type)]))
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = Alpha_dodged, color = Method_Type)) +
    # Error bars for estimates
    geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper),
                  width = error_bar_width,
                  size = 0.7,
                  alpha = 0.6) +
    # Lines connecting estimates across alpha values
    geom_line(aes(y = Mean_Estimate, group = Method_Type),
              size = 0.8,
              alpha = 0.7) +
    # Lines connecting truths across alpha values (dashed)
    geom_line(aes(y = Mean_Truth, group = Method_Type),
              size = 0.6,
              linetype = "dashed",
              alpha = 0.5) +
    # Estimate points (filled circles)
    geom_point(aes(y = Mean_Estimate, fill = Method_Type),
               shape = 21,
               size = estimate_size,
               stroke = 1,
               color = "white") +
    # Truth points (open circles) - SAME X POSITION as estimates
    geom_point(aes(y = Mean_Truth),
               shape = 21,
               size = truth_size,
               fill = "white",
               stroke = 1.2) +
    # Log scale
    scale_x_log10(
      breaks = c(alpha_offset, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5),
      labels = c("0", "0.001", "0.002", "0.005", "0.01", "0.02", "0.05", "0.1", "0.2", "0.5")
    ) +
    # Color scales
    scale_color_brewer(name = "Method", palette = "Dark2") +
    scale_fill_brewer(name = "Method", palette = "Dark2") +
    # Labels
    labs(title = title,
         subtitle = subtitle,
         x = expression(alpha ~ "(log scale)"),
         y = y_label) +
    # Theme
    theme_bw() +
    theme(
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = "gray80"),
      legend.title = element_text(size = 11, face = "bold"),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray50"),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 11)
    ) +
    # Add annotation
    annotate("text",
             x = min(plot_data$Alpha_plot) * 1.5,
             y = min(c(plot_data$CI_Lower, plot_data$Mean_Truth)) - 
               0.03 * diff(range(c(plot_data$CI_Upper, plot_data$CI_Lower, plot_data$Mean_Truth))),
             label = "Filled circles = estimates; Open circles = truth values",
             size = 3,
             color = "gray40",
             fontface = "italic",
             hjust = 0)
  
  return(p)
}

# Usage:
p1 <- plot_tmle_jittered(sim_results)
print(p1)
p1

