library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# Violin plot with alpha value selection
plot_sampling_violin <- function(sim_results,
                                 alpha_values = NULL,  # e.g., c(0, 0.1, 0.01)
                                 method_types = NULL,  # e.g., c("mm", "simple")
                                 show_boxplot = TRUE,
                                 show_mean_points = TRUE,
                                 violin_alpha = 0.7,
                                 ncol = 3,
                                 common_scale = TRUE) {
  
  if(length(sim_results) == 0) {
    cat("No results to plot\n")
    return(NULL)
  }
  
  # Get all method names
  all_method_names <- names(sim_results[[1]])
  
  # Filter methods based on alpha values and method types
  method_names <- all_method_names
  
  # Filter by alpha values if specified
  if(!is.null(alpha_values)) {
    # Create patterns to match
    alpha_patterns <- paste0("_alpha_(", paste(alpha_values, collapse = "|"), ")$")
    
    # Handle alpha_0 specially
    if(0 %in% alpha_values) {
      alpha_patterns <- paste0("(_alpha_0$|", alpha_patterns, ")")
    }
    
    # Also include static if alpha = 0 is requested
    if(0 %in% alpha_values) {
      method_names <- method_names[grepl(alpha_patterns, method_names) | method_names == "static"]
    } else {
      method_names <- method_names[grepl(alpha_patterns, method_names)]
    }
  }
  
  # Filter by method types if specified
  if(!is.null(method_types)) {
    type_patterns <- paste0("^(", paste(method_types, collapse = "|"), ")_alpha_")
    method_names <- method_names[grepl(type_patterns, method_names) | method_names == "static"]
  }
  
  if(length(method_names) == 0) {
    cat("No methods match the specified criteria\n")
    return(NULL)
  }
  
  cat("Plotting methods:", paste(method_names, collapse = ", "), "\n")
  
  # Extract data for selected methods
  all_data <- list()
  
  for(method in method_names) {
    estimates <- sapply(sim_results, function(s) s[[method]]$estimate)
    truths <- sapply(sim_results, function(s) s[[method]]$truth)
    
    # Parse method name to extract alpha value for labeling
    alpha_val <- case_when(
      method == "static" ~ 0,
      grepl("_alpha_0$", method) ~ 0,
      grepl("_alpha_", method) ~ as.numeric(str_extract(method, "(?<=_alpha_)[0-9.]+"))
    )
    
    method_type <- case_when(
      method == "static" ~ "Static",
      grepl("^mm_", method) ~ "MM",
      grepl("^simple_", method) ~ "Simple",
      grepl("^tmle_ltmle_superlearner", method) ~ "TMLE SL",
      TRUE ~ "Other"
    )
    
    method_label <- paste0(method_type, " (α=", alpha_val, ")")
    
    method_data <- data.frame(
      Method = method,
      Method_Label = method_label,
      Alpha = alpha_val,
      Method_Type = method_type,
      Estimate = estimates,
      Truth = truths,
      stringsAsFactors = FALSE
    )
    
    all_data[[method]] <- method_data
  }
  
  # Combine and reshape
  plot_data <- bind_rows(all_data) %>%
    pivot_longer(cols = c(Estimate, Truth),
                 names_to = "Type",
                 values_to = "Value") %>%
    filter(!is.na(Value))
  
  # Order methods by alpha value and type for better visualization
  plot_data <- plot_data %>%
    mutate(Method_Label = factor(Method_Label, 
                                 levels = unique(Method_Label[order(Alpha, Method_Type)])))
  
  # Calculate y-axis limits if using common scale
  if(common_scale) {
    y_min <- min(plot_data$Value, na.rm = TRUE)
    y_max <- max(plot_data$Value, na.rm = TRUE)
    y_range <- y_max - y_min
    y_limits <- c(y_min - 0.05 * y_range, y_max + 0.05 * y_range)
  } else {
    y_limits <- NULL
  }
  
  # Calculate means
  mean_data <- plot_data %>%
    group_by(Method_Label, Type) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create violin plot
  p <- ggplot(plot_data, aes(x = Type, y = Value, fill = Type)) +
    geom_violin(alpha = violin_alpha, 
                scale = "width",
                trim = FALSE)
  
  # Add boxplot if requested
  if(show_boxplot) {
    p <- p + 
      geom_boxplot(width = 0.1, 
                   alpha = 0.5, 
                   outlier.size = 1,
                   outlier.alpha = 0.5,
                   show.legend = FALSE)
  }
  
  # Add mean points if requested
  if(show_mean_points) {
    p <- p +
      geom_point(data = mean_data,
                 aes(x = Type, y = Mean),
                 size = 3,
                 shape = 23,
                 fill = "white",
                 color = "black",
                 stroke = 1.5)
  }
  
  # Apply faceting and scales
  p <- p +
    facet_wrap(~Method_Label, ncol = ncol, scales = if(common_scale) "fixed" else "free_y")
  
  # Set y-limits if using common scale
  if(!is.null(y_limits)) {
    p <- p + scale_y_continuous(limits = y_limits)
  }
  
  # Styling
  p <- p +
    scale_fill_manual(values = c("Estimate" = "steelblue", "Truth" = "darkred")) +
    labs(title = "Sampling Distributions for Selected Alpha Values",
         subtitle = if(common_scale) 
           sprintf("Common scale: [%.3f, %.3f]", y_limits[1], y_limits[2]) 
         else "Free scales",
         x = "",
         y = "Value") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "gray95", color = NA),
      legend.position = "top",
      legend.title = element_blank(),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# Simplified version comparing specific alpha values side by side
plot_violin_alpha_comparison <- function(sim_results,
                                         alpha_values = c(0, 0.1),
                                         method_types = NULL,
                                         show_boxplot = TRUE) {
  
  if(length(sim_results) == 0) {
    cat("No results to plot\n")
    return(NULL)
  }
  
  # Get all method names
  all_method_names <- names(sim_results[[1]])
  
  # Filter methods
  method_names <- all_method_names
  
  # Filter by alpha values
  if(!is.null(alpha_values)) {
    alpha_patterns <- paste0("_alpha_(", paste(alpha_values, collapse = "|"), ")$")
    if(0 %in% alpha_values) {
      alpha_patterns <- paste0("(_alpha_0$|", alpha_patterns, ")")
      method_names <- method_names[grepl(alpha_patterns, method_names) | method_names == "static"]
    } else {
      method_names <- method_names[grepl(alpha_patterns, method_names)]
    }
  }
  
  # Filter by method types
  if(!is.null(method_types)) {
    type_patterns <- paste0("^(", paste(method_types, collapse = "|"), ")_alpha_")
    method_names <- method_names[grepl(type_patterns, method_names) | method_names == "static"]
  }
  
  if(length(method_names) == 0) {
    cat("No methods match the criteria\n")
    return(NULL)
  }
  
  # Extract data
  all_data <- list()
  
  for(method in method_names) {
    estimates <- sapply(sim_results, function(s) s[[method]]$estimate)
    truths <- sapply(sim_results, function(s) s[[method]]$truth)
    
    # Parse method info
    alpha_val <- case_when(
      method == "static" ~ 0,
      grepl("_alpha_0$", method) ~ 0,
      grepl("_alpha_", method) ~ as.numeric(str_extract(method, "(?<=_alpha_)[0-9.]+"))
    )
    
    method_type <- case_when(
      method == "static" ~ "Static",
      grepl("^mm_", method) ~ "MM",
      grepl("^simple_", method) ~ "Simple",
      grepl("^tmle_ltmle_superlearner", method) ~ "TMLE SL",
      TRUE ~ "Other"
    )
    
    method_data <- data.frame(
      Method_Type = method_type,
      Alpha = factor(paste0("α=", alpha_val)),
      Estimate = estimates,
      Truth = truths,
      stringsAsFactors = FALSE
    )
    
    all_data[[method]] <- method_data
  }
  
  # Combine and reshape
  plot_data <- bind_rows(all_data) %>%
    pivot_longer(cols = c(Estimate, Truth),
                 names_to = "Value_Type",
                 values_to = "Value") %>%
    filter(!is.na(Value))
  
  # Create grouped violin plot
  p <- ggplot(plot_data, aes(x = interaction(Alpha, Value_Type), y = Value, 
                             fill = Value_Type)) +
    geom_violin(alpha = 0.7, scale = "width")
  
  if(show_boxplot) {
    p <- p + 
      geom_boxplot(width = 0.1, alpha = 0.5, 
                   outlier.size = 1,
                   position = position_dodge(width = 0.9))
  }
  
  p <- p +
    facet_wrap(~Method_Type, ncol = 2, scales = "fixed") +
    scale_fill_manual(values = c("Estimate" = "steelblue", "Truth" = "darkred")) +
    labs(title = paste("Sampling Distributions: α =", paste(alpha_values, collapse = ", ")),
         x = "",
         y = "Value",
         fill = "Type") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "gray95", color = NA),
      legend.position = "top",
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
      panel.grid.major.x = element_blank()
    )
  
  return(p)
}

# Usage examples:
# Plot only alpha = 0 and 0.1
p1 <- plot_sampling_violin(sim_results,
                          alpha_values = c(0, 0.005,0.02,0.16,1))
print(p1)
