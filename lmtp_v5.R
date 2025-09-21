library(lmtp)
library(data.table)

# Modified Data generation function with specified propensity score pattern
genData_10tp_modified <- function(n, intervene_A = NULL, intervene_C = NULL) {
  
  # Exogenous variables for 10 time points
  for(i in 1:10){
    assign(paste0("U.Lt", i), rnorm(n, mean=0, sd=1))
    assign(paste0("U.At", i), rnorm(n, 0, 1))
    assign(paste0("U.Ct", i), rnorm(n, 0, 1))
  }
  U.Y11 <- rnorm(n, 0, 1)
  
  # Initialize lists
  L <- vector("list", 10)
  A <- vector("list", 10)
  C <- vector("list", 10)
  
  # Modified base intercepts to achieve desired propensity scores
  # Time 1-5: ~60%, Time 6-8: ~10%, Time 9-10: ~60%
  base_intercepts <- c(
    0.4, 0.4, 0.4, 0.4, 0.4,    # Time 1-5: high propensity (~60%)
    -2.2, -2.2, -2.2,           # Time 6-8: low propensity (~10%)
    0.4, 0.4                    # Time 9-10: high propensity (~60%)
  )
  
  for(t in 1:10) {
    if(t == 1) {
      continue_t <- rep(TRUE, n)
    } else {
      continue_t <- !is.na(C[[t-1]]) & C[[t-1]] == 1
    }
    
    # Time-varying covariates
    if(t == 1) {
      L[[t]] = rbinom(n, 1, 0.5)
    } else {
      L[[t]] = rep(NA, n)
      if(sum(continue_t) > 0) {
        L[[t]][continue_t] = rbinom(sum(continue_t), 1, 
                                    plogis(0.3*A[[t-1]][continue_t] + 
                                             0.2*L[[t-1]][continue_t] + 
                                             get(paste0("U.Lt", t))[continue_t]))
      }
    }
    
    # Treatment
    if(is.null(intervene_A)){
      if(t == 1) {
        prob_At <- plogis(base_intercepts[t] + 0.2*L[[t]] + get(paste0("U.At", t)))
        A[[t]] = rbinom(n, 1, prob_At)
      } else {
        A[[t]] = rep(NA, n)
        if(sum(continue_t) > 0) {
          prob_At <- plogis(base_intercepts[t] + 
                              0.2*L[[t]][continue_t] + 
                              0.15*A[[t-1]][continue_t] + 
                              0.1*L[[t-1]][continue_t] + 
                              get(paste0("U.At", t))[continue_t])
          A[[t]][continue_t] = rbinom(sum(continue_t), 1, prob_At)
        }
      }
    } else {
      # Handle intervention - can be a function or a constant
      if(is.function(intervene_A)) {
        # Build data frame up to current time point for the shift function
        temp_data <- data.frame(matrix(NA, nrow = n, ncol = 0))
        for(j in 1:t) {
          if(j <= t) {
            temp_data[[paste0("L", j)]] <- if(j == t) L[[j]] else L[[j]]
            if(j < t) {
              temp_data[[paste0("A", j)]] <- A[[j]]
            }
          }
        }
        A[[t]] <- rep(NA, n)
        if(t == 1) {
          temp_data$A1 <- rbinom(n, 1, plogis(base_intercepts[1] + 0.2*L[[1]] + get("U.At1")))
        } else {
          temp_data[[paste0("A", t)]] <- rep(NA, n)
          if(sum(continue_t) > 0) {
            prob_At <- plogis(base_intercepts[t] + 
                                0.2*L[[t]][continue_t] + 
                                0.15*A[[t-1]][continue_t] + 
                                0.1*L[[t-1]][continue_t] + 
                                get(paste0("U.At", t))[continue_t])
            temp_data[[paste0("A", t)]][continue_t] <- rbinom(sum(continue_t), 1, prob_At)
          }
        }
        A[[t]] <- intervene_A(temp_data, paste0("A", t))
      } else {
        A[[t]] <- rep(NA, n)
        A[[t]][continue_t] <- intervene_A
      }
    }
    
    # Censoring (set to 1 for no censoring)
    if(is.null(intervene_C)){
      C[[t]] = rep(NA, n)
      C[[t]][continue_t] = 1
    } else {
      C[[t]] <- rep(NA, n)
      C[[t]][continue_t] <- intervene_C
    }
  }
  
  # Final outcome
  continue_final <- !is.na(C[[10]]) & C[[10]] == 1
  treatment_effects <- seq(0.5, 2.0, length.out = 10)
  
  Y11 <- rep(NA, n)
  if(sum(continue_final) > 0) {
    outcome_val <- 2
    for(t in 1:10) {
      outcome_val <- outcome_val + ifelse(!is.na(A[[t]][continue_final]), 
                                          treatment_effects[t] * A[[t]][continue_final], 0)
    }
    outcome_val <- outcome_val + 0.5*L[[1]][continue_final] + 
      0.3*ifelse(!is.na(L[[5]][continue_final]), L[[5]][continue_final], 0) + 
      0.4*ifelse(!is.na(L[[10]][continue_final]), L[[10]][continue_final], 0) + 
      U.Y11[continue_final]
    Y11[continue_final] <- outcome_val
  }
  
  # Create dataframe
  ObsData <- data.frame(matrix(NA, nrow = n, ncol = 0))
  
  for(t in 1:10) {
    ObsData[[paste0("L", t)]] <- L[[t]]
    ObsData[[paste0("A", t)]] <- A[[t]]
    ObsData[[paste0("C", t)]] <- C[[t]]
  }
  ObsData$Y11 <- Y11
  
  return(ObsData)
}

# Static shift function
static_shift <- function(data, trt) {
  data <- as.data.frame(data)
  return(ifelse(is.na(data[[trt]]), NA, 1))
}

# Create shift function from fitted models
create_shift_from_models <- function(g_models, shift_type = "simple", alpha = 0.001) {
  
  if(shift_type == "simple") {
    shift_func <- function(data, trt) {
      data <- as.data.frame(data)
      n <- nrow(data)
      t_current <- as.numeric(gsub("A", "", trt))
      
      # Compute cumulative product using pre-fitted models
      cum_prod <- rep(1, n)
      
      for(t in 1:t_current) {
        if(!is.null(g_models[[paste0("A", t)]])) {
          g_t <- rep(1, n)
          g_t[!is.na(data[[paste0("A", t)]])] <- predict(g_models[[paste0("A", t)]], 
                                                         newdata = data[!is.na(data[[paste0("A", t)]]),], 
                                                         type = "response")
          cum_prod <- cum_prod * g_t
        }
      }
      
      result <- ifelse(cum_prod > alpha, 1, data[[trt]])
      return(ifelse(is.na(result), data[[trt]], result))
    }
    
  } else if(shift_type == "markmaya") {
    shift_func <- function(data, trt) {
      data <- as.data.frame(data)
      n <- nrow(data)
      t_current <- as.numeric(gsub("A", "", trt))
      
      # Compute propensity scores using pre-fitted models
      g <- vector("list", t_current)
      
      for(t in 1:t_current) {
        if(!is.null(g_models[[paste0("A", t)]])) {
          g[[t]] <- rep(NA, n)
          g[[t]][!is.na(data[[paste0("A", t)]])] <- predict(g_models[[paste0("A", t)]], 
                                                            newdata = data[!is.na(data[[paste0("A", t)]]),], 
                                                            type = "response")
        } else {
          g[[t]] <- rep(NA, n)
        }
      }
      
      # Apply Mark-Maya logic
      if(t_current == 1) {
        return(ifelse(g[[1]] > alpha, 1, data[[trt]]))
      }
      
      # Build cumulative product selectively
      cum_prod_hist <- rep(1, n)
      
      for(t in 1:(t_current-1)) {
        if(t == 1) {
          cum_prod_hist <- ifelse(g[[1]] > alpha, g[[1]], cum_prod_hist)
        } else {
          cum_prod_hist <- ifelse(!is.na(g[[t]]) & g[[t]] * cum_prod_hist > alpha, 
                                  cum_prod_hist * g[[t]], 
                                  cum_prod_hist)
        }
      }
      
      # Final decision for current time point
      result <- ifelse(!is.na(g[[t_current]]) & g[[t_current]] * cum_prod_hist > alpha, 
                       1, data[[trt]])
      return(ifelse(is.na(result), data[[trt]], result))
    }
  }
  
  return(shift_func)
}

# Fit propensity score models once
fit_propensity_models <- function(data, tau = 10) {
  g_models <- list()
  
  for(t in 1:tau) {
    formula_str <- paste0("A", t, " ~ L", t)
    if(t > 1) {
      formula_str <- paste0("A", t, " ~ L", t, " + A", t-1, " + L", t-1)
    }
    
    tryCatch({
      if(any(!is.na(data[[paste0("A", t)]]))) {
        g_models[[paste0("A", t)]] <- glm(as.formula(formula_str), 
                                          data = data[!is.na(data[[paste0("A", t)]]),], 
                                          family = binomial())
      } else {
        g_models[[paste0("A", t)]] <- NULL
      }
    }, error = function(e) {
      g_models[[paste0("A", t)]] <- NULL
    })
  }
  
  return(g_models)
}

# Enhanced diagnostic function to check both rates
check_intervention_rates_enhanced <- function(data, shift_func, tau = 10) {
  A_vars <- paste0("A", 1:tau)
  rates_after <- numeric(tau)
  rates_intervened <- numeric(tau)
  
  for(i in 1:tau) {
    original <- data[[A_vars[i]]]
    shifted <- shift_func(data, A_vars[i])
    
    # Proportion of A_t = 1 after intervention
    rates_after[i] <- mean(shifted == 1 & !is.na(shifted), na.rm = TRUE)
    
    # Proportion that were actually intervened (changed from 0 to 1)
    rates_intervened[i] <- mean(original == 0 & shifted == 1 & !is.na(shifted), na.rm = TRUE)
  }
  
  result <- data.frame(
    Time = 1:tau,
    Rate_After = rates_after,
    Rate_Intervened = rates_intervened,
    Natural_Rate = sapply(A_vars, function(a) mean(data[[a]] == 1, na.rm = TRUE))
  )
  
  return(result)
}

# Function to calculate replication-specific truth using pre-fitted models
calculate_truth_with_models <- function(g_models, shift_type = "simple", alpha = 0.001, 
                                        mc_size = 50000) {
  
  if(shift_type == "static") {
    # For static, we just need one truth calculation
    data_truth <- genData_10tp_modified(n = mc_size, intervene_A = 1, intervene_C = 1)
    return(mean(data_truth$Y11, na.rm = TRUE))
  }
  
  # Create shift function from fitted models
  shift_func <- create_shift_from_models(g_models, shift_type = shift_type, alpha = alpha)
  
  # Apply to MC sample
  data_truth <- genData_10tp_modified(n = mc_size, 
                                      intervene_A = shift_func, 
                                      intervene_C = 1)
  
  return(mean(data_truth$Y11, na.rm = TRUE))
}

# Main simulation study function with replication-specific truths
run_enhanced_simulation <- function(n_sim = 100, n_obs = 500, tau = 10,
                                    alphas = c(0.001, 0.01),
                                    g_bounds = c(0.999, 0.99)) {
  
  cat("\n=== Starting Enhanced Simulation Study ===\n")
  cat("Parameters: n_sim =", n_sim, ", n_obs =", n_obs, ", tau =", tau, "\n")
  
  # Get static truth once (doesn't depend on data)
  cat("\nCalculating Static Intervention Truth...\n")
  static_truth <- calculate_truth_with_models(NULL, shift_type = "static")
  cat("Static Truth (A=1 for all):", round(static_truth, 3), "\n\n")
  
  # Storage for results
  results <- list()
  truths <- list()
  intervention_rates <- list()
  
  # Progress bar
  pb <- txtProgressBar(min = 0, max = n_sim, style = 3)
  
  for(sim in 1:n_sim) {
    setTxtProgressBar(pb, sim)
    
    # Generate data
    data <- genData_10tp_modified(n = n_obs)
    data <- as.data.frame(data)
    
    # Fit propensity models ONCE for this sample
    g_models <- fit_propensity_models(data, tau = tau)
    
    # Define variables for lmtp
    A <- paste0("A", 1:tau)
    L <- lapply(1:tau, function(t) paste0("L", t))
    Y <- "Y11"
    
    sim_results <- list()
    sim_truths <- list()
    sim_rates <- list()
    
    # 1. Static intervention
    for(gb in g_bounds) {
      result_obj <- tryCatch({
        suppressWarnings({
          result <- lmtp_tmle(
            data = data,
            trt = A,
            outcome = Y,
            time_vary = L,
            cens = NULL,
            shift = static_shift,
            mtp = TRUE,
            outcome_type = "continuous",
            folds = 10,
            learners_trt = "SL.glm",
            learners_outcome = "SL.glm",
            control = lmtp_control(.trim = gb)
          )
        })
        list(estimate = result$estimate@x, 
             ci_lower = result$estimate@conf_int[1],
             ci_upper = result$estimate@conf_int[2])
      }, error = function(e) {
        list(estimate = NA, ci_lower = NA, ci_upper = NA)
      })
      
      method_name <- paste0("static_gb_", gb)
      sim_results[[method_name]] <- result_obj
      sim_truths[[method_name]] <- static_truth
      
      # Calculate intervention rates
      rates <- check_intervention_rates_enhanced(data, static_shift, tau)
      sim_rates[[method_name]] <- rates
    }
    
    # 2. Simple shift with different alphas
    for(a in alphas) {
      # Create shift function from pre-fitted models
      shift_s <- create_shift_from_models(g_models, shift_type = "simple", alpha = a)
      
      # Calculate replication-specific truth using the same models
      rep_truth <- calculate_truth_with_models(g_models, shift_type = "simple", 
                                               alpha = a, mc_size = 50000)
      
      result_obj <- tryCatch({
        suppressWarnings({
          result <- lmtp_tmle(
            data = data,
            trt = A,
            outcome = Y,
            time_vary = L,
            cens = NULL,
            shift = shift_s,
            mtp = TRUE,
            outcome_type = "continuous",
            folds = 10,
            learners_trt = "SL.glm",
            learners_outcome = "SL.glm",
            control = lmtp_control(.trim = 0.999)
          )
        })
        list(estimate = result$estimate@x, 
             ci_lower = result$estimate@conf_int[1],
             ci_upper = result$estimate@conf_int[2])
      }, error = function(e) {
        list(estimate = NA, ci_lower = NA, ci_upper = NA)
      })
      
      method_name <- paste0("simple_alpha_", a)
      sim_results[[method_name]] <- result_obj
      sim_truths[[method_name]] <- rep_truth
      
      # Calculate intervention rates
      rates <- check_intervention_rates_enhanced(data, shift_s, tau)
      sim_rates[[method_name]] <- rates
    }
    
    # 3. Mark-Maya shift with different alphas
    for(a in alphas) {
      # Create shift function from pre-fitted models
      shift_mm <- create_shift_from_models(g_models, shift_type = "markmaya", alpha = a)
      
      # Calculate replication-specific truth using the same models
      rep_truth <- calculate_truth_with_models(g_models, shift_type = "markmaya", 
                                               alpha = a, mc_size = 50000)
      
      result_obj <- tryCatch({
        suppressWarnings({
          result <- lmtp_tmle(
            data = data,
            trt = A,
            outcome = Y,
            time_vary = L,
            cens = NULL,
            shift = shift_mm,
            mtp = TRUE,
            outcome_type = "continuous",
            folds = 10,
            learners_trt = "SL.glm",
            learners_outcome = "SL.glm",
            control = lmtp_control(.trim = 0.999)
          )
        })
        list(estimate = result$estimate@x, 
             ci_lower = result$estimate@conf_int[1],
             ci_upper = result$estimate@conf_int[2])
      }, error = function(e) {
        list(estimate = NA, ci_lower = NA, ci_upper = NA)
      })
      
      method_name <- paste0("mm_alpha_", a)
      sim_results[[method_name]] <- result_obj
      sim_truths[[method_name]] <- rep_truth
      
      # Calculate intervention rates
      rates <- check_intervention_rates_enhanced(data, shift_mm, tau)
      sim_rates[[method_name]] <- rates
    }
    
    # Store this simulation's results
    results[[sim]] <- sim_results
    truths[[sim]] <- sim_truths
    intervention_rates[[sim]] <- sim_rates
  }
  
  close(pb)
  
  return(list(
    static_truth = static_truth,
    results = results,
    truths = truths,
    intervention_rates = intervention_rates
  ))
}

# Enhanced summarize function with replication-specific truths and coverage
summarize_enhanced_results <- function(sim_output) {
  if(length(sim_output$results) == 0) {
    cat("No results to summarize\n")
    return(NULL)
  }
  
  # Get method names
  method_names <- names(sim_output$results[[1]])
  n_sims <- length(sim_output$results)
  n_methods <- length(method_names)
  
  # Create storage matrices
  estimates_matrix <- matrix(NA, nrow = n_sims, ncol = n_methods)
  ci_lower_matrix <- matrix(NA, nrow = n_sims, ncol = n_methods)
  ci_upper_matrix <- matrix(NA, nrow = n_sims, ncol = n_methods)
  truths_matrix <- matrix(NA, nrow = n_sims, ncol = n_methods)
  colnames(estimates_matrix) <- method_names
  colnames(ci_lower_matrix) <- method_names
  colnames(ci_upper_matrix) <- method_names
  colnames(truths_matrix) <- method_names
  
  # Extract values properly
  for(i in 1:n_sims) {
    for(j in 1:n_methods) {
      method <- method_names[j]
      if(method %in% names(sim_output$results[[i]])) {
        result_obj <- sim_output$results[[i]][[method]]
        truth_val <- sim_output$truths[[i]][[method]]
        
        # Extract estimate and CI values from list
        if(!is.null(result_obj) && is.list(result_obj)) {
          estimates_matrix[i, j] <- result_obj$estimate
          ci_lower_matrix[i, j] <- result_obj$ci_lower
          ci_upper_matrix[i, j] <- result_obj$ci_upper
        }
        
        truths_matrix[i, j] <- truth_val
      }
    }
  }
  
  # Calculate metrics with replication-specific truths
  metrics <- data.frame(
    Method = method_names,
    Mean_Estimate = colMeans(estimates_matrix, na.rm = TRUE),
    Mean_Truth = colMeans(truths_matrix, na.rm = TRUE),
    SD_Truth = apply(truths_matrix, 2, sd, na.rm = TRUE),
    Bias = NA,
    Variance = apply(estimates_matrix, 2, var, na.rm = TRUE),
    MSE = NA,
    Coverage = NA,
    CI_Width = NA,
    N_Valid = colSums(!is.na(estimates_matrix))
  )
  
  # Calculate bias, MSE, and coverage using replication-specific truths
  for(j in 1:n_methods) {
    valid_rows <- !is.na(estimates_matrix[, j]) & !is.na(truths_matrix[, j])
    if(sum(valid_rows) > 0) {
      # Bias and MSE
      metrics$Bias[j] <- mean(estimates_matrix[valid_rows, j] - truths_matrix[valid_rows, j])
      metrics$MSE[j] <- mean((estimates_matrix[valid_rows, j] - truths_matrix[valid_rows, j])^2)
      
      # Coverage - proportion of CIs that contain the truth
      valid_ci_rows <- valid_rows & !is.na(ci_lower_matrix[, j]) & !is.na(ci_upper_matrix[, j])
      if(sum(valid_ci_rows) > 0) {
        coverage_indicators <- (ci_lower_matrix[valid_ci_rows, j] <= truths_matrix[valid_ci_rows, j]) & 
          (ci_upper_matrix[valid_ci_rows, j] >= truths_matrix[valid_ci_rows, j])
        metrics$Coverage[j] <- mean(coverage_indicators)
        metrics$CI_Width[j] <- mean(ci_upper_matrix[valid_ci_rows, j] - ci_lower_matrix[valid_ci_rows, j])
      }
    }
  }
  
  # Calculate average intervention rates
  avg_rates_after <- matrix(NA, nrow = n_methods, ncol = 10)
  avg_rates_intervened <- matrix(NA, nrow = n_methods, ncol = 10)
  rownames(avg_rates_after) <- method_names
  rownames(avg_rates_intervened) <- method_names
  
  for(method in method_names) {
    rates_list <- lapply(sim_output$intervention_rates, function(x) x[[method]])
    rates_df <- do.call(rbind, lapply(rates_list, function(x) {
      if(!is.null(x)) {
        return(cbind(x$Rate_After, x$Rate_Intervened))
      } else {
        return(matrix(NA, nrow = 10, ncol = 2))
      }
    }))
    
    if(nrow(rates_df) > 0) {
      # Extract rates for each time point
      for(t in 1:10) {
        idx <- seq(t, nrow(rates_df), by = 10)
        avg_rates_after[method, t] <- mean(rates_df[idx, 1], na.rm = TRUE)
        avg_rates_intervened[method, t] <- mean(rates_df[idx, 2], na.rm = TRUE)
      }
    }
  }
  
  return(list(
    metrics = metrics,
    avg_rates_after = avg_rates_after,
    avg_rates_intervened = avg_rates_intervened
  ))
}

# ============= TEST THE ENHANCED FUNCTIONS =============

# First, test the modified DGP
set.seed(123)
cat("\n=== Testing Modified DGP ===\n")
test_data <- genData_10tp_modified(n = 1000)
test_data <- as.data.frame(test_data)

# Check natural propensity scores
natural_rates <- sapply(1:10, function(t) mean(test_data[[paste0("A", t)]] == 1, na.rm = TRUE))
cat("\nNatural treatment rates by time point:\n")
print(data.frame(Time = 1:10, Natural_Rate = round(natural_rates, 3)))

# Test enhanced intervention rates with pre-fitted models
cat("\n=== Testing Enhanced Intervention Rates ===\n")

# Fit models once
g_models_test <- fit_propensity_models(test_data, tau = 10)

# Static shift
rates_static <- check_intervention_rates_enhanced(test_data, static_shift)
cat("\nStatic shift rates:\n")
print(rates_static)

# Simple shift using pre-fitted models
shift_simple <- create_shift_from_models(g_models_test, shift_type = "simple", alpha = 0.001)
rates_simple <- check_intervention_rates_enhanced(test_data, shift_simple)
cat("\nSimple shift (alpha=0.001) rates:\n")
print(rates_simple)

# Mark-Maya shift using pre-fitted models
shift_mm <- create_shift_from_models(g_models_test, shift_type = "markmaya", alpha = 0.001)
rates_mm <- check_intervention_rates_enhanced(test_data, shift_mm)
cat("\nMark-Maya shift (alpha=0.001) rates:\n")
print(rates_mm)

# Test single estimator runs
cat("\n=== Testing Individual Estimators ===\n")
sl_lib <- c("SL.glm")
A <- paste0("A", 1:10)
L <- lapply(1:10, function(t) paste0("L", t))
Y <- "Y11"

# Test static
cat("\n1. Static shift:\n")
result_static <- lmtp_tmle(
  data = test_data,
  trt = A,
  outcome = Y,
  time_vary = L,
  cens = NULL,
  shift = static_shift,
  mtp = TRUE,
  outcome_type = "continuous",
  folds = 5,
  learners_trt = sl_lib,
  learners_outcome = sl_lib,
  control = lmtp_control(.trim = 0.999)
)
print(result_static)

# Test Mark-Maya with pre-fitted models
cat("\n2. Mark-Maya shift (alpha=0.001):\n")
result_mm <- lmtp_tmle(
  data = test_data,
  trt = A,
  outcome = Y,
  time_vary = L,
  cens = NULL,
  shift = shift_mm,
  mtp = TRUE,
  outcome_type = "continuous",
  folds = 5,
  learners_trt = sl_lib,
  learners_outcome = sl_lib,
  control = lmtp_control(.trim = 0.999)
)
print(result_mm)

# Test Simple shift with pre-fitted models
cat("\n3. Simple shift (alpha=0.001):\n")
result_simple <- lmtp_tmle(
  data = test_data,
  trt = A,
  outcome = Y,
  time_vary = L,
  cens = NULL,
  shift = shift_simple,
  mtp = TRUE,
  outcome_type = "continuous",
  folds = 5,
  learners_trt = sl_lib,
  learners_outcome = sl_lib,
  control = lmtp_control(.trim = 0.999)
)
print(result_simple)

# Run a small simulation to test
cat("\n=== Running Small Test Simulation (n_sim=5) ===\n")
test_sim <- run_enhanced_simulation(n_sim = 5, n_obs = 500)
test_summary <- summarize_enhanced_results(test_sim)

cat("\n=== Simulation Results ===\n")
cat("\nPerformance Metrics:\n")
print(test_summary$metrics)

cat("\nAverage Rates After Intervention (by time point):\n")
print(round(test_summary$avg_rates_after, 3))

cat("\nAverage Rates of Actual Intervention (by time point):\n")
print(round(test_summary$avg_rates_intervened, 3))

# Main simulation (uncomment when ready for full run)
cat("\n=== Ready to run full simulation ===\n")
cat("To run the full simulation with B=100, uncomment and run:\n")
cat("full_sim <- run_enhanced_simulation(n_sim = 100, n_obs = 500)\n")
cat("full_summary <- summarize_enhanced_results(full_sim)\n")
cat("print(full_summary$metrics)\n")
cat("# Save results\n")
cat("saveRDS(full_sim, 'simulation_results_enhanced.rds')\n")
cat("write.csv(full_summary$metrics, 'simulation_metrics_enhanced.csv')\n")
cat("write.csv(full_summary$avg_rates_after, 'avg_rates_after.csv')\n")
cat("write.csv(full_summary$avg_rates_intervened, 'avg_rates_intervened.csv')\n")