# LTMLE implementation: 
# 1. NO CV-TMLE. 
# 2. Only the abar but no deterministic.g.function

library(lmtp)
library(ltmle)
library(data.table)
library(SuperLearner)

# ============= DATA GENERATION FUNCTIONS =============

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
      # Handle intervention
      if(is.function(intervene_A)) {
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
    
    # Censoring
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

# Fit propensity score models
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

# Static shift function
static_shift <- function(data, trt) {
  data <- as.data.frame(data)
  return(ifelse(is.na(data[[trt]]), NA, 1))
}

# Create shift function from fitted models - UPDATED WITH alpha_threshold
create_shift_from_models <- function(g_models, shift_type = "simple", alpha_threshold = 0.001) {
  
  if(shift_type == "simple") {
    shift_func <- function(data, trt) {
      data <- as.data.frame(data)
      n <- nrow(data)
      t_current <- as.numeric(gsub("A", "", trt)) 
      
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
      
      result <- ifelse(cum_prod > alpha_threshold, 1, data[[trt]])
      return(ifelse(is.na(result), data[[trt]], result))
    }
    
  } else if(shift_type == "markmaya") {
    shift_func <- function(data, trt) {
      data <- as.data.frame(data)
      n <- nrow(data)
      t_current <- as.numeric(gsub("A", "", trt))
      
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
      
      if(t_current == 1) {
        return(ifelse(g[[1]] > alpha_threshold, 1, data[[trt]]))
      }
      
      cum_prod_hist <- rep(1, n)
      
      for(t in 1:(t_current-1)) {
        if(t == 1) {
          cum_prod_hist <- ifelse(g[[1]] > alpha_threshold, g[[1]], cum_prod_hist)
        } else {
          cum_prod_hist <- ifelse(!is.na(g[[t]]) & g[[t]] * cum_prod_hist > alpha_threshold, 
                                  cum_prod_hist * g[[t]], 
                                  cum_prod_hist)
        }
      }
      
      result <- ifelse(!is.na(g[[t_current]]) & g[[t_current]] * cum_prod_hist > alpha_threshold, 
                       1, data[[trt]])
      return(ifelse(is.na(result), data[[trt]], result))
    }
  }
  
  return(shift_func)
}

# Calculate truth using Monte Carlo - UPDATED WITH alpha_threshold
calculate_truth_with_models <- function(g_models, shift_type = "static", alpha_threshold = 0.001, 
                                        mc_size = 50000) {
  
  if(shift_type == "static") {
    data_truth <- genData_10tp_modified(n = mc_size, intervene_A = 1, intervene_C = 1)
    return(mean(data_truth$Y11, na.rm = TRUE))
  }
  
  shift_func <- create_shift_from_models(g_models, shift_type = shift_type, alpha_threshold = alpha_threshold)
  data_truth <- genData_10tp_modified(n = mc_size, 
                                      intervene_A = shift_func, 
                                      intervene_C = 1)
  
  return(mean(data_truth$Y11, na.rm = TRUE))
}

# ============= LTMLE-BASED TMLE IMPLEMENTATION (NO CROSS-FITTING) =============

# Function to create shifted data for ltmle
create_shifted_data_ltmle <- function(data, shift_func, tau = 10) {
  data_shifted <- data
  for(t in 1:tau) {
    A_col <- paste0("A", t)
    data_shifted[[A_col]] <- shift_func(data, A_col)
  }
  return(data_shifted)
}

# LTMLE-based TMLE WITHOUT cross-fitting - UPDATED WITH alpha_threshold
tmle_ltmle <- function(data, shift_type = "static", alpha_threshold = 0.001,
                       tau = 10, mc_size = 50000,
                       gbounds = c(0.01, 0.99)) {
  
  n <- nrow(data)
  
  cat("Processing TMLE without cross-fitting...\n")
  
  # Fit propensity models on FULL data
  if(shift_type == "static") {
    shift_func <- static_shift
  } else {
    g_models <- fit_propensity_models(data, tau = tau)
    shift_func <- create_shift_from_models(g_models, 
                                           shift_type = shift_type, 
                                           alpha_threshold = alpha_threshold)
  }
  
  # Prepare data for ltmle format
  # Order should be L1, A1, C1, L2, A2, C2, ..., L10, A10, C10, Y
  col_order <- c()
  for(t in 1:tau) {
    col_order <- c(col_order, paste0("L", t), paste0("A", t), paste0("C", t))
  }
  col_order <- c(col_order, "Y11")
  
  # Reorder data
  data_ordered <- data[, col_order]
  
  # Create shifted version of data
  data_shifted <- create_shifted_data_ltmle(data, shift_func, tau)
  print(data_shifted)
  data_shifted_ordered <- data_shifted[, col_order]
  
  # Define node types for ltmle
  Anodes <- paste0("A", 1:tau)
  Cnodes <- paste0("C", 1:tau)
  Lnodes <- paste0("L", 1:tau)
  Ynodes <- "Y11"
  
  # Create abar matrix (shifted treatments) - must be a matrix, not data.frame
  abar_df <- data_shifted_ordered[, Anodes, drop = FALSE]
  abar <- as.matrix(abar_df)
  
  # SuperLearner library
  SL.library=list("SL.glm","SL.mean")
  attr(SL.library, "return.fit") == TRUE
  
  # Run ltmle on FULL data
  result <- tryCatch({
    ltmle_fit <- ltmle(
      data = data_ordered,
      Anodes = Anodes,
      Cnodes = Cnodes,
      Lnodes = Lnodes,
      Ynodes = Ynodes,
      abar = abar,
      gbounds = gbounds,
      SL.library = SL.library,
      estimate.time = FALSE,
      observation.weights = NULL,
      variance.method = "ic"
    )
    
    # Extract estimate and IC
    psi <- ltmle_fit$estimates["tmle"]
    IC <- ltmle_fit$IC$tmle
    
    # Handle potential NULL IC
    if(is.null(IC)) {
      IC <- rep(NA, 0)
    }
    
    # Calculate standard error
    if(length(IC) > 0 && !all(is.na(IC))) {
      se <- sd(IC) / sqrt(n)
    } else {
      se <- NA
    }
    
    # Calculate truth
    if(shift_type == "static") {
      truth <- calculate_truth_with_models(NULL, shift_type = "static", 
                                           mc_size = mc_size)
    } else {
      truth <- calculate_truth_with_models(g_models, 
                                           shift_type = shift_type,
                                           alpha_threshold = alpha_threshold, 
                                           mc_size = mc_size)
    }
    
    list(
      estimate = psi,
      truth = truth,
      se = se,
      ci_lower = psi - 1.96 * se,
      ci_upper = psi + 1.96 * se,
      IC = IC,
      n_valid = sum(!is.na(data$Y11)),
      ltmle_fit = ltmle_fit
    )
  }, error = function(e) {
    warning(paste("LTMLE error:", e$message))
    list(
      estimate = NA,
      truth = NA,
      se = NA,
      ci_lower = NA,
      ci_upper = NA,
      IC = rep(NA, n),
      n_valid = 0
    )
  })
  
  return(result)
}

# Enhanced simulation function using LTMLE-based TMLE (NO CROSS-FITTING) - UPDATED
run_tmle_simulation <- function(n_sim = 500, n_obs = 500, tau = 10,
                                alpha_thresholds = c(0.001, 0.01),
                                mc_size = 50000,
                                gbounds = c(0.01, 0.99)) {
  
  cat("\n=== Starting LTMLE-based TMLE Simulation Study (No Cross-Fitting) ===\n")
  cat("Parameters: n_sim =", n_sim, ", n_obs =", n_obs, ", tau =", tau, "\n")
  cat("Alpha thresholds:", alpha_thresholds, "\n")
  
  # Storage for results
  results <- list()
  
  # Progress bar
  pb <- txtProgressBar(min = 0, max = n_sim, style = 3)
  
  for(sim in 1:n_sim) {
    setTxtProgressBar(pb, sim)
    
    # Generate data
    data <- genData_10tp_modified(n = n_obs)
    data <- as.data.frame(data)
    
    sim_results <- list()
    
    # 1. Static intervention with TMLE
    result_static <- tmle_ltmle(
      data = data,
      shift_type = "static",
      tau = tau,
      mc_size = mc_size,
      gbounds = gbounds
    )
    
    sim_results[["static"]] <- list(
      estimate = result_static$estimate,
      truth = result_static$truth,
      ci_lower = result_static$ci_lower,
      ci_upper = result_static$ci_upper,
      se = result_static$se
    )
    
    # 2. Simple shift with different alpha thresholds
    for(a in alpha_thresholds) {
      result_simple <- tmle_ltmle(
        data = data,
        shift_type = "simple",
        alpha_threshold = a,
        tau = tau,
        mc_size = mc_size,
        gbounds = gbounds
      )
      
      method_name <- paste0("simple_alpha_", a)
      sim_results[[method_name]] <- list(
        estimate = result_simple$estimate,
        truth = result_simple$truth,
        ci_lower = result_simple$ci_lower,
        ci_upper = result_simple$ci_upper,
        se = result_simple$se
      )
    }
    
    # 3. Mark-Maya shift with different alpha thresholds
    for(a in alpha_thresholds) {
      result_mm <- tmle_ltmle(
        data = data,
        shift_type = "markmaya",
        alpha_threshold = a,
        tau = tau,
        mc_size = mc_size,
        gbounds = gbounds
      )
      
      method_name <- paste0("mm_alpha_", a)
      sim_results[[method_name]] <- list(
        estimate = result_mm$estimate,
        truth = result_mm$truth,
        ci_lower = result_mm$ci_lower,
        ci_upper = result_mm$ci_upper,
        se = result_mm$se
      )
    }
    
    results[[sim]] <- sim_results
  }
  
  close(pb)
  
  return(results)
}

# Function to summarize TMLE results (no CV)
summarize_tmle_results <- function(sim_results) {
  if(length(sim_results) == 0) {
    cat("No results to summarize\n")
    return(NULL)
  }
  
  # Get method names
  method_names <- names(sim_results[[1]])
  n_sims <- length(sim_results)
  
  # Extract values into matrices
  estimates <- sapply(method_names, function(m) {
    sapply(sim_results, function(s) s[[m]]$estimate)
  })
  
  truths <- sapply(method_names, function(m) {
    sapply(sim_results, function(s) s[[m]]$truth)
  })
  
  ci_lowers <- sapply(method_names, function(m) {
    sapply(sim_results, function(s) s[[m]]$ci_lower)
  })
  
  ci_uppers <- sapply(method_names, function(m) {
    sapply(sim_results, function(s) s[[m]]$ci_upper)
  })
  
  ses <- sapply(method_names, function(m) {
    sapply(sim_results, function(s) s[[m]]$se)
  })
  
  # Calculate metrics
  metrics <- data.frame(
    Method = method_names,
    Mean_Estimate = colMeans(estimates, na.rm = TRUE),
    Mean_Truth = colMeans(truths, na.rm = TRUE),
    SD_Truth = apply(truths, 2, sd, na.rm = TRUE),
    Bias = colMeans(estimates - truths, na.rm = TRUE),
    Variance = apply(estimates, 2, var, na.rm = TRUE),
    MSE = colMeans((estimates - truths)^2, na.rm = TRUE),
    Coverage = NA,
    Mean_SE = colMeans(ses, na.rm = TRUE),
    Empirical_SE = apply(estimates, 2, sd, na.rm = TRUE),
    CI_Width = colMeans(ci_uppers - ci_lowers, na.rm = TRUE),
    N_Valid = colSums(!is.na(estimates))
  )
  
  # Calculate coverage
  for(j in 1:length(method_names)) {
    valid_rows <- !is.na(estimates[, j]) & !is.na(truths[, j]) & 
      !is.na(ci_lowers[, j]) & !is.na(ci_uppers[, j])
    if(sum(valid_rows) > 0) {
      coverage_indicators <- (ci_lowers[valid_rows, j] <= truths[valid_rows, j]) & 
        (ci_uppers[valid_rows, j] >= truths[valid_rows, j])
      metrics$Coverage[j] <- mean(coverage_indicators)
    }
  }
  
  return(metrics)
}

# ============= TEST THE LTMLE-BASED TMLE (NO CROSS-FITTING) =============

cat("\n=== Testing LTMLE-based TMLE Implementation (No Cross-Fitting) ===\n")
cat("Running small test simulation (n_sim = 100)...\n\n")

set.seed(123)
test_results <- run_tmle_simulation(
  n_sim = 100, 
  n_obs = 500, 
  alpha_thresholds = c(0.01, 0.05),
  mc_size = 10000,  # Smaller for testing
  gbounds = c(0.001, 0.999)
)

test_summary <- summarize_tmle_results(test_results)

cat("\n=== LTMLE TMLE Test Results (No Cross-Fitting) ===\n")
print(test_summary)

# cat("\n=== Standard Error Comparison ===\n")
# cat("Mean SE (from IC): ", round(test_summary$Mean_SE, 4), "\n")
# cat("Empirical SE:     ", round(test_summary$Empirical_SE, 4), "\n")
# cat("Ratio (Mean/Emp): ", round(test_summary$Mean_SE / test_summary$Empirical_SE, 2), "\n")
# 
# # Main simulation code instructions
# cat("\n=== Ready for Full LTMLE TMLE Simulation (No Cross-Fitting) ===\n")
# cat("To run the full simulation with B=500, run:\n")
# cat("full_results <- run_tmle_simulation(\n")
# cat("  n_sim = 500,\n")
# cat("  n_obs = 500,\n")
# cat("  alpha_thresholds = c(0.001, 0.01),\n")
# cat("  mc_size = 50000,\n")
# cat("  gbounds = c(0.01, 0.99)\n")
# cat(")\n")
# cat("full_summary <- summarize_tmle_results(full_results)\n")
# cat("print(full_summary)\n")
# cat("# Save results\n")
# cat("saveRDS(full_results, 'tmle_ltmle_results_no_cv.rds')\n")
# cat("write.csv(full_summary, 'tmle_ltmle_summary_no_cv.csv')\n")