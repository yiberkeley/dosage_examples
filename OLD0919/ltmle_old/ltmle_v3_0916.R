# LTMLE implementation: 
# 1. NO CV-TMLE. 
# 2. abar with deterministic.g.function to handle the clever covariates properly.


library(lmtp)
library(ltmle)
library(data.table)
library(SuperLearner)

# ============= DATA GENERATION FUNCTIONS (UNCHANGED) =============
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

# Create shift function from fitted models
create_shift_from_models <- function(g_models, shift_type = "simple", alpha = 0.001) {
  
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
      
      result <- ifelse(cum_prod > alpha, 1, data[[trt]])
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
        return(ifelse(g[[1]] > alpha, 1, data[[trt]]))
      }
      
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
      
      result <- ifelse(!is.na(g[[t_current]]) & g[[t_current]] * cum_prod_hist > alpha, 
                       1, data[[trt]])
      return(ifelse(is.na(result), data[[trt]], result))
    }
  }
  
  return(shift_func)
}

# ============= NEW: Track intervention decisions separately =============

# Function to determine which units are intervened upon
determine_intervention_status <- function(data, g_models, shift_type, alpha, tau = 10) {
  n <- nrow(data)
  intervention_status <- matrix(FALSE, nrow = n, ncol = tau)
  colnames(intervention_status) <- paste0("A", 1:tau)
  
  if(shift_type == "simple") {
    for(t_current in 1:tau) {
      cum_prod <- rep(1, n)
      
      for(t in 1:t_current) {
        if(!is.null(g_models[[paste0("A", t)]])) {
          g_t <- rep(1, n) 
          valid_idx <- !is.na(data[[paste0("A", t)]])
          if(any(valid_idx)) {
            g_t[valid_idx] <- predict(g_models[[paste0("A", t)]], 
                                      newdata = data[valid_idx, , drop = FALSE], 
                                      type = "response")
          }
          cum_prod <- cum_prod * g_t
        }
      }
      
      # Intervene (shift to 1) when cum_prod > alpha
      intervention_status[, t_current] <- cum_prod > alpha & !is.na(data[[paste0("A", t_current)]])
    }
    
  } else if(shift_type == "markmaya") {
    for(t_current in 1:tau) {
      g <- vector("list", t_current)
      
      for(t in 1:t_current) {
        if(!is.null(g_models[[paste0("A", t)]])) {
          g[[t]] <- rep(NA, n)
          valid_idx <- !is.na(data[[paste0("A", t)]])
          if(any(valid_idx)) {
            g[[t]][valid_idx] <- predict(g_models[[paste0("A", t)]], 
                                         newdata = data[valid_idx, , drop = FALSE], 
                                         type = "response")
          }
        } else {
          g[[t]] <- rep(NA, n)
        }
      }
      
      if(t_current == 1) {
        intervention_status[, 1] <- !is.na(g[[1]]) & g[[1]] > alpha
      } else {
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
        
        intervention_status[, t_current] <- !is.na(g[[t_current]]) & 
          g[[t_current]] * cum_prod_hist > alpha &
          !is.na(data[[paste0("A", t_current)]])
      }
    }
  }
  
  return(intervention_status)
}

# Calculate truth using Monte Carlo
calculate_truth_with_models <- function(g_models, shift_type = "static", alpha = 0.001, 
                                        mc_size = 50000) {
  
  if(shift_type == "static") {
    data_truth <- genData_10tp_modified(n = mc_size, intervene_A = 1, intervene_C = 1)
    return(mean(data_truth$Y11, na.rm = TRUE))
  }
  
  shift_func <- create_shift_from_models(g_models, shift_type = shift_type, alpha = alpha)
  data_truth <- genData_10tp_modified(n = mc_size, 
                                      intervene_A = shift_func, 
                                      intervene_C = 1)
  
  return(mean(data_truth$Y11, na.rm = TRUE))
}

# ============= LTMLE-BASED TMLE WITH PROPER STOCHASTIC SHIFT HANDLING =============

tmle_ltmle <- function(data, shift_type = "static", alpha = 0.001,
                       tau = 10, mc_size = 50000,
                       gbounds = c(0.01, 0.99)) {
  
  n <- nrow(data)
  
  cat("Processing TMLE without cross-fitting...\n")
  
  # Prepare data for ltmle format
  col_order <- c()
  for(t in 1:tau) {
    col_order <- c(col_order, paste0("L", t), paste0("A", t), paste0("C", t))
  }
  col_order <- c(col_order, "Y11")
  
  # Reorder data
  data_ordered <- data[, col_order]
  
  # Define node types for ltmle
  Anodes <- paste0("A", 1:tau)
  Cnodes <- paste0("C", 1:tau)
  Lnodes <- paste0("L", 1:tau)
  Ynodes <- "Y11"
  
  # SuperLearner library
  SL.library <- list("SL.glm", "SL.mean")
  
  if(shift_type == "static") {
    # For static intervention, all treated
    abar <- matrix(1, nrow = n, ncol = tau)
    colnames(abar) <- Anodes
    
    # Run ltmle
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
      
      if(is.null(IC)) {
        IC <- rep(NA, n)
      }
      
      # Calculate standard error
      if(length(IC) > 0 && !all(is.na(IC))) {
        se <- sd(IC) / sqrt(n)
      } else {
        se <- NA
      }
      
      # Calculate truth
      truth <- calculate_truth_with_models(NULL, shift_type = "static", 
                                           mc_size = mc_size)
      
      list(
        estimate = psi,
        truth = truth,
        se = se,
        ci_lower = psi - 1.96 * se,
        ci_upper = psi + 1.96 * se,
        IC = IC,
        n_valid = sum(!is.na(data$Y11))
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
    
  } else {
    # For stochastic shifts
    
    # Fit propensity models
    g_models <- fit_propensity_models(data, tau = tau)
    
    # Create shift function
    shift_func <- create_shift_from_models(g_models, 
                                           shift_type = shift_type, 
                                           alpha = alpha)
    
    # Determine which units are actually intervened upon
    intervention_status <- determine_intervention_status(data, g_models, 
                                                         shift_type, alpha, tau)
    
    # Create abar matrix
    abar <- matrix(NA, nrow = n, ncol = tau)
    colnames(abar) <- Anodes
    for(t in 1:tau) {
      A_col <- paste0("A", t)
      abar[, t] <- shift_func(data, A_col)
    }
    
    # Modified deterministic.g.function that passes g-values for ALL units
    det_g_func <- function(data_ltmle, current.node, nodes) {
      # Get column name
      col_name <- names(data_ltmle)[current.node]
      
      # Only handle A nodes
      if(!grepl("^A", col_name)) {
        return(NULL)
      }
      
      # Get time index
      t_idx <- as.numeric(gsub("A", "", col_name))
      
      # Units where we DON'T intervene (keep natural value)
      not_intervened <- !intervention_status[, t_idx]
      
      # Get observed A values
      observed_A <- data_ltmle[[current.node]]
      
      # Get the number of observations
      n <- length(observed_A)
      
      # Create prob1 vector for ALL units
      prob1 <- rep(NA, n)
      
      # For deterministic units (not intervened): probability = observed value
      # (if A=1, prob(A=1)=1; if A=0, prob(A=1)=0)
      prob1[not_intervened & !is.na(observed_A)] <- observed_A[not_intervened & !is.na(observed_A)]
      
      # For non-deterministic units (intervened): use the fitted g(1|L) from our models
      intervened <- !not_intervened
      if(sum(intervened & !is.na(observed_A)) > 0) {
        # Get predictions from our pre-fitted g_models
        if(!is.null(g_models[[paste0("A", t_idx)]])) {
          # Predict g(1|L) for intervened units
          pred_idx <- which(intervened & !is.na(observed_A))
          if(length(pred_idx) > 0) {
            g_preds <- predict(g_models[[paste0("A", t_idx)]], 
                               newdata = data[pred_idx, , drop = FALSE], 
                               type = "response")
            prob1[pred_idx] <- g_preds
          }
        }
      }
      
      # Remove NA values and corresponding indices
      valid_idx <- !is.na(prob1)
      
      if(sum(valid_idx) == 0) {
        return(NULL)
      }
      
      # Return: is.deterministic for BOTH deterministic and non-deterministic units
      # but prob1 contains g-values for BOTH deterministic and non-deterministic units
      return(list(
        is.deterministic = valid_idx, 
        prob1 = prob1[valid_idx]
      ))
    }   
    # Run ltmle with deterministic.g.function
    result <- tryCatch({
      ltmle_fit <- ltmle(
        data = data_ordered,
        Anodes = Anodes,
        Cnodes = Cnodes,
        Lnodes = Lnodes,
        Ynodes = Ynodes,
        abar = abar,
        deterministic.g.function = det_g_func,
        gbounds = gbounds,
        SL.library = SL.library,
        estimate.time = FALSE,
        observation.weights = NULL,
        variance.method = "ic"
      )
      
      # Extract estimate and IC
      psi <- ltmle_fit$estimates["tmle"]
      IC <- ltmle_fit$IC$tmle
      
      if(is.null(IC)) {
        IC <- rep(NA, n)
      }
      
      # Calculate standard error
      if(length(IC) > 0 && !all(is.na(IC))) {
        se <- sd(IC) / sqrt(n)
      } else {
        se <- NA
      }
      
      # Calculate truth
      truth <- calculate_truth_with_models(g_models, 
                                           shift_type = shift_type,
                                           alpha = alpha, 
                                           mc_size = mc_size)
      
      list(
        estimate = psi,
        truth = truth,
        se = se,
        ci_lower = psi - 1.96 * se,
        ci_upper = psi + 1.96 * se,
        IC = IC,
        n_valid = sum(!is.na(data$Y11))
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
  }
  
  return(result)
}

# ============= SIMULATION FUNCTIONS =============

# Enhanced simulation function
run_tmle_simulation <- function(n_sim = 500, n_obs = 500, tau = 10,
                                alphas = c(0.001, 0.01),
                                mc_size = 50000,
                                gbounds = c(0.01, 0.99)) {
  
  cat("\n=== Starting LTMLE-based TMLE Simulation Study ===\n")
  cat("Parameters: n_sim =", n_sim, ", n_obs =", n_obs, ", tau =", tau, "\n")
  
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
    
    # 2. Simple shift with different alphas
    for(a in alphas) {
      result_simple <- tmle_ltmle(
        data = data,
        shift_type = "simple",
        alpha = a,
        tau = tau,
        mc_size = mc_size,
        gbounds = c(0.001, 0.999)
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
    
    # 3. Mark-Maya shift with different alphas
    for(a in alphas) {
      result_mm <- tmle_ltmle(
        data = data,
        shift_type = "markmaya",
        alpha = a,
        tau = tau,
        mc_size = mc_size,
        gbounds = c(0.001, 0.999)
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

# Function to summarize TMLE results
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

# ============= TEST THE IMPLEMENTATION =============

# cat("\n=== Testing LTMLE-based TMLE Implementation ===\n")
# cat("Running small test simulation (n_sim = 5)...\n\n")
# 
# set.seed(123)
# test_results <- run_tmle_simulation(
#   n_sim = 100, 
#   n_obs = 500, 
#   alphas = c(0.01, 0.05),
#   mc_size = 1000,  # Smaller for testing
#   gbounds = c(0.01, 0.99) #This is for the static, the shift intervention is set to c(0.001, 0.999) by default
# )
# 
# test_summary <- summarize_tmle_results(test_results)
# 
# cat("\n=== LTMLE TMLE Test Results ===\n")
# print(test_summary)

# cat("\n=== Standard Error Comparison ===\n")
# cat("Mean SE (from IC): ", round(test_summary$Mean_SE, 4), "\n")
# cat("Empirical SE:     ", round(test_summary$Empirical_SE, 4), "\n")
# cat("Ratio (Mean/Emp): ", round(test_summary$Mean_SE / test_summary$Empirical_SE, 2), "\n")
# 

# Main simulation code
cat("\n=== Ready for Full LTMLE TMLE Simulation ===\n")
cat("To run the full simulation with B=500, run:\n")
full_results <- run_tmle_simulation(
  n_sim = 20,
  n_obs = 500,
  alphas = c(0,1),
  mc_size = 10000,
  gbounds = c(0.001, 0.999) #This is for the static, the shift intervention is set to c(0.001, 0.999) by default
)
full_summary <- summarize_tmle_results(full_results)
print(full_summary)
# Save results
# saveRDS(full_results, 'tmle_ltmle_results.rds')
# write.csv(full_summary, 'tmle_ltmle_summary.csv')



