# LTMLE implementation: 
# 1. NO CV-TMLE. 
# 2. abar with deterministic.g.function to handle the clever covariates properly.
# 3. With censoring.
# 4. SuperLearner for propensity score models
# 5. No Bootstrap-based variance estimation

library(lmtp)
library(ltmle)
library(data.table)
library(SuperLearner)

source("DGPs.R")
# ============= DATA GENERATION FUNCTIONS (UNCHANGED) =============
# Modified Data generation function with actual censoring mechanism
genData_10tp_modified <- WeakSignalDGP

# ============= UPDATED: Fit propensity score models using SuperLearner with use of S3 method Wrapper =============
fit_propensity_models <- function(data, tau = 10, SL.library = NULL) {
  # Default SuperLearner library if not specified
  if(is.null(SL.library)) {
    SL.library <- c("SL.glm", "SL.mean")  # Simple default for faster computation
    # For better performance, consider: c("SL.glm", "SL.gam", "SL.ranger", "SL.mean")
  }
  
  g_models <- list()
  
  for(t in 1:tau) {
    # Prepare predictor variables
    if(t == 1) {
      predictor_vars <- paste0("L", t)
    } else {
      predictor_vars <- c(paste0("L", t), paste0("A", t-1), paste0("L", t-1))
    }
    
    tryCatch({
      # Get indices of non-missing outcomes
      valid_idx <- !is.na(data[[paste0("A", t)]])
      
      if(any(valid_idx)) {
        # Prepare X and Y for SuperLearner
        X <- data[valid_idx, predictor_vars, drop = FALSE]
        Y <- data[valid_idx, paste0("A", t)]
        
        # Fit SuperLearner model
        sl_fit <- SuperLearner(
          Y = Y,
          X = X,
          family = binomial(),
          SL.library = SL.library,
          verbose = FALSE,
          cvControl = list(V = 10)  # 10-fold cross-validation
        )
        
        # Create a wrapper object that mimics GLM predict behavior
        # This ensures compatibility with existing predict calls in your code
        g_models[[paste0("A", t)]] <- list(
          sl_model = sl_fit,
          predictor_vars = predictor_vars,
          type = "SuperLearner"
        )
        
        # Add a custom predict method to maintain compatibility
        class(g_models[[paste0("A", t)]]) <- c("sl_wrapper", class(g_models[[paste0("A", t)]]))
        
      } else {
        g_models[[paste0("A", t)]] <- NULL
      }
    }, error = function(e) {
      warning(paste("Error fitting SuperLearner for A", t, ":", e$message))
      g_models[[paste0("A", t)]] <- NULL
    })
  }
  
  return(g_models)
}

# Define a predict method for our SuperLearner wrapper
# This ensures seamless integration with existing code
predict.sl_wrapper <- function(object, newdata, type = "response", ...) {
  if(is.null(object) || is.null(object$sl_model)) {
    return(rep(NA, nrow(newdata)))
  }
  
  # Extract predictor variables
  X_new <- newdata[, object$predictor_vars, drop = FALSE]
  
  # Get predictions from SuperLearner
  preds <- predict(object$sl_model, newdata = X_new, onlySL = TRUE)$pred
  
  # Return predictions (already in probability scale for binomial family)
  return(as.vector(preds))
}

# ========================== SHIFT FUNCTIONS ==========================

# Static shift function
static_shift <- function(data, trt) {
  data <- as.data.frame(data)
  return(ifelse(is.na(data[[trt]]), NA, 1))
}

# Create shift function from fitted models (UNCHANGED - works with SuperLearner via S3 method)
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
      
      result <- ifelse(cum_prod >= alpha, 1, data[[trt]])
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
        return(ifelse(g[[1]] >= alpha, 1, data[[trt]]))
      }
      
      cum_prod_hist <- rep(1, n)
      
      for(t in 1:(t_current-1)) {
        if(t == 1) {
          cum_prod_hist <- ifelse(g[[1]] >= alpha, g[[1]], cum_prod_hist)
        } else {
          cum_prod_hist <- ifelse(!is.na(g[[t]]) & g[[t]] * cum_prod_hist >= alpha, 
                                  cum_prod_hist * g[[t]], 
                                  cum_prod_hist)
        }
      }
      
      result <- ifelse(!is.na(g[[t_current]]) & g[[t_current]] * cum_prod_hist >= alpha, 
                       1, data[[trt]])
      return(ifelse(is.na(result), data[[trt]], result))
    }
  }
  
  return(shift_func)
}

# Function to determine which units are intervened upon (UNCHANGED - works with SuperLearner via S3 method)
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
      intervention_status[, t_current] <- cum_prod >= alpha & !is.na(data[[paste0("A", t_current)]])
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
        intervention_status[, 1] <- !is.na(g[[1]]) & g[[1]] >= alpha
      } else {
        cum_prod_hist <- rep(1, n)
        
        for(t in 1:(t_current-1)) {
          if(t == 1) {
            cum_prod_hist <- ifelse(g[[1]] >= alpha, g[[1]], cum_prod_hist)
          } else {
            cum_prod_hist <- ifelse(!is.na(g[[t]]) & g[[t]] * cum_prod_hist >= alpha, 
                                    cum_prod_hist * g[[t]], 
                                    cum_prod_hist)
          }
        }
        
        intervention_status[, t_current] <- !is.na(g[[t_current]]) & 
          g[[t_current]] * cum_prod_hist >= alpha &
          !is.na(data[[paste0("A", t_current)]])
      }
    }
  }
  
  return(intervention_status)
}

# Calculate truth using Monte Carlo
calculate_truth_with_models <- function(g_models, shift_type = "static", alpha = 0.001, 
                                        mc_size = 50000, SL.library = NULL) {
  
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
                       gbounds = c(0.01, 0.99),
                       SL.library.g = NULL) {  # Added parameter for propensity model SL library
  
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
  
  # SuperLearner library for ltmle's internal models
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
    
    # Fit propensity models with SuperLearner
    g_models <- fit_propensity_models(data, tau = tau, SL.library = SL.library.g)
    
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
    
    # Modified deterministic.g.function that passes g-values for ALL units (UNCHANGED - works with SuperLearner via S3 method)
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
      # The choice of the values to be deterministic is very very important:
      
      # 1.
      # When we have no censoring, we still choose every individual's At to 
      # be deterministic according to the prefitted g hat if it is being intervened.
      # The reason is that in such a way we can also pool since otherwise only the undeterministic 
      # one will be used to fit the g_t. Thus, in our simulation, there maybe too few observations that
      # are really intervened at later time points and yield poor g fit and even not possible if too few.
      
      # 2.
      # When we have censoring, the choice is even more subtle. We do not only need to handle 
      # the A_t's we intervened, we also need to handle the A_t's that are censored.
      # Our valid index will return TRUE OR FALSE depending on prob1 is not NA or NA respectively.
      # So, if we use the below code:
      # prob1 <- rep(NA, n),
      # we are setting the censored one as non deterministic and the ltmle will try to learn g_t out
      # of the censored A_t's and thus creating trouble.
      # To get around it, we will set the censored A_t's to be also deterministic and give it a default value.
      # But this value, due to the censoring and corresponding filtering to do the tmle update, is really
      # not being used at all. I have tried prob1 <- rep(1, n) and prob1 <- rep(0, n) and under the same seed control,
      # they give same result. In this way, we get arround the ltmle requirment and 
      # avoid implicit fitting of g using non deterministic A_t'S.
      
      # So that all are valid and for the censored default to one.
      prob1 <- rep(1, n)
      
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
                                           mc_size = mc_size,
                                           SL.library = SL.library.g)
      
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
library(parallel)
library(pbapply)

# Modified simulation function with parallel processing
run_tmle_simulation <- function(n_sim = 500, n_obs = 500, tau = 10,
                                alphas = c(0.001, 0.01),
                                mc_size = 50000,
                                gbounds = c(0.01, 0.99),
                                SL.library.g = NULL,
                                n_cores = NULL) {  # Added n_cores parameter
  
  cat("\n=== Starting LTMLE-based TMLE Simulation Study ===\n")
  cat("Parameters: n_sim =", n_sim, ", n_obs =", n_obs, ", tau =", tau, "\n")
  
  # Determine number of cores
  if(is.null(n_cores)) {
    n_cores <- detectCores() - 1
    n_cores <- min(n_cores, n_sim)
  }
  cat("Using", n_cores, "cores for parallel processing\n")
  
  if(!is.null(SL.library.g)) {
    cat("Using SuperLearner for propensity models with library:", paste(SL.library.g, collapse = ", "), "\n")
  } else {
    cat("Using default SuperLearner library for propensity models\n")
  }
  
  # Create cluster
  cl <- makeCluster(n_cores)
  
  # Export necessary objects to cluster workers
  clusterExport(cl, c(
    "genData_10tp_modified",
    "fit_propensity_models", 
    "predict.sl_wrapper",
    "static_shift",
    "create_shift_from_models",
    "determine_intervention_status",
    "calculate_truth_with_models",
    "tmle_ltmle",
    "n_obs", "tau", "alphas", "mc_size", "gbounds", "SL.library.g"
  ), envir = environment())
  
  # Load required libraries on each worker
  clusterEvalQ(cl, {
    library(lmtp)
    library(ltmle)
    library(data.table)
    library(SuperLearner)
  })
  
  # Set random seeds for reproducibility
  clusterSetRNGStream(cl, 123)
  
  # Define function to run single simulation
  run_single_sim <- function(sim_id) {
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
      gbounds = gbounds,
      SL.library.g = SL.library.g
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
        gbounds = c(0, 1),
        SL.library.g = SL.library.g
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
        gbounds = c(0, 1),
        SL.library.g = SL.library.g
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
    
    return(sim_results)
  }
  
  # REPLACE THE FOR LOOP WITH pblapply
  # This maintains the progress bar and runs in parallel
  results <- pblapply(1:n_sim, run_single_sim, cl = cl)
  
  # Clean up cluster
  stopCluster(cl)
  
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

# cat("\n=== Testing LTMLE-based TMLE Implementation with SuperLearner ===\n")
# cat("Running small test simulation (n_sim = 5)...\n\n")
# 
# # set.seed(123)
# 
# # Define SuperLearner library for propensity models
# # Start with simple library for testing
# SL.lib.g <- c("SL.glm", "SL.mean")
# # For better performance, consider:
# # SL.lib.g <- c("SL.glm", "SL.gam", "SL.ranger", "SL.mean")
# 
# # Test run
# test_results <- run_tmle_simulation(
#   n_sim = 5, 
#   n_obs = 500, 
#   alphas = c(0, 0.01, 0.05, 1),
#   mc_size = 10000,  # Smaller for testing
#   gbounds = c(0.001, 0.999),
#   SL.library.g = SL.lib.g  # Pass the SuperLearner library
# )
# 
# test_summary <- summarize_tmle_results(test_results)
# 
# cat("\n=== LTMLE TMLE Test Results with SuperLearner ===\n")
# print(test_summary)
# 
# cat("\n=== Standard Error Comparison ===\n")
# cat("Mean SE (from IC): ", round(test_summary$Mean_SE, 4), "\n")
# cat("Empirical SE:     ", round(test_summary$Empirical_SE, 4), "\n")
# cat("Ratio (Mean/Emp): ", round(test_summary$Mean_SE / test_summary$Empirical_SE, 2), "\n")

cat("\n=== Ready for Full LTMLE TMLE Simulation ===\n")
cat("To run the full simulation with B=500 and more advanced SuperLearner library, use:\n")
cat("SL.lib.g <- c('SL.glm', 'SL.gam', 'SL.ranger', 'SL.xgboost', 'SL.mean')\n")
cat("full_results <- run_tmle_simulation(n_sim=500, n_obs=500, alphas=c(0.001, 0.01), mc_size=50000, SL.library.g=SL.lib.g)\n")

# Uncomment to run full simulation:
set.seed(09182025)
SL.lib.g <- c("SL.glm", "SL.mean")  # Or more advanced library
full_results <- run_tmle_simulation(
  n_sim = 5,
  n_obs = 500,
  alphas =  c(0, 0.005, 0.01, 0.02, 5/(sqrt(500)*log(500)), 0.05, 0.10, 0.20, 0.40, 0.60, 1),
  mc_size = 50000,
  gbounds = c(0.001, 0.999),
  SL.library.g = SL.lib.g
)
full_summary <- summarize_tmle_results(full_results)
print(full_summary)
saveRDS(full_results, 'tmle_ltmle_superlearner_results_seed091825_weak.rds')
write.csv(full_summary, 'tmle_ltmle_superlearner_summary_seed091825_weak.csv')
