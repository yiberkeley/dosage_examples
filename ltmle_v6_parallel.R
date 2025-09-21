# LTMLE implementation with Bootstrap Variance Estimation - PARALLELIZED VERSION
# 1. NO CV-TMLE. 
# 2. abar with deterministic.g.function to handle the clever covariates properly.
# 3. With censoring.
# 4. SuperLearner for propensity score models
# 5. Bootstrap-based variance estimation
# 6. Parallel processing for simulations

library(lmtp)
library(ltmle)
library(data.table)
library(SuperLearner)
library(parallel)
library(pbapply)

# ============= DATA GENERATION FUNCTIONS (UNCHANGED) =============
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
  
  censoring_intercepts <- rep(2.2, 10)
  
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
      if(sum(continue_t) > 0) {
        if(t == 1) {
          prob_not_censored <- plogis(censoring_intercepts[t] - 
                                        0.3*A[[t]][continue_t] - 
                                        0.2*L[[t]][continue_t] + 
                                        get(paste0("U.Ct", t))[continue_t])
        } else {
          prob_not_censored <- plogis(censoring_intercepts[t] - 
                                        0.3*A[[t]][continue_t] - 
                                        0.2*L[[t]][continue_t] -
                                        0.15*A[[t-1]][continue_t] -
                                        0.1*L[[t-1]][continue_t] + 
                                        get(paste0("U.Ct", t))[continue_t])
        }
        C[[t]][continue_t] = rbinom(sum(continue_t), 1, prob_not_censored)
      }
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

# ============= Fit propensity score models using SuperLearner =============
fit_propensity_models <- function(data, tau = 10, SL.library = NULL) {
  if(is.null(SL.library)) {
    SL.library <- c("SL.glm", "SL.mean")
  }
  
  g_models <- list()
  
  for(t in 1:tau) {
    if(t == 1) {
      predictor_vars <- paste0("L", t)
    } else {
      predictor_vars <- c(paste0("L", t), paste0("A", t-1), paste0("L", t-1))
    }
    
    tryCatch({
      valid_idx <- !is.na(data[[paste0("A", t)]])
      
      if(any(valid_idx)) {
        X <- data[valid_idx, predictor_vars, drop = FALSE]
        Y <- data[valid_idx, paste0("A", t)]
        
        sl_fit <- SuperLearner(
          Y = Y,
          X = X,
          family = binomial(),
          SL.library = SL.library,
          verbose = FALSE,
          cvControl = list(V = 10)
        )
        
        g_models[[paste0("A", t)]] <- list(
          sl_model = sl_fit,
          predictor_vars = predictor_vars,
          type = "SuperLearner"
        )
        
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

# Define predict method for SuperLearner wrapper
predict.sl_wrapper <- function(object, newdata, type = "response", ...) {
  if(is.null(object) || is.null(object$sl_model)) {
    return(rep(NA, nrow(newdata)))
  }
  
  X_new <- newdata[, object$predictor_vars, drop = FALSE]
  preds <- predict(object$sl_model, newdata = X_new, onlySL = TRUE)$pred
  
  return(as.vector(preds))
}

# ========================== SHIFT FUNCTIONS (UNCHANGED) ==========================
static_shift <- function(data, trt) {
  data <- as.data.frame(data)
  return(ifelse(is.na(data[[trt]]), NA, 1))
}

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

# ============= NEW: Helper function to run single LTMLE =============
run_single_ltmle <- function(data, shift_type, alpha, tau, gbounds, 
                             g_models = NULL, intervention_status = NULL,
                             return_ic = FALSE) {  # Added parameter to optionally return IC
  
  n <- nrow(data)
  
  # Prepare data for ltmle format
  col_order <- c()
  for(t in 1:tau) {
    col_order <- c(col_order, paste0("L", t), paste0("A", t), paste0("C", t))
  }
  col_order <- c(col_order, "Y11")
  
  data_ordered <- data[, col_order]
  
  # Define node types
  Anodes <- paste0("A", 1:tau)
  Cnodes <- paste0("C", 1:tau)
  Lnodes <- paste0("L", 1:tau)
  Ynodes <- "Y11"
  
  # SuperLearner library for ltmle's internal models
  SL.library <- list("SL.glm", "SL.mean")
  
  if(shift_type == "static") {
    # For static intervention
    abar <- matrix(1, nrow = n, ncol = tau)
    colnames(abar) <- Anodes
    
    # Run ltmle
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
    
    if(return_ic) {
      return(list(
        estimate = ltmle_fit$estimates["tmle"],
        IC = ltmle_fit$IC$tmle
      ))
    } else {
      return(ltmle_fit$estimates["tmle"])
    }
    
  } else {
    # For stochastic shifts
    
    # Create shift function
    shift_func <- create_shift_from_models(g_models, 
                                           shift_type = shift_type, 
                                           alpha = alpha)
    
    # Create abar matrix
    abar <- matrix(NA, nrow = n, ncol = tau)
    colnames(abar) <- Anodes
    for(t in 1:tau) {
      A_col <- paste0("A", t)
      abar[, t] <- shift_func(data, A_col)
    }
    
    # Deterministic g function
    det_g_func <- function(data_ltmle, current.node, nodes) {
      col_name <- names(data_ltmle)[current.node]
      
      if(!grepl("^A", col_name)) {
        return(NULL)
      }
      
      t_idx <- as.numeric(gsub("A", "", col_name))
      not_intervened <- !intervention_status[, t_idx]
      observed_A <- data_ltmle[[current.node]]
      n <- length(observed_A)
      
      prob1 <- rep(1, n)
      prob1[not_intervened & !is.na(observed_A)] <- observed_A[not_intervened & !is.na(observed_A)]
      
      intervened <- !not_intervened
      if(sum(intervened & !is.na(observed_A)) > 0) {
        if(!is.null(g_models[[paste0("A", t_idx)]])) {
          pred_idx <- which(intervened & !is.na(observed_A))
          if(length(pred_idx) > 0) {
            g_preds <- predict(g_models[[paste0("A", t_idx)]], 
                               newdata = data[pred_idx, , drop = FALSE], 
                               type = "response")
            prob1[pred_idx] <- g_preds
          }
        }
      }
      
      valid_idx <- !is.na(prob1)
      
      if(sum(valid_idx) == 0) {
        return(NULL)
      }
      
      return(list(
        is.deterministic = valid_idx, 
        prob1 = prob1[valid_idx]
      ))
    }
    
    # Run ltmle
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
    
    if(return_ic) {
      return(list(
        estimate = ltmle_fit$estimates["tmle"],
        IC = ltmle_fit$IC$tmle
      ))
    } else {
      return(ltmle_fit$estimates["tmle"])
    }
  }
}

# ============= MODIFIED: LTMLE WITH BOTH IC AND BOOTSTRAP VARIANCE =============
tmle_ltmle_bootstrap <- function(data, shift_type = "static", alpha = 0.001,
                                 tau = 10, mc_size = 50000,
                                 gbounds = c(0.01, 0.99),
                                 SL.library.g = NULL,
                                 n_bootstrap = 200,
                                 verbose = FALSE) {  # Added verbose parameter
  
  n <- nrow(data)
  
  if(verbose) {
    cat("Processing TMLE with IC and bootstrap variance estimation...\n")
  }
  
  # Step 1: Get point estimate AND IC-based variance from original data
  if(shift_type == "static") {
    # Static intervention
    result <- tryCatch({
      ltmle_result <- run_single_ltmle(data, shift_type, alpha, tau, gbounds, 
                                       return_ic = TRUE)
      psi_original <- ltmle_result$estimate
      IC_original <- ltmle_result$IC
      
      # Calculate IC-based SE and CI
      if(!is.null(IC_original) && length(IC_original) > 0 && !all(is.na(IC_original))) {
        se_ic <- sd(IC_original) / sqrt(n)
        ci_lower_ic <- psi_original - 1.96 * se_ic
        ci_upper_ic <- psi_original + 1.96 * se_ic
      } else {
        se_ic <- NA
        ci_lower_ic <- NA
        ci_upper_ic <- NA
      }
      
      # Calculate truth
      truth <- calculate_truth_with_models(NULL, shift_type = "static", 
                                           mc_size = mc_size)
      
      list(psi = psi_original, truth = truth,
           se_ic = se_ic, ci_lower_ic = ci_lower_ic, ci_upper_ic = ci_upper_ic)
    }, error = function(e) {
      warning(paste("LTMLE error in original data:", e$message))
      list(psi = NA, truth = NA, se_ic = NA, ci_lower_ic = NA, ci_upper_ic = NA)
    })
    
  } else {
    # Stochastic shifts - fit models on original data
    g_models_original <- fit_propensity_models(data, tau = tau, SL.library = SL.library.g)
    intervention_status_original <- determine_intervention_status(data, g_models_original, 
                                                                  shift_type, alpha, tau)
    
    result <- tryCatch({
      ltmle_result <- run_single_ltmle(data, shift_type, alpha, tau, gbounds,
                                       g_models_original, intervention_status_original,
                                       return_ic = TRUE)
      psi_original <- ltmle_result$estimate
      IC_original <- ltmle_result$IC
      
      # Calculate IC-based SE and CI
      if(!is.null(IC_original) && length(IC_original) > 0 && !all(is.na(IC_original))) {
        se_ic <- sd(IC_original) / sqrt(n)
        ci_lower_ic <- psi_original - 1.96 * se_ic
        ci_upper_ic <- psi_original + 1.96 * se_ic
      } else {
        se_ic <- NA
        ci_lower_ic <- NA
        ci_upper_ic <- NA
      }
      
      # Calculate truth
      truth <- calculate_truth_with_models(g_models_original, 
                                           shift_type = shift_type,
                                           alpha = alpha, 
                                           mc_size = mc_size,
                                           SL.library = SL.library.g)
      
      list(psi = psi_original, truth = truth,
           se_ic = se_ic, ci_lower_ic = ci_lower_ic, ci_upper_ic = ci_upper_ic,
           g_models = g_models_original, 
           intervention_status = intervention_status_original)
    }, error = function(e) {
      warning(paste("LTMLE error in original data:", e$message))
      list(psi = NA, truth = NA, se_ic = NA, ci_lower_ic = NA, ci_upper_ic = NA)
    })
  }
  
  if(is.na(result$psi)) {
    return(list(
      estimate = NA,
      truth = NA,
      se_ic = NA,
      ci_lower_ic = NA,
      ci_upper_ic = NA,
      se = NA,
      ci_lower_normal = NA,
      ci_upper_normal = NA,
      ci_lower_percentile = NA,
      ci_upper_percentile = NA,
      n_valid = 0
    ))
  }
  
  # Step 2: Bootstrap for variance estimation
  bootstrap_estimates <- numeric(n_bootstrap)
  
  for(b in 1:n_bootstrap) {
    # Resample with replacement
    boot_indices <- sample(1:n, n, replace = TRUE)
    boot_data <- data[boot_indices, ]
    
    if(shift_type == "static") {
      # Static intervention - no need to refit models
      boot_est <- tryCatch({
        run_single_ltmle(boot_data, shift_type, alpha, tau, gbounds, return_ic = FALSE)
      }, error = function(e) {
        NA
      })
      
    } else {
      # Stochastic shifts - refit models on bootstrap sample
      g_models_boot <- fit_propensity_models(boot_data, tau = tau, SL.library = SL.library.g)
      intervention_status_boot <- determine_intervention_status(boot_data, g_models_boot, 
                                                                shift_type, alpha, tau)
      
      boot_est <- tryCatch({
        run_single_ltmle(boot_data, shift_type, alpha, tau, gbounds,
                         g_models_boot, intervention_status_boot, return_ic = FALSE)
      }, error = function(e) {
        NA
      })
    }
    
    bootstrap_estimates[b] <- boot_est
  }
  
  # Remove NA bootstrap estimates
  valid_boot <- !is.na(bootstrap_estimates)
  if(sum(valid_boot) < 10) {
    warning("Too few valid bootstrap samples")
    return(list(
      estimate = result$psi,
      truth = result$truth,
      se_ic = result$se_ic,
      ci_lower_ic = result$ci_lower_ic,
      ci_upper_ic = result$ci_upper_ic,
      se = NA,
      ci_lower_normal = NA,
      ci_upper_normal = NA,
      ci_lower_percentile = NA,
      ci_upper_percentile = NA,
      n_valid = sum(!is.na(data$Y11))
    ))
  }
  
  # Calculate bootstrap SE
  se_boot <- sd(bootstrap_estimates[valid_boot])
  
  # Calculate confidence intervals - BOTH bootstrap methods
  # Method 1: Normal approximation
  ci_lower_normal <- result$psi - 1.96 * se_boot
  ci_upper_normal <- result$psi + 1.96 * se_boot
  
  # Method 2: Percentile bootstrap
  ci_lower_percentile <- quantile(bootstrap_estimates[valid_boot], 0.025)
  ci_upper_percentile <- quantile(bootstrap_estimates[valid_boot], 0.975)
  
  return(list(
    estimate = result$psi,
    truth = result$truth,
    se_ic = result$se_ic,                    # IC-based SE
    ci_lower_ic = result$ci_lower_ic,        # IC-based CI lower
    ci_upper_ic = result$ci_upper_ic,        # IC-based CI upper
    se = se_boot,                            # Bootstrap SE
    ci_lower_normal = ci_lower_normal,       # Bootstrap normal CI lower
    ci_upper_normal = ci_upper_normal,       # Bootstrap normal CI upper
    ci_lower_percentile = ci_lower_percentile,  # Bootstrap percentile CI lower
    ci_upper_percentile = ci_upper_percentile,  # Bootstrap percentile CI upper
    n_valid = sum(!is.na(data$Y11)),
    n_boot_valid = sum(valid_boot)
  ))
}

# ============= PARALLELIZED SIMULATION FUNCTION WITH BOOTSTRAP =============
run_tmle_simulation_bootstrap_parallel <- function(n_sim = 500, n_obs = 500, tau = 10,
                                                   alphas = c(0.001, 0.01),
                                                   mc_size = 50000,
                                                   gbounds = c(0.01, 0.99),
                                                   SL.library.g = NULL,
                                                   n_bootstrap = 200,
                                                   n_cores = NULL) {
  
  cat("\n=== Starting Parallelized LTMLE-based TMLE Simulation with Bootstrap Variance ===\n")
  cat("Parameters: n_sim =", n_sim, ", n_obs =", n_obs, ", tau =", tau, 
      ", n_bootstrap =", n_bootstrap, "\n")
  
  # Determine number of cores
  if(is.null(n_cores)) {
    n_cores <- detectCores() - 1
    n_cores <- min(n_cores, n_sim)
  }
  cat("Using", n_cores, "cores for parallel processing\n")
  
  if(!is.null(SL.library.g)) {
    cat("Using SuperLearner for propensity models with library:", 
        paste(SL.library.g, collapse = ", "), "\n")
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
    "run_single_ltmle",
    "tmle_ltmle_bootstrap",
    "n_obs", "tau", "alphas", "mc_size", "gbounds", "SL.library.g", "n_bootstrap"
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
  
  # Define function to run single simulation with bootstrap
  run_single_sim_bootstrap <- function(sim_id) {
    # Generate data
    data <- genData_10tp_modified(n = n_obs)
    data <- as.data.frame(data)
    
    sim_results <- list()
    
    # 1. Static intervention with TMLE
    result_static <- tmle_ltmle_bootstrap(
      data = data,
      shift_type = "static",
      tau = tau,
      mc_size = mc_size,
      gbounds = gbounds,
      SL.library.g = SL.library.g,
      n_bootstrap = n_bootstrap,
      verbose = FALSE  # Turn off verbose in parallel
    )
    
    sim_results[["static"]] <- list(
      estimate = result_static$estimate,
      truth = result_static$truth,
      ci_lower_ic = result_static$ci_lower_ic,           # IC-based CI
      ci_upper_ic = result_static$ci_upper_ic,           # IC-based CI
      ci_lower_normal = result_static$ci_lower_normal,   # Bootstrap normal CI
      ci_upper_normal = result_static$ci_upper_normal,   # Bootstrap normal CI
      ci_lower_percentile = result_static$ci_lower_percentile,  # Bootstrap percentile CI
      ci_upper_percentile = result_static$ci_upper_percentile,  # Bootstrap percentile CI
      se_ic = result_static$se_ic,  # IC-based SE
      se = result_static$se          # Bootstrap SE
    )
    
    # 2. Simple shift with different alphas
    for(a in alphas) {
      result_simple <- tmle_ltmle_bootstrap(
        data = data,
        shift_type = "simple",
        alpha = a,
        tau = tau,
        mc_size = mc_size,
        gbounds = c(0, 1),
        SL.library.g = SL.library.g,
        n_bootstrap = n_bootstrap,
        verbose = FALSE
      )
      
      method_name <- paste0("simple_alpha_", a)
      sim_results[[method_name]] <- list(
        estimate = result_simple$estimate,
        truth = result_simple$truth,
        ci_lower_ic = result_simple$ci_lower_ic,
        ci_upper_ic = result_simple$ci_upper_ic,
        ci_lower_normal = result_simple$ci_lower_normal,
        ci_upper_normal = result_simple$ci_upper_normal,
        ci_lower_percentile = result_simple$ci_lower_percentile,
        ci_upper_percentile = result_simple$ci_upper_percentile,
        se_ic = result_simple$se_ic,
        se = result_simple$se
      )
    }
    
    # 3. Mark-Maya shift with different alphas
    for(a in alphas) {
      result_mm <- tmle_ltmle_bootstrap(
        data = data,
        shift_type = "markmaya",
        alpha = a,
        tau = tau,
        mc_size = mc_size,
        gbounds = c(0, 1),
        SL.library.g = SL.library.g,
        n_bootstrap = n_bootstrap,
        verbose = FALSE
      )
      
      method_name <- paste0("mm_alpha_", a)
      sim_results[[method_name]] <- list(
        estimate = result_mm$estimate,
        truth = result_mm$truth,
        ci_lower_ic = result_mm$ci_lower_ic,
        ci_upper_ic = result_mm$ci_upper_ic,
        ci_lower_normal = result_mm$ci_lower_normal,
        ci_upper_normal = result_mm$ci_upper_normal,
        ci_lower_percentile = result_mm$ci_lower_percentile,
        ci_upper_percentile = result_mm$ci_upper_percentile,
        se_ic = result_mm$se_ic,
        se = result_mm$se
      )
    }
    
    return(sim_results)
  }
  
  # Run simulations in parallel with progress bar
  results <- pblapply(1:n_sim, run_single_sim_bootstrap, cl = cl)
  
  # Clean up cluster
  stopCluster(cl)
  
  return(results)
}

# ============= COMPREHENSIVE SUMMARY FUNCTION WITH ALL METRICS =============
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
  
  # Extract all CI types
  ci_lowers_ic <- sapply(method_names, function(m) {
    sapply(sim_results, function(s) s[[m]]$ci_lower_ic)
  })
  
  ci_uppers_ic <- sapply(method_names, function(m) {
    sapply(sim_results, function(s) s[[m]]$ci_upper_ic)
  })
  
  ci_lowers_normal <- sapply(method_names, function(m) {
    sapply(sim_results, function(s) s[[m]]$ci_lower_normal)
  })
  
  ci_uppers_normal <- sapply(method_names, function(m) {
    sapply(sim_results, function(s) s[[m]]$ci_upper_normal)
  })
  
  ci_lowers_percentile <- sapply(method_names, function(m) {
    sapply(sim_results, function(s) s[[m]]$ci_lower_percentile)
  })
  
  ci_uppers_percentile <- sapply(method_names, function(m) {
    sapply(sim_results, function(s) s[[m]]$ci_upper_percentile)
  })
  
  # Extract both SE types
  ses_ic <- sapply(method_names, function(m) {
    sapply(sim_results, function(s) s[[m]]$se_ic)
  })
  
  ses_boot <- sapply(method_names, function(m) {
    sapply(sim_results, function(s) s[[m]]$se)
  })
  
  # Calculate empirical SE for oracle
  empirical_ses <- apply(estimates, 2, sd, na.rm = TRUE)
  
  # Initialize comprehensive metrics dataframe
  metrics <- data.frame(
    Method = method_names,
    # Core estimation metrics
    Mean_Estimate = colMeans(estimates, na.rm = TRUE),
    Mean_Truth = colMeans(truths, na.rm = TRUE),
    SD_Truth = apply(truths, 2, sd, na.rm = TRUE),
    Bias = colMeans(estimates - truths, na.rm = TRUE),
    Variance = apply(estimates, 2, var, na.rm = TRUE),
    MSE = colMeans((estimates - truths)^2, na.rm = TRUE),
    # Coverage for all CI types
    Coverage_Oracle = NA,
    Coverage_IC = NA,
    Coverage_Normal = NA,
    Coverage_Percentile = NA,
    # Standard errors
    Mean_SE_IC = colMeans(ses_ic, na.rm = TRUE),
    Mean_SE_Boot = colMeans(ses_boot, na.rm = TRUE),
    Empirical_SE = empirical_ses,
    # SE ratios
    SE_Ratio_IC_Emp = colMeans(ses_ic, na.rm = TRUE) / empirical_ses,
    SE_Ratio_Boot_Emp = colMeans(ses_boot, na.rm = TRUE) / empirical_ses,
    # CI widths
    CI_Width_Oracle = 2 * 1.96 * empirical_ses,
    CI_Width_IC = colMeans(ci_uppers_ic - ci_lowers_ic, na.rm = TRUE),
    CI_Width_Normal = colMeans(ci_uppers_normal - ci_lowers_normal, na.rm = TRUE),
    CI_Width_Percentile = colMeans(ci_uppers_percentile - ci_lowers_percentile, na.rm = TRUE),
    # Width ratios relative to oracle
    Width_Ratio_IC_Oracle = NA,
    Width_Ratio_Normal_Oracle = NA,
    Width_Ratio_Percentile_Oracle = NA,
    # Number of valid simulations
    N_Valid = colSums(!is.na(estimates))
  )
  
  # Calculate width ratios
  metrics$Width_Ratio_IC_Oracle <- metrics$CI_Width_IC / metrics$CI_Width_Oracle
  metrics$Width_Ratio_Normal_Oracle <- metrics$CI_Width_Normal / metrics$CI_Width_Oracle
  metrics$Width_Ratio_Percentile_Oracle <- metrics$CI_Width_Percentile / metrics$CI_Width_Oracle
  
  # Calculate coverage for all CI types
  for(j in 1:length(method_names)) {
    # Oracle CI coverage (using empirical SE)
    valid_rows <- !is.na(estimates[, j]) & !is.na(truths[, j])
    if(sum(valid_rows) > 0 && !is.na(empirical_ses[j])) {
      oracle_ci_lower <- estimates[valid_rows, j] - 1.96 * empirical_ses[j]
      oracle_ci_upper <- estimates[valid_rows, j] + 1.96 * empirical_ses[j]
      coverage_indicators_oracle <- (oracle_ci_lower <= truths[valid_rows, j]) & 
        (oracle_ci_upper >= truths[valid_rows, j])
      metrics$Coverage_Oracle[j] <- mean(coverage_indicators_oracle)
    }
    
    # IC-based CI coverage
    valid_rows_ic <- !is.na(estimates[, j]) & !is.na(truths[, j]) & 
      !is.na(ci_lowers_ic[, j]) & !is.na(ci_uppers_ic[, j])
    if(sum(valid_rows_ic) > 0) {
      coverage_indicators_ic <- (ci_lowers_ic[valid_rows_ic, j] <= truths[valid_rows_ic, j]) & 
        (ci_uppers_ic[valid_rows_ic, j] >= truths[valid_rows_ic, j])
      metrics$Coverage_IC[j] <- mean(coverage_indicators_ic)
    }
    
    # Normal approximation CI coverage
    valid_rows_normal <- !is.na(estimates[, j]) & !is.na(truths[, j]) & 
      !is.na(ci_lowers_normal[, j]) & !is.na(ci_uppers_normal[, j])
    if(sum(valid_rows_normal) > 0) {
      coverage_indicators_normal <- (ci_lowers_normal[valid_rows_normal, j] <= truths[valid_rows_normal, j]) & 
        (ci_uppers_normal[valid_rows_normal, j] >= truths[valid_rows_normal, j])
      metrics$Coverage_Normal[j] <- mean(coverage_indicators_normal)
    }
    
    # Percentile CI coverage
    valid_rows_percentile <- !is.na(estimates[, j]) & !is.na(truths[, j]) & 
      !is.na(ci_lowers_percentile[, j]) & !is.na(ci_uppers_percentile[, j])
    if(sum(valid_rows_percentile) > 0) {
      coverage_indicators_percentile <- (ci_lowers_percentile[valid_rows_percentile, j] <= truths[valid_rows_percentile, j]) & 
        (ci_uppers_percentile[valid_rows_percentile, j] >= truths[valid_rows_percentile, j])
      metrics$Coverage_Percentile[j] <- mean(coverage_indicators_percentile)
    }
  }
  
  return(metrics)
}


# ============= TEST THE IMPLEMENTATION =============

# cat("\n=== Testing Parallelized LTMLE with IC, Bootstrap, and Oracle Variance Estimation ===\n")
# cat("    Calculating FOUR CI types: Oracle, IC-based, Bootstrap Normal, Bootstrap Percentile\n")
# cat("Running small test simulation (n_sim = 10, n_bootstrap = 50, n_cores = 2)...\n\n")
# 
# set.seed(123)
# 
# # Define SuperLearner library for propensity models
# SL.lib.g <- c("SL.glm", "SL.mean")
# 
# # Test run with bootstrap and parallel processing
# test_results <- run_tmle_simulation_bootstrap_parallel(
#   n_sim = 10,  # Small number for testing
#   n_obs = 500, 
#   alphas = c(0, 0.01, 0.05, 1),
#   mc_size = 10000,
#   gbounds = c(0.001, 0.999),
#   SL.library.g = SL.lib.g,
#   n_bootstrap = 50,  # Smaller for testing
#   n_cores = 2  # Use only 2 cores for testing
# )
# 
# test_summary <- summarize_tmle_results(test_results)


# Uncomment to run full simulation:
set.seed(0919)
SL.lib.g <- c("SL.glm", "SL.mean")
full_results <- run_tmle_simulation_bootstrap_parallel(
  n_sim = 150,
  n_obs = 500,
  alphas = c(0, 0.005, 0.01, 0.02, 5/(sqrt(500)*log(500)), 0.05, 0.10, 0.20, 0.40, 0.60, 1),
  mc_size = 50000,
  gbounds = c(0.001, 0.999),
  SL.library.g = SL.lib.g,
  n_bootstrap = 200,
  n_cores = NULL  # Will automatically use detectCores() - 1
)
full_summary <- summarize_tmle_results(full_results)

# Save all results
saveRDS(full_results, 'tmle_ltmle_all_ci_methods_results_seed0919.rds')
write.csv(full_summary, 'tmle_ltmle_all_ci_methods_summary_seed0919.csv')


