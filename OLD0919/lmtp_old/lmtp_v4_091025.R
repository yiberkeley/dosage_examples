library(lmtp)
library(data.table)

# Data generation function for 10 time points
genData_10tp <- function(n, intervene_A = NULL, intervene_C = NULL) {
  
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
  
  # Time points with gradually increasing propensity scores
  base_intercepts <- c(-4.2, -4.0, -3.0, -2.8, -2.0, -1.5, 0, 0.5, 1.5, 2.0)
  
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
        prob_At <- plogis(base_intercepts[t] + 0.2*L[[t]])
        A[[t]] = rbinom(n, 1, prob_At)
      } else {
        A[[t]] = rep(NA, n)
        if(sum(continue_t) > 0) {
          prob_At <- plogis(base_intercepts[t] + 
                              0.2*L[[t]][continue_t] + 
                              0.15*A[[t-1]][continue_t] + 
                              0.1*L[[t-1]][continue_t])
          A[[t]][continue_t] = rbinom(sum(continue_t), 1, prob_At)
        }
      }
    } else {
      A[[t]] <- rep(NA, n)
      A[[t]][continue_t] <- intervene_A
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

# Mark-Maya shift function with closure
create_mark_maya_shift <- function(alpha = 0.001) {
  function(data, trt) {
    data <- as.data.frame(data)
    n <- nrow(data)
    
    # Extract time point
    t_current <- as.numeric(gsub("A", "", trt))
    
    # Compute ALL propensity scores up to current time
    g <- vector("list", t_current)
    
    for(t in 1:t_current) {
      formula_str <- paste0("A", t, " ~ L", t)
      if(t > 1) {
        #OVER-SPECIFICATION
        # prev_L <- paste0("L", 1:(t-1), collapse = " + ")
        # prev_A <- paste0("A", 1:(t-1), collapse = " + ")
        # formula_str <- paste0(formula_str, " + ", prev_L, " + ", prev_A)
        
        #CORRECT-SPECIFICATION
        formula_str <- paste0("A", t, " ~ L", t, " + A", t-1, " + L", t-1)
      }
      
      tryCatch({
        if(any(!is.na(data[[paste0("A", t)]]))) {
          fit <- glm(as.formula(formula_str), 
                     data = data[!is.na(data[[paste0("A", t)]]),], 
                     family = binomial())
          g[[t]] <- rep(NA, n)
          g[[t]][!is.na(data[[paste0("A", t)]])] <- predict(fit, 
                                                            newdata = data[!is.na(data[[paste0("A", t)]]),], 
                                                            type = "response")
        } else {
          g[[t]] <- rep(NA, n)
        }
      }, error = function(e) {
        g[[t]] <- rep(NA, n)
      })
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

# Simple shift function with closure
create_simple_shift <- function(alpha = 0.001) {
  function(data, trt) {
    data <- as.data.frame(data)
    n <- nrow(data)
    
    t_current <- as.numeric(gsub("A", "", trt))
    
    # Compute cumulative product of ALL propensities
    cum_prod <- rep(1, n)
    
    for(t in 1:t_current) {
      formula_str <- paste0("A", t, " ~ L", t)
      if(t > 1) {
        #OVER-SPECIFICATION
        # prev_L <- paste0("L", 1:(t-1), collapse = " + ")
        # prev_A <- paste0("A", 1:(t-1), collapse = " + ")
        # formula_str <- paste0(formula_str, " + ", prev_L, " + ", prev_A)
        
        #CORRECT-SPECIFICATION
        formula_str <- paste0("A", t, " ~ L", t, " + A", t-1, " + L", t-1)
      }
      
      tryCatch({
        if(any(!is.na(data[[paste0("A", t)]]))) {
          fit <- glm(as.formula(formula_str), 
                     data = data[!is.na(data[[paste0("A", t)]]),], 
                     family = binomial())
          g_t <- rep(1, n)
          g_t[!is.na(data[[paste0("A", t)]])] <- predict(fit, 
                                                         newdata = data[!is.na(data[[paste0("A", t)]]),], 
                                                         type = "response")
          cum_prod <- cum_prod * g_t
        }
      }, error = function(e) {
        # Keep cum_prod unchanged
      })
    }
    
    result <- ifelse(cum_prod > alpha, 1, data[[trt]])
    return(ifelse(is.na(result), data[[trt]], result))
  }
}

# Diagnostic function to check intervention rates
check_intervention_rates <- function(data, shift_func, tau = 10) {
  A_vars <- paste0("A", 1:tau)
  rates <- numeric(tau)
  
  for(i in 1:tau) {
    original <- data[[A_vars[i]]]
    shifted <- shift_func(data, A_vars[i])
    rates[i] <- mean(shifted == 1 & !is.na(shifted), na.rm = TRUE)
  }
  
  names(rates) <- A_vars
  return(rates)
}


# Main simulation study function
run_simulation_study <- function(n_sim = 10, n_obs = 500, tau = 10,
                                 alphas = c(0.0001, 0.001, 0.01),
                                 g_bounds = c(0.999, 0.99, 0.95)) {
  
  cat("\n=== Starting Simulation Study ===\n")
  cat("Parameters: n_sim =", n_sim, ", n_obs =", n_obs, ", tau =", tau, "\n")
  
  # Get Monte Carlo truth
  cat("\nCalculating Monte Carlo truth...\n")
  data_truth <- genData_10tp(n = 50000, intervene_A = 1, intervene_C = 1)
  data_truth <- as.data.frame(data_truth)
  mc_truth <- mean(data_truth$Y11, na.rm = TRUE)
  cat("Monte Carlo Truth (A=1 for all):", round(mc_truth, 3), "\n\n")
  
  # Storage
  results <- list()
  
  # Progress bar
  pb <- txtProgressBar(min = 0, max = n_sim, style = 3)
  
  for(sim in 1:n_sim) {
    setTxtProgressBar(pb, sim)
    
    # Generate data and convert to data.frame
    data <- genData_10tp(n = n_obs)
    data <- as.data.frame(data)
    
    # Define variables for lmtp
    A <- paste0("A", 1:tau)
    L <- lapply(1:tau, function(t) paste0("L", t))
    Y <- "Y11"
    
    sim_results <- list()
    
    # 1. Static intervention with different g-bounds
    for(gb in g_bounds) {
      result_val <- tryCatch({
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
        result$estimate  # Return the theta value
      }, error = function(e) {
        NA  # Return NA on error
      })
      sim_results[[paste0("static_gb_", gb)]] <- result_val
    }
    
    # 2. Mark-Maya shift with different alphas
    for(a in alphas) {
      result_val <- tryCatch({
        shift_mm <- create_mark_maya_shift(alpha = a)
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
        result$estimate   # Return the theta value
      }, error = function(e) {
        NA  # Return NA on error
      })
      sim_results[[paste0("mm_alpha_", a)]] <- result_val
    }
    
    # 3. Simple shift with different alphas
    for(a in alphas) {
      result_val <- tryCatch({
        shift_s <- create_simple_shift(alpha = a)
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
        result$estimate   # Return the theta value
      }, error = function(e) {
        NA  # Return NA on error
      })
      sim_results[[paste0("simple_alpha_", a)]] <- result_val
    }
    
    # Store this simulation's results
    results[[sim]] <- sim_results
  }
  
  close(pb)
  
  return(list(
    truth = mc_truth,
    results = results
  ))
}


# Updated summarize results function with confidence intervals
summarize_results <- function(sim_output) {
  # Check if we have any results
  if(length(sim_output$results) == 0 || all(sapply(sim_output$results, length) == 0)) {
    cat("No results to summarize\n")
    return(NULL)
  }
  
  # Get all method names from first non-empty result
  method_names <- NULL
  for(i in 1:length(sim_output$results)) {
    if(length(sim_output$results[[i]]) > 0) {
      method_names <- names(sim_output$results[[i]])
      break
    }
  }
  
  if(is.null(method_names)) {
    cat("No valid results found\n")
    return(NULL)
  }
  
  # Create matrices for estimates and confidence intervals
  n_sims <- length(sim_output$results)
  n_methods <- length(method_names)
  
  estimates_matrix <- matrix(NA, nrow = n_sims, ncol = n_methods)
  ci_lower_matrix <- matrix(NA, nrow = n_sims, ncol = n_methods)
  ci_upper_matrix <- matrix(NA, nrow = n_sims, ncol = n_methods)
  colnames(estimates_matrix) <- method_names
  colnames(ci_lower_matrix) <- method_names
  colnames(ci_upper_matrix) <- method_names
  
  # Extract values from S4 objects
  for(i in 1:n_sims) {
    if(length(sim_output$results[[i]]) > 0) {
      for(j in 1:n_methods) {
        method <- method_names[j]
        if(method %in% names(sim_output$results[[i]])) {
          result_obj <- sim_output$results[[i]][[method]]
          if(!is.na(result_obj) && !is.null(result_obj)) {
            # Extract from S4 object
            estimates_matrix[i, j] <- result_obj@x
            ci <- result_obj@conf_int
            if(length(ci) == 2) {
              ci_lower_matrix[i, j] <- ci[1]
              ci_upper_matrix[i, j] <- ci[2]
            }
          }
        }
      }
    }
  }
  
  # Calculate coverage (proportion of CIs that contain the truth)
  coverage <- numeric(n_methods)
  ci_width <- numeric(n_methods)
  
  for(j in 1:n_methods) {
    valid_rows <- !is.na(ci_lower_matrix[, j]) & !is.na(ci_upper_matrix[, j])
    if(sum(valid_rows) > 0) {
      coverage[j] <- mean(ci_lower_matrix[valid_rows, j] <= sim_output$truth & 
                            ci_upper_matrix[valid_rows, j] >= sim_output$truth)
      ci_width[j] <- mean(ci_upper_matrix[valid_rows, j] - ci_lower_matrix[valid_rows, j])
    } else {
      coverage[j] <- NA
      ci_width[j] <- NA
    }
  }
  
  # Calculate metrics
  metrics <- data.frame(
    Method = method_names,
    Mean = colMeans(estimates_matrix, na.rm = TRUE),
    Bias = colMeans(estimates_matrix, na.rm = TRUE) - sim_output$truth,
    Variance = apply(estimates_matrix, 2, var, na.rm = TRUE),
    MSE = colMeans((estimates_matrix - sim_output$truth)^2, na.rm = TRUE),
    Coverage = coverage,
    CI_Width = ci_width,
    N_Valid = colSums(!is.na(estimates_matrix)),
    Truth = sim_output$truth
  )
  
  rownames(metrics) <- NULL
  return(metrics)
}

# ============= RUN THE ANALYSIS =============

# Test single run with diagnostics
set.seed(123)
cat("Generating test data...\n")
data <- genData_10tp(n = 500)
data <- as.data.frame(data)  # Convert to data.frame

A <- paste0("A", 1:10) #CAN BE MULTIPLE A's at each time point, handled through list as L.
L <- lapply(1:10, function(t) paste0("L", t))
Y <- "Y11"

# Check intervention rates for different shifts
cat("\n=== Checking Intervention Rates ===\n")

# Static shift
rates_static <- check_intervention_rates(data, static_shift)
cat("Static shift intervention rates:\n")
print(round(rates_static, 3))

# Mark-Maya shift with alpha = 0.01
shift_mm <- create_mark_maya_shift(alpha = 0.01)
rates_mm <- check_intervention_rates(data, shift_mm)
cat("\nMark-Maya shift (alpha=0.01) intervention rates:\n")
print(round(rates_mm, 3))

# Mark-Maya shift with alpha = 0.001
shift_mm <- create_mark_maya_shift(alpha = 0.001)
rates_mm <- check_intervention_rates(data, shift_mm)
cat("\nMark-Maya shift (alpha=0.001) intervention rates:\n")
print(round(rates_mm, 3))


# Simple shift with alpha = 0.001
shift_simple <- create_simple_shift(alpha = 0.001)
rates_simple <- check_intervention_rates(data, shift_simple)
cat("\nSimple shift (alpha=0.001) intervention rates:\n")
print(round(rates_simple, 3))



# Test individual estimators
cat("\n=== Testing Individual Estimators ===\n")

sl_lib <- c("SL.glm")

# Test static
cat("\n1. Static shift:\n")
result_static <- lmtp_tmle(
  data = data,
  trt = A,
  outcome = Y,
  time_vary = L,
  cens = NULL,
  shift = static_binary_on,
  mtp = TRUE,
  outcome_type = "continuous",
  folds = 1,
  learners_trt = sl_lib,
  learners_outcome = sl_lib,
  control = lmtp_control(.trim = 0.999)
)
print(result_static)

# Test Mark-Maya
cat("\n2. Mark-Maya shift:\n")
result_mm <- lmtp_tmle(
  data = data,
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

# Test Simple shift
cat("\n2. simple shift:\n")
result_simple <- lmtp_tmle(
  data = data,
  trt = A,
  outcome = Y,
  time_vary = L,
  cens = NULL,
  shift = shift_simple,
  mtp = TRUE,
  outcome_type = "continuous",
  folds = 1, #CV-TMLE fold number
  learners_trt = sl_lib,
  learners_outcome = sl_lib,
  control = lmtp_control(.trim = 0.999, 
                         .learners_outcome_folds = 5,
                         .learners_trt_folds = 5,
                         .return_full_fits = TRUE)
)
print(result_simple)



# Run full simulation (USE THE glm as learner for speed, but has some issues with fitted value 0 and 1 for propensity score)
sim_output <- run_simulation_study(n_sim = 5, n_obs = 500)
metrics <- summarize_results(sim_output)
print(metrics)



## Takeaway 1
#m as outcome regression, r as density ratio.

## Takeaway 2
result_simple$estimate #This is the ife object that has autodiff as operators on them based on the ICs.
result_simple$estimate-result_static$estimate

result_simple$estimate@std_error #Std, maybe round to 0 in the diplay (now the default is 3 decimal)

## Takeaway 3
#Theorem 2,Non-parametric causal effects based on longitudinal modified treatment policies, Ivan Diaz
#Point: Assume that d does not depend on P. So, that for the shift intervention based on the propensity score is not acceptable right now.


