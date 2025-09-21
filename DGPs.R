
# DGP with storng Signal, Effect vary from 15 to 8 as alpha increases.
StrongSignalDGP <- function(n, intervene_A = NULL, intervene_C = NULL) {
  
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
  
  # Censoring base intercepts to achieve ~10% censoring probability
  # Since we want P(C=0) ≈ 0.1, and C=0 means censored, we need P(C=1) ≈ 0.9
  # logit(0.9) ≈ 2.2
  censoring_intercepts <- rep(2.2, 10)  # This gives roughly 90% probability of NOT being censored
  
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
    
    # Censoring with actual mechanism (~10% censoring probability when not intervened)
    if(is.null(intervene_C)){
      C[[t]] = rep(NA, n)
      if(sum(continue_t) > 0) {
        # Generate censoring based on covariates and treatment history
        # P(C=1|history) means probability of NOT being censored
        if(t == 1) {
          # For first time point, simpler model
          prob_not_censored <- plogis(censoring_intercepts[t] - 
                                        0.3*A[[t]][continue_t] - 
                                        0.2*L[[t]][continue_t] + 
                                        get(paste0("U.Ct", t))[continue_t])
        } else {
          # For subsequent time points, include history
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
      # When intervening on censoring
      C[[t]] <- rep(NA, n)
      C[[t]][continue_t] <- intervene_C  # Usually set to 1 to prevent censoring
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



# DGP with weak Signal, Effect vary from to as alpha increases.
WeakSignalDGP <- function(n, intervene_A = NULL, intervene_C = NULL) {
  
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
  
  # Censoring base intercepts to achieve ~10% censoring probability
  # Since we want P(C=0) ≈ 0.1, and C=0 means censored, we need P(C=1) ≈ 0.9
  # logit(0.9) ≈ 2.2
  censoring_intercepts <- rep(2.2, 10)  # This gives roughly 90% probability of NOT being censored
  
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
    
    # Censoring with actual mechanism (~10% censoring probability when not intervened)
    if(is.null(intervene_C)){
      C[[t]] = rep(NA, n)
      if(sum(continue_t) > 0) {
        # Generate censoring based on covariates and treatment history
        # P(C=1|history) means probability of NOT being censored
        if(t == 1) {
          # For first time point, simpler model
          prob_not_censored <- plogis(censoring_intercepts[t] - 
                                        0.3*A[[t]][continue_t] - 
                                        0.2*L[[t]][continue_t] + 
                                        get(paste0("U.Ct", t))[continue_t])
        } else {
          # For subsequent time points, include history
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
      # When intervening on censoring
      C[[t]] <- rep(NA, n)
      C[[t]][continue_t] <- intervene_C  # Usually set to 1 to prevent censoring
    }
  }
  
  # Final outcome - MODIFIED SECTION
  continue_final <- !is.na(C[[10]]) & C[[10]] == 1
  
  # CHANGED: Much smaller treatment effects that sum to about 1.5
  treatment_effects <- seq(0.05, 0.25, length.out = 10)  # Sum ≈ 1.5
  
  Y11 <- rep(NA, n)
  if(sum(continue_final) > 0) {
    # CHANGED: Start from -0.75 instead of 2
    outcome_val <- -0.75 
    
    for(t in 1:10) {
      outcome_val <- outcome_val + ifelse(!is.na(A[[t]][continue_final]), 
                                          treatment_effects[t] * A[[t]][continue_final], 0)
    }
    
    # CHANGED: Reduced covariate effects to keep variance reasonable
    outcome_val <- outcome_val + 
      0.1*L[[1]][continue_final] +  # Reduced from 0.5 to 0.1
      0.05*ifelse(!is.na(L[[5]][continue_final]), L[[5]][continue_final], 0) +  # Reduced from 0.3 to 0.05
      0*ifelse(!is.na(L[[10]][continue_final]), L[[10]][continue_final], 0) +  # Reduced from 0.4 to 0
      0.1*U.Y11[continue_final]  # Scaled down the random noise
    
    Y11[continue_final] <- outcome_val
  }
  
  # [Rest of the code remains the same]
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