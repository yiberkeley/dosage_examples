library(lmtp)
library(data.table)

# Enhanced Data generation function with 4 time points, censoring, and controlled propensity scores
# Enhanced Data generation function with controlled propensity score progression
genData <- function(n, p, intervene_A = NULL, intervene_C = NULL) {
  
  # Exogenous variables for 4 time points
  for(i in 1:4){
    assign(paste0("U.Lt", i), rnorm(n, mean=0, sd=1))
    assign(paste0("U.At", i), rnorm(n, 0, 1))
    assign(paste0("U.Ct", i), rnorm(n, 0, 1))
  }
  U.Y5 <- rnorm(n, 0, 1)  # Continuous outcome error term
  
  # Initialize lists to store variables
  L <- vector("list", 4)
  A <- vector("list", 4)
  C <- vector("list", 4)
  
  # Time point 1 - VERY LOW propensity scores (~0.01-0.02)
  L[[1]] = rbinom(n, 1, 0.5)  # Simple 50/50 split for L1
  
  if(is.null(intervene_A)){
    # Extremely low propensity scores for A1 (around 0.01-0.02)
    prob_A1 <- plogis(-4.2 + 0.2*L[[1]])  # Intercept -4.2 gives extremely low baseline probability
    A[[1]] = rbinom(n, 1, prob_A1)
  } else {
    A[[1]] <- rep(intervene_A, n)
  }
  
  # Censoring at time 1 (1 = not censored, 0 = censored)
  if(is.null(intervene_C)){
    C[[1]] = as.numeric(U.Ct1 < plogis(0.5 + 0.3*L[[1]]))
  } else {
    C[[1]] <- rep(intervene_C, n)
  }
  
  # Time point 2 - Low propensity scores (~0.04-0.06)
  continue_t2 <- C[[1]] == 1
  
  L[[2]] = ifelse(continue_t2, 
                  as.numeric(plogis(A[[1]] + L[[1]] + U.Lt2) > 0.3), 
                  NA)
  
  if(is.null(intervene_A)){
    # Very low propensity scores for A2 (around 0.04-0.06)
    prob_A2 <- ifelse(continue_t2, 
                      plogis(-3.0 + 0.2*L[[2]] + 0.2*A[[1]] + 0.1*L[[1]]), 
                      NA)
    A[[2]] = ifelse(continue_t2, rbinom(sum(continue_t2), 1, prob_A2[continue_t2]), NA)
  } else {
    A[[2]] <- ifelse(continue_t2, intervene_A, NA)
  }
  
  if(is.null(intervene_C)){
    C[[2]] = ifelse(continue_t2,
                    as.numeric(U.Ct2 < plogis(0.5 + 0.3*L[[2]] - 0.2*A[[2]] + 0.1*A[[1]])),
                    NA)
  } else {
    C[[2]] <- ifelse(continue_t2, intervene_C, NA)
  }
  
  # Time point 3 - VERY HIGH propensity scores (~0.8-0.9)
  # MODIFIED: L3 and A3 now depend only on A1, L1 (not A2, L2)
  continue_t3 <- continue_t2 & C[[2]] == 1 & !is.na(C[[2]])
  
  # L3 depends only on A1, L1, and noise (removed A2, L2 dependence)
  L[[3]] = ifelse(continue_t3,
                  as.numeric(plogis(A[[1]] + L[[1]] + U.Lt3) > 0.3),
                  NA)
  
  if(is.null(intervene_A)){
    # A3 depends only on L3, A1, L1 (removed A2, L2 dependence)
    prob_A3 <- ifelse(continue_t3,
                      plogis(2.2 + 0.4*L[[3]] + 0.3*A[[1]] + 0.2*L[[1]]),
                      NA)
    A[[3]] = ifelse(continue_t3, rbinom(sum(continue_t3), 1, prob_A3[continue_t3]), NA)
  } else {
    A[[3]] <- ifelse(continue_t3, intervene_A, NA)
  }
  
  if(is.null(intervene_C)){
    # C3 can still depend on recent variables and A1, L1
    C[[3]] = ifelse(continue_t3,
                    as.numeric(U.Ct3 < plogis(0.5 + 0.3*L[[3]] - 0.2*A[[3]] + 0.1*A[[1]])),
                    NA)
  } else {
    C[[3]] <- ifelse(continue_t3, intervene_C, NA)
  }
  
  # Time point 4 - VERY HIGH propensity scores (~0.8-0.9)
  # MODIFIED: L4 and A4 also independent of A2, L2
  continue_t4 <- continue_t3 & C[[3]] == 1 & !is.na(C[[3]])
  
  # L4 depends only on A3, L3, and noise (A2, L2 already removed from A3, L3)
  L[[4]] = ifelse(continue_t4,
                  as.numeric(plogis(A[[3]] + L[[3]] + U.Lt4) > 0.3),
                  NA)
  
  if(is.null(intervene_A)){
    # A4 depends only on L4, A3, L3 (which are already independent of A2, L2)
    prob_A4 <- ifelse(continue_t4,
                      plogis(2.5 + 0.4*L[[4]] + 0.3*A[[3]] + 0.2*L[[3]]),
                      NA)
    A[[4]] = ifelse(continue_t4, rbinom(sum(continue_t4), 1, prob_A4[continue_t4]), NA)
  } else {
    A[[4]] <- ifelse(continue_t4, intervene_A, NA)
  }
  
  if(is.null(intervene_C)){
    C[[4]] = ifelse(continue_t4,
                    as.numeric(U.Ct4 < plogis(0.5 + 0.3*L[[4]] - 0.2*A[[4]] + 0.1*A[[3]])),
                    NA)
  } else {
    C[[4]] <- ifelse(continue_t4, intervene_C, NA)
  }
  
  # Final continuous outcome Y5 (only for those not censored at time 4)
  continue_final <- continue_t4 & C[[4]] == 1 & !is.na(C[[4]])
  
  Y5 <- ifelse(continue_final,
               # Continuous outcome with treatment effects
               2 + 0.5*L[[1]] + 8*A[[1]] +
                 0.7*L[[3]] + 1.0*A[[3]] + 0.8*L[[4]] + 1.1*A[[4]] + U.Y5,
               NA)
  
  # Create dataframe with proper naming convention
  ObsData <- data.table(
    L1 = L[[1]], A1 = A[[1]], C1 = C[[1]],
    L2 = L[[2]], A2 = A[[2]], C2 = C[[2]], 
    L3 = L[[3]], A3 = A[[3]], C3 = C[[3]],
    L4 = L[[4]], A4 = A[[4]], C4 = C[[4]],
    Y5 = Y5
  )
  
  return(ObsData)
}

# Generate data
set.seed(123)
data <- genData(n=1000, p=4, intervene_A=NULL, intervene_C=1)
alpha <- 0.0013  # most likely skip the second and continue the third and fourth in the new dynamic shift 

#############################################################################
# MARK AND MAYA'S PROPOSAL: CUMULATIVE PROPENSITY INTERVENTION
#
# The intervention strategy:
# 1. Calculate propensity scores for treatment at each time point
# 2. Track cumulative product of propensity scores
# 3. Intervene (set treatment to 1) only when cumulative propensity > alpha
# 4. Skip intervention when "very unlikely to intervene" (cumulative propensity <= alpha)
# 5. Also intervene to prevent censoring (set all C nodes to 1)
#############################################################################

# Mark and Maya's cumulative propensity shift function
mark_maya_shift <- function(data, trt) {
  data <- as.data.frame(data)
  n <- nrow(data)
  
  # Compute propensity scores for all treatment time points
  tryCatch({
    fit1 <- glm(A1 ~ L1, data = data, family = binomial())
    g1 <- predict(fit1, type = "response")
    
    # Only compute later propensity scores if we have non-missing data
    if(any(!is.na(data$A2))) {
      fit2 <- glm(A2 ~ L1 + A1 + L2, data = data[!is.na(data$A2),], family = binomial())
      g2 <- rep(NA, n)
      g2[!is.na(data$A2)] <- predict(fit2, newdata = data[!is.na(data$A2),], type = "response")
    } else {
      g2 <- rep(NA, n)
    }
    
    if(any(!is.na(data$A3))) {
      fit3 <- glm(A3 ~ L1 + A1 + L2 + A2 + L3, data = data[!is.na(data$A3),], family = binomial())
      g3 <- rep(NA, n)
      g3[!is.na(data$A3)] <- predict(fit3, newdata = data[!is.na(data$A3),], type = "response")
    } else {
      g3 <- rep(NA, n)
    }
    
    if(any(!is.na(data$A4))) {
      fit4 <- glm(A4 ~ L1 + A1 + L2 + A2 + L3 + A3 + L4, data = data[!is.na(data$A4),], family = binomial())
      g4 <- rep(NA, n)
      g4[!is.na(data$A4)] <- predict(fit4, newdata = data[!is.na(data$A4),], type = "response")
    } else {
      g4 <- rep(NA, n)
    }
  }, error = function(e) {
    cat("Error in propensity score estimation:", e$message, "\n")
    return(data[[trt]])
  })
  
  # Mark and Maya's cumulative propensity intervention logic
  if (trt == "A1") {
    # At time 1: intervene if g1 > alpha
    return(ifelse(g1 > alpha, 1, data$A1))
    
  } else if (trt == "A2") {
    # At time 2: intervene if cumulative product g1*g2 > alpha
    # Skip if g1 <= alpha (very unlikely to intervene at baseline)
    cum_prod_hist <- 1
    cum_prod_hist <- ifelse(g1 > alpha, g1, cum_prod_hist)  # Start cumulative product
    result <- ifelse(!is.na(g2) & g2 * cum_prod_hist > alpha, 1, data$A2)
    
    return(ifelse(is.na(result), data$A2, result))
    
  } else if (trt == "A3") {
    # At time 3: track cumulative product through g1, g2, g3
    cum_prod_hist <- 1
    
    # Update cumulative product based on g1
    cum_prod_hist <- ifelse(g1 > alpha, cum_prod_hist * g1, cum_prod_hist)
    
    # Update cumulative product based on g2
    cum_prod_hist <- ifelse(!is.na(g2) & g2 * cum_prod_hist > alpha, cum_prod_hist * g2, cum_prod_hist)
    
    # Intervene if g3 * cumulative_product > alpha
    result <- ifelse(!is.na(g3) & g3 * cum_prod_hist > alpha, 1, data$A3)
    return(ifelse(is.na(result), data$A3, result))
    
  } else if (trt == "A4") {
    # At time 4: track cumulative product through all previous time points
    cum_prod_hist <- 1
    
    # Update cumulative product based on g1
    cum_prod_hist <- ifelse(g1 > alpha, cum_prod_hist * g1, cum_prod_hist)
    
    # Update cumulative product based on g2
    cum_prod_hist <- ifelse(!is.na(g2) & g2 * cum_prod_hist > alpha, cum_prod_hist * g2, cum_prod_hist)
    
    # Update cumulative product based on g3
    cum_prod_hist <- ifelse(!is.na(g3) & g3 * cum_prod_hist > alpha, cum_prod_hist * g3, cum_prod_hist)
    
    # Intervene if g4 * cumulative_product > alpha
    result <- ifelse(!is.na(g4) & g4 * cum_prod_hist > alpha, 1, data$A4)
    return(ifelse(is.na(result), data$A4, result))
    
  } else {
    return(data[[trt]])  # For any other variable, return as is
  }
}

# Simple cumulative propensity shift function
simple_shift <- function(data, trt) {
  data <- as.data.frame(data)
  n <- nrow(data)
  
  # Compute propensity scores for all treatment time points
  tryCatch({
    fit1 <- glm(A1 ~ L1, data = data, family = binomial())
    g1 <- predict(fit1, type = "response")
    
    # Only compute later propensity scores if we have non-missing data
    if(any(!is.na(data$A2))) {
      fit2 <- glm(A2 ~ L1 + A1 + L2, data = data[!is.na(data$A2),], family = binomial())
      g2 <- rep(NA, n)
      g2[!is.na(data$A2)] <- predict(fit2, newdata = data[!is.na(data$A2),], type = "response")
    } else {
      g2 <- rep(NA, n)
    }
    
    if(any(!is.na(data$A3))) {
      fit3 <- glm(A3 ~ L1 + A1 + L2 + A2 + L3, data = data[!is.na(data$A3),], family = binomial())
      g3 <- rep(NA, n)
      g3[!is.na(data$A3)] <- predict(fit3, newdata = data[!is.na(data$A3),], type = "response")
    } else {
      g3 <- rep(NA, n)
    }
    
    if(any(!is.na(data$A4))) {
      fit4 <- glm(A4 ~ L1 + A1 + L2 + A2 + L3 + A3 + L4, data = data[!is.na(data$A4),], family = binomial())
      g4 <- rep(NA, n)
      g4[!is.na(data$A4)] <- predict(fit4, newdata = data[!is.na(data$A4),], type = "response")
    } else {
      g4 <- rep(NA, n)
    }
  }, error = function(e) {
    cat("Error in propensity score estimation:", e$message, "\n")
    return(data[[trt]])
  })
  
  # simple cumulative propensity intervention logic
  if (trt == "A1") {
    # At time 1: intervene if g1 > alpha
    return(ifelse(g1 > alpha, 1, data$A1))
    
  } else if (trt == "A2") {
    
    result <- ifelse(!is.na(g2) & g2 *g1  > alpha, 1, data$A2)
    return(ifelse(is.na(result), data$A2, result))
  
  } else if (trt == "A3") {
   
    result <- ifelse(!is.na(g3) & g3 * g2 *g1 > alpha, 1, data$A3)
    return(ifelse(is.na(result), data$A3, result))
    
  } else if (trt == "A4") {
    
    result <- ifelse(!is.na(g4) & g4 * g3 * g2 *g1 > alpha, 1, data$A4)
    return(ifelse(is.na(result), data$A4, result))
    
  } else {
    return(data[[trt]])  # For any other variable, return as is
  }
}


# Static shift
static_shift <- function(data, trt) {
  data <- as.data.frame(data)
  return(ifelse(is.na(data[[trt]]), data[[trt]], 1))
}


# Run analysis with 4 time points, censoring, and Mark & Maya's intervention
data <- as.data.frame(data)

# Define variables for lmtp
A <- c("A1", "A2", "A3", "A4")  # Treatment variables
L <- list(c("L1"), c("L2"), c("L3"), c("L4"))  # Time-varying covariates
C <- c("C1", "C2", "C3", "C4")  # Censoring indicators
Y <- "Y5"  # Single continuous outcome

# Check propensity score distributions
cat("Checking propensity score distributions...\n")
fit1 <- glm(A1 ~ L1, data = data, family = binomial())
g1 <- predict(fit1, type = "response")

fit2 <- glm(A2 ~ L1 + A1 + L2, data = data[!is.na(data$A2),], family = binomial())
g2 <- predict(fit2, type = "response")

fit3 <- glm(A3 ~ L1 + A1 + L2 + A2 + L3, data = data[!is.na(data$A3),], family = binomial())
g3 <- predict(fit3, type = "response")

fit4 <- glm(A4 ~ L1 + A1 + L2 + A2 + L3 + A3 + L4, data = data[!is.na(data$A4),], family = binomial())
g4 <- predict(fit4, type = "response")

cat("Proportion with g1 < 0.05:", mean(g1 < 0.05), "\n")
cat("Proportion with g2 > 0.3:", mean(g2 > 0.3, na.rm = TRUE), "\n")
cat("Proportion with g3 > 0.3:", mean(g3 > 0.3, na.rm = TRUE), "\n")
cat("Proportion with g4 > 0.3:", mean(g4 > 0.3, na.rm = TRUE), "\n")

cat("Running LMTP analysis with Mark & Maya's cumulative propensity intervention...\n")
cat("Data structure:\n")
cat("- Treatment nodes:", paste(A, collapse=", "), "\n")
cat("- Censoring nodes:", paste(C, collapse=", "), "\n") 
cat("- Outcome: Y5 (continuous)\n")
cat("- Time-varying covariates: L1, L2, L3, L4\n")
cat("- Alpha threshold:", alpha, "\n\n")

# Check data structure
cat("Sample of generated data:\n")
print(head(data[, c("L1", "A1", "C1", "L2", "A2", "C2", "L3", "A3", "C3", "L4", "A4", "C4", "Y5")]))

# Run the analysis
# ?lmtp_control
result_mark_maya <- lmtp_tmle(
  data = data,
  trt = A,
  outcome = Y,
  time_vary = L,
  cens = NULL,
  shift = mark_maya_shift,
  mtp = TRUE,
  outcome_type = "continuous",  # Changed from "survival" to "continuous"
  folds = 5,
  learners_trt = c("SL.glm"),
  learners_outcome = c("SL.glm"),
  control=lmtp_control(.trim = 1)
)

print(result_mark_maya)

result_simple <- lmtp_tmle(
  data = data,
  trt = A,
  outcome = Y,
  time_vary = L,
  cens = NULL,
  shift = simple_shift,
  mtp = TRUE,
  outcome_type = "continuous",  # Changed from "survival" to "continuous"
  folds = 5,
  learners_trt = c("SL.glm"),
  learners_outcome = c("SL.glm"),
  control=lmtp_control(.trim = 1)
)

print(result_simple)

result_static<- lmtp_tmle(
  data = data,
  trt = A,
  outcome = Y,
  time_vary = L,
  cens = NULL,
  shift = static_shift,
  mtp = TRUE,
  outcome_type = "continuous",  # Changed from "survival" to "continuous"
  folds = 5,
  learners_trt = c("SL.glm"),
  learners_outcome = c("SL.glm"),
  control=lmtp_control(.trim = 1)
)

print(result_static)

#Static Truth
data1 <- genData(n=100000, p=4, intervene_A=1, intervene_C=1)
mean(data1$Y5)
#############################################################################
# SUMMARY OF MODIFICATIONS:
#
# 1. CHANGED TO CONTINUOUS OUTCOME:
#    - Removed Y2, Y3, Y4 (survival time points)
#    - Y5 is now a continuous outcome with treatment effects
#    - Changed outcome_type to "continuous"
#
# 2. CONTROLLED PROPENSITY SCORE DISTRIBUTIONS:
#    - Modified L1 generation to create bimodal distribution
#    - Adjusted treatment probabilities to achieve:
#      * ~50% with g1 < 0.05 (low baseline propensity)
#      * All with g2, g3, g4 > 0.3 (higher later propensities)
#    - Set alpha = 0.05
#
# 3. ENHANCED DATA GENERATION:
#    - More explicit control over treatment probabilities
#    - Stronger effects for later treatments
#    - Continuous outcome with cumulative treatment effects
#
# The result estimates the effect of Mark & Maya's dynamic intervention
# policy on the continuous outcome Y5.
#############################################################################

