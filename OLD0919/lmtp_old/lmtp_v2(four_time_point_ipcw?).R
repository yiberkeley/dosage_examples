library(lmtp)
library(data.table)

# Enhanced Data generation function with 4 time points and censoring
genData <- function(n, p, intervene_A = NULL, intervene_C = NULL) {
  
  # Exogenous variables for 4 time points
  for(i in 1:4){
    assign(paste0("U.Lt", i), rnorm(n, mean=0, sd=1))
    assign(paste0("U.At", i), rnorm(n, 0, 1))
    assign(paste0("U.Ct", i), rnorm(n, 0, 1))
    if(i <= 4) {
      assign(paste0("U.Yt", i+1), rnorm(n, 0, 1))
    }
  }
  
  # Initialize lists to store variables
  L <- vector("list", 4)
  A <- vector("list", 4)
  C <- vector("list", 4)
  Y <- vector("list", 5)
  
  # Time point 1
  L[[1]] = as.numeric(U.Lt1 < 0.4) # Most people controlled at baseline
  
  # Treatment at time 1
  if(is.null(intervene_A)){
    A[[1]] = as.numeric(U.At1 + 1 < plogis(L[[1]]))
  } else {
    A[[1]] <- rep(intervene_A, n)
  }
  
  # Censoring at time 1 (1 = not censored, 0 = censored)
  if(is.null(intervene_C)){
    C[[1]] = as.numeric(U.Ct1 < plogis(0.5 + 0.3*L[[1]] - 0.2*A[[1]]))
  } else {
    C[[1]] <- rep(intervene_C, n)
  }
  
  # Early outcome at time 2 (only for non-censored)
  Y[[2]] = ifelse(C[[1]] == 1, 
                  as.numeric(U.Yt2 > plogis(L[[1]] + (A[[1]]*2))), 
                  NA)
  
  # Time point 2 (only for those not censored and no early outcome)
  continue_t2 <- C[[1]] == 1 & (is.na(Y[[2]]) | Y[[2]] == 0)
  
  L[[2]] = ifelse(continue_t2, 
                  as.numeric(plogis(A[[1]] + L[[1]] + U.Lt2) > 0.4), 
                  NA)
  
  if(is.null(intervene_A)){
    A[[2]] = ifelse(continue_t2,
                    as.numeric(U.At2 + 1 < plogis(L[[2]] + A[[1]] + L[[1]])),
                    NA)
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
  
  Y[[3]] = ifelse(continue_t2 & C[[2]] == 1,
                  as.numeric(U.Yt3 + 1 < plogis(L[[1]] + A[[1]] + L[[2]] + A[[2]])),
                  ifelse(continue_t2 & C[[2]] == 0, NA, Y[[2]]))
  
  # Time point 3
  continue_t3 <- continue_t2 & C[[2]] == 1 & (is.na(Y[[3]]) | Y[[3]] == 0)
  
  L[[3]] = ifelse(continue_t3,
                  as.numeric(plogis(A[[2]] + L[[2]] + U.Lt3) > 0.4),
                  NA)
  
  if(is.null(intervene_A)){
    A[[3]] = ifelse(continue_t3,
                    as.numeric(U.At3 + 1 < plogis(L[[3]] + A[[2]] + L[[2]])),
                    NA)
  } else {
    A[[3]] <- ifelse(continue_t3, intervene_A, NA)
  }
  
  if(is.null(intervene_C)){
    C[[3]] = ifelse(continue_t3,
                    as.numeric(U.Ct3 < plogis(0.5 + 0.3*L[[3]] - 0.2*A[[3]] + 0.1*A[[2]])),
                    NA)
  } else {
    C[[3]] <- ifelse(continue_t3, intervene_C, NA)
  }
  
  Y[[4]] = ifelse(continue_t3 & C[[3]] == 1,
                  as.numeric(U.Yt4 + 1 < plogis(L[[2]] + A[[2]] + L[[3]] + A[[3]])),
                  ifelse(continue_t3 & C[[3]] == 0, NA, Y[[3]]))
  
  # Time point 4
  continue_t4 <- continue_t3 & C[[3]] == 1 & (is.na(Y[[4]]) | Y[[4]] == 0)
  
  L[[4]] = ifelse(continue_t4,
                  as.numeric(plogis(A[[3]] + L[[3]] + U.Lt4) > 0.4),
                  NA)
  
  if(is.null(intervene_A)){
    A[[4]] = ifelse(continue_t4,
                    as.numeric(U.At4 + 1 < plogis(L[[4]] + A[[3]] + L[[3]])),
                    NA)
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
  
  # Final outcome at time 5
  Y[[5]] = ifelse(continue_t4 & C[[4]] == 1,
                  as.numeric(U.Yt5 + 1 < plogis(L[[3]] + A[[3]] + L[[4]] + A[[4]])),
                  ifelse(continue_t4 & C[[4]] == 0, NA, Y[[4]]))
  
  # Create dataframe with proper naming convention
  ObsData <- data.table(
    L1 = L[[1]], A1 = A[[1]], C1 = C[[1]],
    L2 = L[[2]], A2 = A[[2]], C2 = C[[2]], 
    L3 = L[[3]], A3 = A[[3]], C3 = C[[3]],
    L4 = L[[4]], A4 = A[[4]], C4 = C[[4]],
    Y2 = Y[[2]], Y3 = Y[[3]], Y4 = Y[[4]], Y5 = Y[[5]]
  )
  
  # For survival analysis, carry forward last observation
  for(i in 3:5) {
    y_col <- paste0("Y", i)
    prev_y_col <- paste0("Y", i-1)
    if(i > 3) {
      ObsData[[y_col]] <- ifelse(is.na(ObsData[[y_col]]) & !is.na(ObsData[[prev_y_col]]),
                                 ObsData[[prev_y_col]], ObsData[[y_col]])
    }
  }
  
  return(ObsData)
}

# Generate data
set.seed(123)
data <- genData(n=1000, p=4, intervene_A=NULL, intervene_C=NULL)
alpha <- 0.1

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

# Run analysis with 4 time points, censoring, and Mark & Maya's intervention
data <- as.data.frame(data)

# Define variables for lmtp
A <- c("A1", "A2", "A3", "A4")  # Treatment variables
L <- list(c("L1"), c("L2"), c("L3"), c("L4"))  # Time-varying covariates
C <- c("C1", "C2", "C3", "C4")  # Censoring indicators
Y <- c("Y2", "Y3", "Y4", "Y5")  # Outcome variables

cat("Running LMTP analysis with Mark & Maya's cumulative propensity intervention...\n")
cat("Data structure:\n")
cat("- Treatment nodes:", paste(A, collapse=", "), "\n")
cat("- Censoring nodes:", paste(C, collapse=", "), "\n") 
cat("- Outcome nodes:", paste(Y, collapse=", "), "\n")
cat("- Time-varying covariates: L1, L2, L3, L4\n")
cat("- Alpha threshold:", alpha, "\n\n")

# Check data structure
cat("Sample of generated data:\n")
print(head(data[, c("L1", "A1", "C1", "L2", "A2", "C2", "L3", "A3", "C3", "L4", "A4", "C4", "Y5")]))

# Run the analysis
result_mark_maya <- lmtp_tmle(
  data = data,
  trt = A,
  outcome = Y,
  time_vary = L,
  cens = C,
  shift = mark_maya_shift,
  mtp = TRUE,
  outcome_type = "survival",
  folds = 5,
  learners_trt = c("SL.glm"),
  learners_outcome = c("SL.glm")
)

print(result_mark_maya)

#############################################################################
# SUMMARY OF MODIFICATIONS:
#
# 1. EXTENDED TO 4 TIME POINTS:
#    - Structure: L1->A1->C1->L2->A2->C2->L3->A3->C3->L4->A4->C4->Y5
#    - Each time point has proper dependency structure
#
# 2. ADDED CENSORING NODES (C1, C2, C3, C4):
#    - C_t indicates if observation continues after time t
#    - Affects subsequent L, A, and Y variables
#    - Intervention sets all C nodes to 1 (no censoring)
#
# 3. IMPLEMENTED MARK & MAYA'S PROPOSAL:
#    - Tracks cumulative product of propensity scores
#    - Only intervenes when cumulative propensity > alpha
#    - Skips intervention when "very unlikely" (≤ alpha)
#    - Handles missing data appropriately
#
# 4. PROPER SURVIVAL OUTCOME STRUCTURE:
#    - Multiple outcome time points (Y2, Y3, Y4, Y5)
#    - Last observation carried forward for survival analysis
#    - Handles censoring in outcome generation
#
# The result estimates the effect of Mark & Maya's dynamic intervention
# policy combined with a no-censoring intervention.
#############################################################################