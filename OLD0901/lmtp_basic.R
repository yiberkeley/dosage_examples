library(lmtp)
library(data.table)
library(dplyr)

# Generate data
set.seed(123)
data <- genData(n=1000, p=2, intervene=NULL)

#lmtp needs the data frame
data<-as.data.frame(data)

# Create a vector of the names of the treatment variables ordered by time
A <- c("A1", "A2")



# Create a list of vectors with the covariates at each time point
# with the list ordered by time
L <- list(c("L1"), c("L2"))

# # Create a vector of the names of the censoring variables ordered by time
# C <- c("C1", "C2")

# Create a vector of the names of the outcome variables ordered by time
Y <- c("Y2", "Y3")

# Create a shift function, called f
# Shift functions have two arguments: data and trt
f_1 <- function(data, trt) {
  1
}

f_0 <- function(data, trt) {
  0
}

# Create a vector indicating what algorithms should be
# used in the SuperLearner
sl_lib <- c("SL.glm", "SL.earth")

# Run lmtp_tmle() using  the data, the treatment vector,
# the time-varying covariate list, and the shift function f
y_1<-lmtp_tmle(data, 
          trt=A, 
          outcome=Y, 
          baseline = NULL,
          time_vary = L, 
          cens = NULL,
          shift = f_1,
          mtp = TRUE,
          outcome_type ="survival",
          folds = 5,
          learners_trt = sl_lib, learners_outcome = sl_lib
         )


y_0<-lmtp_tmle(data, 
               trt=A, 
               outcome=Y, 
               baseline = NULL,
               time_vary = L, 
               cens = NULL,
               shift = f_0,
               mtp = TRUE,
               outcome_type ="survival",
               folds = 5,
               learners_trt = sl_lib, learners_outcome = sl_lib)

# Estimate the causal risk ratio/additive difference
lmtp_contrast(y_1, ref = y_0, type = "rr")
lmtp_contrast(y_1, ref = y_0, type = "additive")

##################################################
####### Dynamic treatment
##################################################


f_d <- function(data, trt) {
  ifelse(data[[sub("A", "L", trt)]] == 0,
         1,
         data[[trt]])
}

y_d<-lmtp_tmle(data, 
               trt=A, 
               outcome=Y, 
               baseline = NULL,
               time_vary = L, 
               cens = NULL,
               shift = f_d,
               mtp = TRUE,
               outcome_type ="survival",
               folds = 1,
               learners_trt = sl_lib, learners_outcome = sl_lib)

#NOTE:ONE CAN ALSO MODIFY THE TREATMENT COLUMN DIRECTLY IN THE DATA AND USE THE SHIFTED ARGUMENT
gen_abar_SI <- function(data, Anodes, Lnodes){
  
  #need to create a matrix the length of the dataset, one obs per patient
  # and the width of the Anodes
  abar_mat <- matrix(NA,nrow=nrow(data),ncol=length(Anodes))
  
  #for Anodes, STOCHASTIC INTERVENTION:
  #set to 1, UNLESS they get an adverse event (L nodes), then set to natural value 
  for (j in 1:length(Anodes)) {
    for (i in 1:nrow(data)) {
      if (data[i, get(Lnodes[j])] == 0) {
        # if prior L node for person i is 0 (i.e., no adverse event)
        #set person i's abar to 1 deterministically
        abar_mat[i, j] <- 1
        
      } else if (data[i, get(Lnodes[j])] == 1) {
        #if prior L node for person i is 1 (i.e., adverse event), then allow to vary naturally 
        abar_mat[i, j] <- data[i, get(Anodes[j])]
        
      }
    }
  }
  return(abar_mat)
}

#Change to data table for the counterfactual specification
data<-as.data.table(data)
SI_abar1 <- gen_abar_SI(data, Anodes=c('A1','A2'), Lnodes= c('L1','L2'))

shifted_data<-data
shifted_data$A1<-SI_abar1[,1]
shifted_data$A2<-SI_abar1[,2]

#Change it to the data frame
data<-as.data.frame(data)
shifted_data<-as.data.frame(shifted_data)

y_s<-lmtp_tmle(data, 
          trt=A, 
          outcome=Y, 
          baseline = NULL,
          time_vary = L, 
          cens = NULL,
          shifted = shifted_data,
          mtp = TRUE,
          outcome_type ="survival",
          folds = 1,
          learners_trt = sl_lib, learners_outcome = sl_lib)

#note that gen_abar_SI <- function(data, Anodes, Lnodes) can be more complex

#Stochastic: https://arxiv.org/pdf/2202.03513 under review

#Todo:
#1.Specifying control in the lmtp for the bound using the help(lmtp)
#2.Extract the time varying information like the propensity estimated,