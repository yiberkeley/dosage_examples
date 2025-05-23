source("./functions.R")
library(data.table)
library(ltmle)

data <- genData(1000,0.5)
str(data)
##############################
##### stochastic intervention: set to 1 unless adverse event occurs
##### based on a time-varying covariate value
##############################

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

SI_abar1 <- gen_abar_SI(data, Anodes=c('A1','A2'), Lnodes= c('L1','L2'))
head(SI_abar1)
head(data)

res <- run_analysis(data,Lnodes=c("L1","L2"),Cnodes=NULL,
             str_name="stochastic_abar",abar=list(SI_abar1,matrix(0,nrow=nrow(data),ncol=2)))
#need to add in censoring?
