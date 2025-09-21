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



##################################################
####### WORKING MSM EXAMPLE: time on max dose 
##################################################
time.points = 2 #number of anodes
regime.matrix <- as.matrix(expand.grid(rep(list(0:1), time.points)))# matrix of all possible regimes given number of anodes 
num.regimes <- nrow(regime.matrix)
regimes <- array(dim = c(nrow(data), time.points, num.regimes)) # nrow(data) x time.points x num.regimes (array to summarize over)
summary.measures <- array(dim = c(num.regimes, 1, 1)) #how we will summarize the the regimes from regime-matrix
for (i in 1:num.regimes) {
  regimes[, , i] <- matrix(regime.matrix[i, ], byrow = TRUE, nrow = nrow(data), ncol = time.points)
  summary.measures[i, 1, 1] <- sum(regime.matrix[i, ])
}
colnames(summary.measures) <- "time.on.treatment"
#remove other Y node
data <- data %>% select(-Y2)
resultMSM <- ltmleMSM(data, 
                      Anodes=c('A1','A2'), 
                      Lnodes=c('L1','L2'), 
                      Ynodes='Y3',
                    survivalOutcome=FALSE,
                    regimes=regimes,
                    summary.measures=summary.measures, 
                    working.msm="Y ~ time.on.treatment")

(summary(resultMSM))



##################################################
####### stochatic/dynamic treatment based on propensity score
##################################################
#remove other Y node
data <- data %>% select(-Y2)

res1 <- ltmle(data,
                  Anodes=c("A1","A2"),
                  Lnodes=c('L1','L2'),
                  Cnodes=NULL,
                  Ynodes='Y3',
                  abar= list(c(1,1),c(0,0)),
                  survivalOutcome=F)

cumgA1 <- as.matrix(res1$cum.g.unbounded[, , 1])
table(cumgA1<0.05)
##these are all common, so overwriting a few for the sake of the example
cumgA1[2,1] <- 0.01
cumgA1[3,1] <- 0.01
table(cumgA1<0.05)

#create abar matrices
abarA1 <- matrix(rep(1,length=nrow(cumgA1)*ncol(cumgA1)),nrow=nrow(cumgA1),ncol=ncol(cumgA1))
abarA0 <- matrix(rep(0,length=nrow(cumgA1)*ncol(cumgA1)),nrow=nrow(cumgA1),ncol=ncol(cumgA1))

#overwrite abar matrix with natural value in cases where cumA1<0.05
observedA <- as.matrix(data[,c('A1','A2')])
table(observedA<0.05)
abarA1[cumgA1<0.05] <- observedA[cumgA1<0.05]
abarA0[cumgA1<0.05] <- observedA[cumgA1<0.05]

#check
observedA[2,]
observedA[2,]
cumgA1[2,]
cumgA1[3,]
abarA1[2,]
abarA1[3,]

##now run analysis with updated values
res2 <- ltmle(data,
              Anodes=c("A1","A2"),
              Lnodes=c('L1','L2'),
              Cnodes=NULL,
              Ynodes='Y3',
              abar= list(abarA1,abarA0),
              survivalOutcome=F)

summary(res2)
#compare to pre:
#summary(res1)

##################################################
####### WORKING MSM EXAMPLE: dose threshold theta
##################################################
# here we wish to define a marginal structural model of threshold dose level, theta. The regimes will be defined as 
#d_θ: Set the dosage level to at least θ at the start of follow-up (after the given lag period), 
#and continue this dosage throughout follow-up, or until the occurrence of a loss-to-management event,
#at which point allow the physician to intervene as observed.

