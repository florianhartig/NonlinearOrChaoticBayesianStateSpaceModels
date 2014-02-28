################################################################################
# Code for reply to Perretti et al
#
# Author: Florian Hartig, http://florianhartig.wordpress.com
#
# This script creates the test data for fit and validation
#
################################################################################

rm(list=ls(all=TRUE))

#setwd("D:/home/Projekte/Papers_in_Progress/13-ResponseToPerrettiEtAl/temp")
setwd("C:/temp")

set.seed(1)


############################################################################
# Creation of chaotic population data with known model parameters

trueR <- 3.7
trueRArray <- seq(from=3.2, to=3.9, length.out=50) # for reproducing Fig S.1
trueK <- 1
trueN0 <- 0.5
sampleSize <- 50
trueprocerror <- 0.005
trueObsErrSD <- sqrt(1/25.5) # approx. 0.1980295 



logistic <- function(K,r,N0,samplesize){
  pop = rep(NA,(samplesize))
  pop[1] = N0
  for (i in 2:(samplesize)){
    pop[i] = pop[i-1] * r * (1 - pop[i-1]/K) * rlnorm(1,0,trueprocerror)
  } 
  return(pmax(0,pop))
}


truepops <- matrix(nrow = 50, ncol = 2*sampleSize)
observedpops <- matrix(nrow = 50, ncol = 2*sampleSize)
truepops2 <- matrix(nrow = 50, ncol = 2*sampleSize)
observedpops2 <- matrix(nrow = 50, ncol = 2*sampleSize)


# 50 population dynamics with fixed r=3.7, to reproduce results of Figs.S5,S6, 
# as well as all results regarding predictive accuracy of the logistic model such 
# as Fig 1A  
for (i in 1:50){
  pop <- logistic(trueK,trueR,trueN0,2*sampleSize)
  truepops[i,] <- pop 
  observedpops[i,] =rlnorm(length(pop), mean = log(pop) , sd = trueObsErrSD)
}

# 50 population dynamics with varying R, to reproduce Fig S1, S2
for (i in 1:50){
  pop <- logistic(trueK,trueRArray[i],trueN0,2*sampleSize)
  truepops2[i,] <- pop 
  observedpops2[i,] =rlnorm(length(pop), mean = log(pop) , sd = trueObsErrSD)
}


#Testing
pdf("populationdynamics.pdf")

pop <- truepops[1,]
observed <- observedpops[1,]
plot(observed, xlab = "Time interval", ylab = "Population size / K", main = "", ylim=c(0,2), cex.lab = 1.2, cex.axis = 1.2)
lines(pop)
dev.off()



# save data for systematic analysis 
save.image(file = "preparation.RData")

