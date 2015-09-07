################################################################################
# 
#
# Author: Florian Hartig, http://florianhartig.wordpress.com
#
# logistic model
#
################################################################################


logisticModel <- function(K,r,N0,samplesize, processerror){
  pop = rep(NA,(samplesize))
  pop[1] = N0
  for (i in 2:(samplesize)){
    pop[i] = pop[i-1] * r * (1 - pop[i-1]/K) * rlnorm(1,0,processerror)
  } 
  return(pmax(0,pop))
}


logisticModelC <- cmpfun(logisticModel)
