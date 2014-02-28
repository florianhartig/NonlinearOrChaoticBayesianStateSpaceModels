################################################################################
# Code for reply to Perretti et al
#
# Author: Florian Hartig, http://florianhartig.wordpress.com
#
# This scripts runs the fits for the 50 test datasets and evaluates
# the posterior mode and and some other summary statistics 
#
# Requires data creation and modelSpecification to be run before
#
################################################################################

rm(list=ls(all=TRUE))

library(rjags) 
library(modeest)

setwd("C:/temp")
load(file = "preparation.RData")


# lumpchains combines the chains in mcmc.list
# copied from emdbook package by Ben Bolker

lumpchains <- function(x){
  x2 <- do.call("rbind", x)
  mcpars <- sapply(x, attr, "mcpar")
  class(x2) <- "mcmc"
  if (var(mcpars[1, ]) > 0 || var(mcpars[3, ]) > 0) 
    stop("can't combine chains with unequal start/thin")
  attr(x2, "mcpar") <- c(mcpars[1, 1], sum(mcpars[2, ]), 
                         mcpars[3,1])
  x2
}


# inits - for optional use later

inits <- list(
  list("r" = 2.2, "K"=trueK, "N0" = 0.5, "ObservationPrecision"=4),
  list("r" = 2.5, "K"=trueK, "N0" = 0.5, "ObservationPrecision"=4),
  list("r" = 2.8, "K"=trueK, "N0" = 0.5, "ObservationPrecision"=4)
)


# Functions to create estimates from multiple samples

getestimates <- function(observed, model, params = c("r", "K", "obsSD", "pop[1]"), 
    alternative = 1, adapt = 5000, update = 5000, mainiter = 100000, thin = 100, chains = 3, useinits = F){
  

  
  lenght <- dim(observed)[1]
  datapoints <- dim(observed)[2]
  posteriormodes <- matrix(nrow = lenght, ncol = 4)
  posteriormedians <- matrix(nrow = lenght, ncol = 4)
  posteriorrandom <- matrix(nrow = lenght, ncol = 4)
  psrf <- rep(NA,lenght)
  sel <- sample(1:(mainiter / thin * chains),100,replace=T) # selecting a random parameter from the posterior, for the posterior predictive 
  
  for (i in 1:lenght){
    ptm <- proc.time()
    if(alternative == 1){
      data = list(observed = observed[i,], nObs=datapoints)      
    }
    if(alternative == 2){
      data = list(observed = matrix(observed[i,], nrow=5), steps=10, substeps = 5)      
    }
    if(alternative == 3){
      data = list(observed = log(observed[i,]), nObs=datapoints)      
    }
    if (useinits == T) jagsModel = jags.model( file= model, data=data, n.chains = chains, n.adapt= adapt, inits = inits)
    else jagsModel = jags.model( file= model, data=data, n.chains = chains, n.adapt= adapt)
    update(jagsModel, update)
    results = coda.samples( jagsModel , variable.names=params,n.iter=mainiter, thin = thin)
    png(paste("l", run,"results.",i,".png", sep=""), width= 1000, height = 1000)
    plot(results)
    dev.off()
    png(paste("l", run,"results",i,"gp.png", sep=""), width= 1000, height = 1000)
    gelman.plot(results)
    mtext(paste("m-psrf",round(gelman.diag(results)$mpsrf, digits=3)), side = 3, outer = F, padj=-3)
    dev.off()
    lumped <- lumpchains(results)
    posteriormodes[i,] <- apply(lumped,2,venter)
    posteriormedians[i,] <- apply(lumped,2,median)
    posteriorrandom[i,] <- lumped[sel[i], ]
    psrf[i]=round(gelman.diag(results)$mpsrf, digits=3)
    proc.time() - ptm
  }
  return(list(modes = posteriormodes, medians = posteriormedians, random = posteriorrandom, mpsrf = psrf))
}



##################################################
# Fitting the models

run = 122 # to chose which runs to evaluate


if (run == 111 | run == 112) model2resultsfixed <- getestimates(observedpops[,1:50], "model1.txt", 
    adapt = 500000, update = 500000, mainiter = 5000000, thin = 500, useinits = T)
  
if (run == 121 | run == 122) model2resultsvary <- getestimates(observedpops2[,1:50], "model1.txt",
    adapt = 500000, update = 500000, mainiter = 5000000, thin = 500, useinits = T)


if (run == 21){
  model3resultsfixed <- getestimates(observedpops[,1:50], "model2.txt", params = c("r", "K", "obsSD", "pop[1,1]"), alternative = 2,
    adapt = 20000, update = 100000, mainiter = 300000, thin = 30, useinits = T)}

if (run == 22){
  model3resultsvary <- getestimates(observedpops2[,1:50], "model2.txt", params = c("r", "K", "obsSD", "pop[1,1]"), alternative = 2,
    adapt = 20000, update = 100000, mainiter = 300000, thin = 30, useinits = T)
}



save.image(file = paste("results",run,".RData",sep=""))



