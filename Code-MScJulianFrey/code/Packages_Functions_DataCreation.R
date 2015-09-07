setwd("Z:\\Masterarbeit\\GitHub\\code")

#--------------------------------------------------------
# Packages
#--------------------------------------------------------


# install.packages("rstan", repos = "http://cran.rstudio.com",  dependencies = TRUE)

#install.packages(c("modeest","IDPmisc","compiler","devtools","ks","hydroGOF","compiler","GenSA","foreach","doParallel","nloptr", "synlik", "EasyABC"),dep=T)

library(rstan) 
rstan_options(auto_write = TRUE)
options(mc.cores = 1)#parallel::detectCores())

#library(modeest)
library(IDPmisc)
library(compiler)
library(devtools)
#library(ks)
library(hydroGOF)
library(compiler)
library(GenSA)
library("foreach")
library("doParallel")
#library(nloptr)
#library(zoo)
#library(fNonlinear)
library(synlik)
library(rjags)
library(R2jags)
library(coda)
library(mgcv)
#library(EasyABC)
library("nwfscNLTS")


#--------------------------------------------------------
# Functions
#--------------------------------------------------------


# function to extract median-values of parameters from STAN
median_stan <- function(model, parameters){
  ex_m <- extract(model, permuted=T, pars=parameters)
  results <- list()
  for(i in parameters){
    results[i] <- median(ex_m[i][[1]])
  }
  return(results)
}

# function to calculate the states from a stanfit with segmentation method
states <- function(stanFit=stanlist[[1]], FUN="logistic"){
  pop <- as.numeric(median_stan(stanFit, paste0("pop[",1:stanFit@par_dims$pop,"]")))
  pars <- as.numeric(median_stan(stanFit, stanFit@model_pars))
  procERR <- matrix(as.numeric(median_stan(stanlist[[1]], paste0("procerr[",rep(1:stanFit@par_dims$procerr[1],each=stanFit@par_dims$procerr[2]),",",seq(1,stanFit@par_dims$procerr[2]),"]"))), nrow=segmentSize-1, byrow=T)
  states <- numeric()
  if(FUN=="logistic") f <- expression(max(0, p  *  pars[1]  * (1-(p/pars[2])) * procERR[j,i]))
  
  for(i in 1:length(pop)) {
    p <- pop[i]
    
    for(j in 1:(stanFit@par_dims$procerr[1])) {
      states <- c(states,p)
      p <- eval(f)
    } }
  return(states)
}

# Better Pairs from Florian Hartig
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="blue4", ...)
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, method = "spearman"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex * r)
}

betterPairs <- function(YourData){
  return(pairs(YourData, lower.panel=function(...) {par(new=TRUE);ipanel.smooth(...)}, diag.panel=panel.hist, upper.panel=panel.cor))
}

# lyapunov exponent 

lyapunov <- function(r, K = 1, fn="logistic"){
  
  burn      <- 100000
  T         <- 1000
  x.burn    <- vector(mode="numeric", length=burn)
  x.burn[1] <- 0.5
  x.a       <- vector(mode="numeric", length=T)
  x.b       <- vector(mode="numeric", length=T)
  lyap      <- vector(mode="numeric", length=T-1)
  avlyap    <- vector(mode="numeric", length=T-1)
  tau.proc  <- 40000
  sd.proc   <- sqrt(1/tau.proc)
  d0        <- 10^-12 # Distance to separate orbits
  
  if(fn == "logistic"){fun <- function(r, K, p){max(0,r*p*(1-p/K))}}
  else if(fn == "ricker"){fun <- function(r, K, p){max(0,exp(r)*p*exp(-p/K))}}
  else{fun <- fn}
  
  
  ## Process model
  # Burn-in
  for(t in 1:(burn-1)){
    x.burn[t+1] <- fun(r,K,p=x.burn[t])
  }
  
  # Calculate lyapunov exponent
  x.a[1] <- x.burn[length(x.burn)] # Initiate on attractor
  x.b[1] <- x.burn[length(x.burn)] + d0
  
  for(t in 1:(T-1)){
    
    x.a[t+1] <- fun(r,K,p=x.a[t]) #max(0,exp(r)*x.a[t]*(exp(-x.a[t]/K)))
    x.b[t+1] <- fun(r,K,p=x.b[t]) #max(0,exp(r)*x.b[t]*(exp(-x.b[t]/K)))
    
    d1        <- sqrt((x.a[t+1] - x.b[t+1])^2) # Calculate distance between orbits
    lyap[t]   <- log(abs(d1/d0))
    avlyap[t] <- mean(lyap[1:t])
    
    x.b[t+1] <- x.a[t+1] + d0*(x.b[t+1]-x.a[t+1])/d1
  }
  
  #plot(as.ts(x.b),plot.type='single')
  #plot(as.ts(lyap))
  #plot(as.ts(avlyap))
  
  return(mean(lyap[(T-1-300):(T-1)]))
  
}

lyapunov <- cmpfun(lyapunov)

############# Basic function to create population time series
# r = growth rate
# K = carriing capathity
# p0 = initial population counts/size with length=lag
# lag = time lag from one generation to another
# nObs = number of observations (excl. burn-in)
# burn = burn in phase
# rf = realization factor optinal parameter if lag is a vector to take into account that offspring of multiple time steps is fertile at once 
population <- function(r, K, p0=runif(1), procSD=0, burn=0, nObs=50, lag=1, FUN="logistic", rf=1){
  l <- max(lag)
  if(length(rf) != length(lag)) rf <- 1/length(lag)
  if(FUN=="logistic") f <- expression(max(0,
                                          
                                          sum(Obs[(i-lag[length(lag)]):(i-lag[1])]*rf)  *  r  * (1-((sum(Obs[(i-lag[length(lag)]):(i-lag[1])]*rf)/K)))) * rlnorm(1,0,procSD))
  
  else if(FUN=="ricker") f <- expression(max(0,sum(Obs[(i-lag[length(lag)]):(i-lag[1])]*rf)*exp(r)*exp(-((sum(Obs[(i-lag[length(lag)]):(i-lag[1])]*rf)/K)))) * rlnorm(1,0,procSD))
  else f <- FUN
  if(length(p0)<l) p0 <- rep(p0,l)[1:l]
  Obs <- p0
  for(i in (l+1):(burn+nObs+1)){
    Obs <- c(Obs,eval(f))
  }
  return(Obs[(burn+1):(burn+nObs)])
}
population <- cmpfun(population)


### cluster settings

registerCluster <- function(n = "cores"){
  
  library(parallel)
  library(foreach)
  library(doParallel)
  if(n == "cores") n <- detectCores()
  else numCores <- n
  cl <- makeCluster(numCores)
  clusterExport(cl, c(lsf.str(), "population" ))
  registerDoParallel(cl)
  
  return(cl)
  
}


#function to calculate the likelihood of logistic map with split time-series 
likelihood <- function(r, K=1, p0=0.5, TrueObs=trueObs, segmentSize=nObs, sd=obsSD){ 

  if(length(TrueObs!=segmentSize)) {SPLIT <- matrix(TrueObs[1:(length(TrueObs)- (length(TrueObs) %% segmentSize))], nrow=segmentSize)} else SPLIT <- TrueObs

  if(length(p0)<ncol(SPLIT)) {p0 <- c(p0, SPLIT[nrow(SPLIT),length(p0):ncol(SPLIT)])}
  results <- numeric()
  
  for(i in 1:ncol(SPLIT)){
    tO <- SPLIT[,i]
    Obs <- px <- p0[i]
    f <- expression(r*px*(1-px/K))
    for(i in 1:(length(tO)-1)){
      px <- max(0,eval(f))
      Obs <- c(Obs,px)
    }
    results <- sum(c(results, sum(dlnorm(tO, log(Obs), log=T, sd=sd))))
  }
  results
}

likelihood <- cmpfun(likelihood)

# blowflys simulation function

blowfunc <- function(nObs= 320, burn= 200, sigma_p = 0.2, sigma_d = 0.3, delta = 0.1, tau = 17, p0 = 600, P = 8, step = 2){
  
  nObs <- nObs*step
  R_t <- expression( rpois(1,P * p_tau * exp(-p_tau/p0) * (e))) 
  S_t <- expression( rbinom(1,p,exp(-delta*(epsilon)))) 
  
  
  
  l <- nObs
  sigma_p <- 1+sigma_p^2
  sigma_d <- sigma_d^2
  shape <- sigma_p^2/sigma_d
  scale <- sigma_d/sigma_p
  
  
  ev <- rgamma(l+burn+tau,shape=shape,scale=scale)
  epsilonv <- rgamma(l+burn+tau,shape=shape,scale=scale)
  p <- p0 <- 600
  
  pop <- numeric()
  for(i in 1:tau){epsilon <- epsilonv[i]; p <- eval(S_t); pop <- c(pop,p) }
  
  
  for(i in tau:(l+burn)){
    e <- ev[i]
    epsilon <- epsilonv[i]
    p_tau <- pop[i-tau+1]
    pop <- c(pop,p)
    p <- eval(R_t)+eval(S_t)
  }
  pop <- pop[(burn):(l+burn)]
  return(pop[1:(step)==step])
}

# stan models

# logistic map model without segmentation
sink("model1.stan")

cat(
  "data{
	real<lower=0> observed[50];
  int nObs;
  }
  parameters{
  real<lower=0,upper=15> r;
  real<lower=0,upper=10> K;
  real<lower=0, upper=200> obsSD;
  real<lower=0, upper=10> pop;
  vector<lower=0,upper=40>[nObs-1] procerr;
  }
  
  
  transformed parameters{
  real procSD;
  procSD <- 0.005;
  }
  
  model{
  vector[nObs] logpop;
  r ~ uniform(0.0,4.5);
  K ~ uniform(0.0,10);
  obsSD ~ uniform(0,5);
  
  pop ~ uniform(0,1.5);  
  observed[1] ~ lognormal(log(pop), obsSD);
  logpop[1] <- pop;
  
  
  for(i in 1:(nObs-1)){
  procerr[i] ~ lognormal(0,procSD);
  logpop[i+1] <- fmax(0,r * logpop[i] *  (1 - ( logpop[i] / K ) ) * procerr[i]);
  observed[i+1] ~ lognormal(log(fmax(0.000000001,logpop[i+1])), obsSD);
  
  }
  
  }
  "
  ,fill = TRUE)

sink()


# logistic map model with segmentation
sink("model2.stan")

cat(
"data {
	int steps;
int substeps;
matrix<lower=0>[substeps, steps] observed;
}
parameters {
real<lower=0,upper=20> r;
real<lower=0,upper=15> K;	
real<lower=0,upper=200> obsSD;
vector<lower=0,upper=45>[steps] pop;
matrix<lower=0,upper=45>[substeps-1,steps] procerr;
}
transformed parameters{
real procSD;
procSD <- 0.005;
}

model{
matrix[substeps,steps] logpop;
r ~ uniform(0.0,4.5);
K ~ uniform(0.0,10);
obsSD ~ uniform(0,5);


# Likelihood

for(i in 1:steps){

pop[i] ~ uniform(0,1.5);  
observed[1,i] ~ lognormal(log(pop[i]), obsSD);
logpop[1,i] <- pop[i];


for(j in 1:(substeps-1)){
procerr[j,i] ~ lognormal(0,procSD);
logpop[j+1,i] <- fmax(0,r * logpop[j,i] *  (1 - ( logpop[j,i] / K ) ) * procerr[j,i]);
observed[j+1,i] ~ lognormal(log(fmax(0.000000001,logpop[1+j,i])), obsSD);

}
}
}"
,fill = TRUE)

sink()

sink("model2_blowfly.jags")

cat(
  "model{ 
  #priors and transformations
  sigma_p ~ dunif(0,0.4)
  sigma_d ~ dunif(0,0.4)
  delta  ~ dunif(0,0.2)
  
  p0 ~ dunif(300,800)
  P ~ dunif(0.0,10)
  obsSD ~ dunif(0,0.5)
  
  shape <- (1+sigma_p^2)^2/sigma_d^2
  scale <- sigma_d^2/(1+sigma_p^2)
  
  epsilon ~ dgamma(shape, 1/scale)
  e ~ dgamma(shape, 1/scale)
  
  # Likelihood
  
  for(i in 1:steps){
  
  
  for(l in 1:(tau)){
  logpop[l,i] ~ dunif(0,10000)
  logpop1[l,i] ~ dunif(0,10000)
  }
  logpop[tau+1,i] ~ dunif(0,10000)
  for(j in (tau+1):(substeps+tau)){
  
  R_t1[j,i] ~ dpois(P * round(logpop[j-tau,i]) * exp(-round(logpop[j-tau,i])/p0) * e)
  S_t1[j,i] ~ dbin(exp(-delta*(epsilon)),round(logpop[j,i]))
  logpop1[j,i] <- max(0.0001, S_t1[j,i] + R_t1[j,i])
  R_t[j,i] ~ dpois(P * round(logpop1[j-tau,i]) * exp(-round(logpop1[j-tau,i])/p0) * e)
  S_t[j,i] ~ dbin(exp(-delta*(epsilon)),round(logpop1[j,i]))
  logpop[j+1,i] <- max(0.0001, S_t[j,i] + R_t[j,i])
  }
  for(k in 1:substeps){
  observed[k,i] ~ dlnorm(log(logpop[tau+k,i]), pow(obsSD,2))
  }
  }
  }
  "
  
,fill = TRUE)

sink()


#--------------------------------------------------------
# Data Creation
#--------------------------------------------------------
set.seed(1337)
obsSD <- sqrt(1/25.5)
r <- 3.7
K <- 1
nObs = 50

truePop <- population(r=r, K= K, p0 = 0.5, procSD = 0.005, burn = 100, nObs=nObs)
trueObs <- rlnorm(length(truePop), mean = log(truePop) , sd = obsSD)