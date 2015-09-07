source("Packages_Functions_DataCreation.R")

##--------------------------------
## Blowfly case study
##--------------------------------
blowfly <- read.csv("blowfly97I.txt")
plot(blowfly$total, type="b")

plot(head(blowfly$total,-1), tail(blowfly$total,-1))
lines(smooth.spline(head(blowfly$total,-1), tail(blowfly$total,-1)))




# Synthetic Likelihood
# Mostly adapted from Matheo2014


blowflySimulator <- function(param , nsim, extraArgs, ...){matrix(rlnorm(extraArgs$nObs*nsim,log(blowSimul(param =param[-7], extraArgs=extraArgs, nsim=nsim)),exp(param[7])),nrow=nsim)}

blow_sl <- synlik(simulator = blowflySimulator,
                  summaries = blowStats,
                  param = log( c( "delta" = 0.16, "P" = 6.5, "N0" = 400, 
                                  "var.p" = 0.1, "tau" = 14, "var.d" = 0.1, "obsSD"=0.15)  ),
                  extraArgs = list("nObs" = length(blowfly$total), "nBurn" = 200, "steps" = 2),
                  plotFun = function(input, ...){ 
                    plot(drop(input), type = 'l', ylab = "Pop", xlab = "Time", ...)
                  }
)

blow_sl@data <- rlnorm(360, log(blowfly$total[-361]), 0.2)
blow_sl@extraArgs$obsData <- blow_sl@data

checkNorm(blow_sl)


blow_smcmc <- smcmc(blow_sl, 
                    initPar = log( c( "delta" = 0.1, "P" = 8, "N0" = 600,  "sig.p" = 0.2, "tau" = 17, "sig.d" = 0.3, "obsSD"=0.1)  ),
                    niter = 50000, 
                    burn = 0,
                    propCov = diag(rep(0.001, 7)),
                    nsim = 250, 
                    prior = function(input, ...){
                      sum(input) +
                        dunif(input[4], log(0.01), log(1), log = TRUE) +
                        dunif(input[6], log(0.01), log(1), log = TRUE)
                    },
                    targetRate = 0.15,
                    multicore = T
)

save(blow_smcmc, file="2015-08-12_blow_smcmc_obs_err")

data(blow_smcmc)
tmpTrans <- rep("exp", 6)
names(tmpTrans) <- names(blow_smcmc@param)
plot(blow_smcmc, trans = tmpTrans)

plot(blowfly$total, type="l"); lines(blowSimul(param=apply(blow_smcmc@chains, 2, median), extraArgs=list("nObs"=length(blowfly$total), "nBurn"=200, "steps"=2), nsim=1),type="l", col=2)

plot(seq(from = 1, 2*length(blowfly$total), by=2), blowfly$total, type="b")
lines(bf1)


#### Segmentation Method


# Fit Blowfly Model with jags

obs <-  rlnorm(360, log(blowfly$total[-361]), 0.2)#blowfly$total[-361]
nObs <- length(obs) 
segmentSize <- 30
inits <- list(list(sigma_p = 0.12, sigma_d = 0.07, delta = 0.15,  p0 = 400, P = 5, obsSD= 0.1),
              list(sigma_p = 0.3, sigma_d = 0.13, delta = 0.2,  p0 = 800, P = 7, obsSD= 0.3))
data = list(observed = matrix(obs, nrow=segmentSize), steps=nObs/segmentSize, substeps = segmentSize, tau=14)

cl <- registerCluster(n = 15)
jagslist <- foreach(tau = 1:15, .packages="R2jags") %dopar%{
  
  data = list(observed = matrix(obs, nrow=segmentSize), steps=nObs/segmentSize, substeps = segmentSize, tau=tau)
  r2jags_samples <- jags(data= data, inits = inits, model.file = "model2_blowfly.jags", parameters.to.save = c("sigma_p", "sigma_d", "delta", "p0", "P" , "obsSD"), n.chains = 2, n.burnin = 500000, n.iter = 1000000, n.thin = 30, DIC = TRUE)
}

stopCluster(cl)

tau <- which.min(sapply(jagslist, function(x) x$BUGSoutput$DIC))*2
r2jags_samples <- jagslist[[tau/2]]
plot(r2jags_samples)
r2jags_samples
traceplot(r2jags_samples)

plot(obs, type="l");lines(blowfunc(
  P= r2jags_samples$BUGSoutput$summary[1,1], 
  delta=r2jags_samples$BUGSoutput$summary[2,1],
  p0=r2jags_samples$BUGSoutput$summary[5,1], 
  sigma_d = r2jags_samples$BUGSoutput$summary[6,1], 
  sigma_p = r2jags_samples$BUGSoutput$summary[7,1],
  tau = tau,
  nObs=360, step=2, burn=0), col=2) ;lines(blowfly$total[-361], col=3)


# Plot

ylim=c(0,15000)
par(mfrow=c(3,1))
plot(seq(1,by=2, length.out=length(blowfly$total)), blowfly$total, main="Nicholsons original Dataset", xlab="day", ylab=expression("p"[t]), type="l", ylim=ylim)

plot(seq(1,by=2, length.out=length(blowfly$total)), blowSimul(param=apply(blow_smcmc@chains, 2, median), extraArgs=list("nObs"=length(blowfly$total), "nBurn"=200, "steps"=2), nsim=1), main="Synthetic Likelihood Model", xlab="day", ylab=expression("p"[t]), type="l", ylim=ylim)

plot(seq(1,by=2, length.out=length(blowfly$total)),blowfunc(
  P= r2jags_samples$BUGSoutput$summary[1,1], 
  delta=r2jags_samples$BUGSoutput$summary[2,1],
  p0=r2jags_samples$BUGSoutput$summary[5,1], 
  sigma_d = r2jags_samples$BUGSoutput$summary[6,1], 
  sigma_p = r2jags_samples$BUGSoutput$summary[7,1],
  tau = tau,
  nObs=360, step=2, burn=0), main="Segmentation state-space Model", xlab="day", ylab=expression("p"[t]), type="l", ylim=ylim)

par(mfrow=c(1,1))


# GAM 

library(mgcv)
tau <- 1
blowfly2 <- data.frame(total=rlnorm(360, log(blowfly$total[-361]), 0.2))
blowfly2$total2 <- c(tail(blowfly2$total,-tau), rep(NA,tau))

blowfly_gam <- gam(total ~ s(total2), data= blowfly2[,])
blowfly_gam_predict <- predict(blowfly_gam, newdata = data.frame(total2=blowfly$total))
plot(blowfly$total, type="l")
lines(1:length(blowfly$total), blowfly_gam_predict, col="green")

blowfly_gam_response <- predict(blowfly_gam, newdata=data.frame(total2=seq(1,15000,10)))
plot(seq(1,15000,10),blowfly_gam_response, type="l")

p <- runif(1,1,15000)
for(i in 1:length(blowfly2$total)){
  p <- c(p, predict(blowfly_gam, newdata = data.frame(total2=tail(p,1))))
}
lines(p)

