source("Packages_Functions_DataCreation.R")

#--------------------------------------------------------
# Fitting the logistic map STAN model without segmentation
#--------------------------------------------------------
nIter <- 50000
warm <- 25000


data = list(observed = trueObs[1:nObs], nObs=50)

inits1 <- list(
  list("r" = 2.5, chain_id=1, "K"=1, "obsSD"=0.1, "pop"=0.5, "procerr"=rep(1, nObs-1)),
  list("r" = 3.0, chain_id=2, "K"=1, "obsSD"=0.1, "pop"=0.5, "procerr"=rep(1, nObs-1)),
  list("r" = 3.5, chain_id=3, "K"=1, "obsSD"=0.1, "pop"=0.5, "procerr"=rep(1, nObs-1))
)

sM1 <- stan_model(file= "model1.stan")
stanModel1 = sampling( sM1 , data=data, chains = 3, iter=nIter, warmup=warm, init=inits1, thin=30) 

params1 <- c("r", "K", "ObservationPrecision", "pop[1]")
stanModel1

plot(stanModel1)

rstan::traceplot(stanModel1, params1)

# extract median parameter estimates
median_stan(stanModel1, params1)


par(mfrow=(c(2,3)))
matplot(seq(1,by=30, length.out = (((nIter+warm)/30)+1)/3),matrix(extract(stanModel1)$r, ncol=3), type="l", xlab="Iteration", ylab="r", lty=1, col=c("#88888888"))
matplot(seq(1,by=30, length.out = (((nIter+warm)/30)+1)/3),matrix(extract(stanModel1)$K, ncol=3), type="l", xlab="Iteration", ylab="K", lty=1, col=c("#88888888"))
matplot(seq(1,by=30, length.out = (((nIter+warm)/30)+1)/3),matrix(extract(stanModel1)$K, ncol=3), type="l", xlab="Iteration", ylab=expression(sigma^'obs'), lty=1, col=c("#88888888"))

hist(extract(stanModel1)$r, breaks = 10, main="", xlab="r", xlim=c(0,4)) ; abline(v=3.7, lwd=2, lty=2)
hist(extract(stanModel1)$K, breaks = 10, main="", xlab="K", xlim=c(0,4)) ; abline(v=1, lwd=2, lty=2)
hist(extract(stanModel1)$obsSD, breaks = 10, main="", xlab=expression(sigma^obs), xlim=c(0,1)) ; abline(v=0.2, lwd=2, lty=2)

par(mfrow=(c(1,1)))

betterPairs(stanModel1,  pars = c('r', "K"))

# calculate kernel density
v1 <- kernel_density_aprox(stanModel1, pars=c("r", "K"))


#--------------------------------------------------------
# fitting the logistic map STAN model with segmentation
#--------------------------------------------------------

sM <- stan_model("model2.stan")
segmentSize <- 5

inits <- list(
  list("r" = 2, chain_id=1, "K"=1.5, "obsSD"=0.3, "pop"=rep(0.5, nObs/segmentSize), "procerr"=matrix(rep(1, nObs-nObs/segmentSize), nrow=segmentSize-1)),
  list("r" = 3, chain_id=2, "K"=1, "obsSD"=0.3,  "pop"=rep(0.5, nObs/segmentSize), "procerr"=matrix(rep(1, nObs-nObs/segmentSize), nrow=segmentSize-1)),
  list("r" = 4, chain_id=3, "K"=1, "obsSD"=0.3, "pop"=rep(0.5, nObs/segmentSize), "procerr"=matrix(rep(1, nObs-nObs/segmentSize), nrow=segmentSize-1))
)


data = list(observed = matrix(trueObs, nrow=segmentSize), steps=nObs/segmentSize, substeps = segmentSize)

ptm <- proc.time()

stanfit = sampling(sM, data=data, iter=nIter, init=inits,chains = 3, seed=1337, thin= 30)
(StanTime <- proc.time() - ptm)

params2 <- c("r", "K", "obsSD", "pop[1]")

stanfit 

par(mfrow=c(1,1))
plot(stanfit)

rstan::traceplot(stanfit, pars=params2)


median_stan(stanfit, params2)


par(mfrow=(c(2,3)))
matplot(seq(1,by=30, length.out = (((nIter+warm)/30)+1)/3),matrix(extract(stanfit)$r, ncol=3), type="l", xlab="Iteration", ylab="r", lty=1, col=c("#88888888"))
matplot(seq(1,by=30, length.out = (((nIter+warm)/30)+1)/3),matrix(extract(stanfit)$K, ncol=3), type="l", xlab="Iteration", ylab="K", lty=1, col=c("#88888888"))
matplot(seq(1,by=30, length.out = (((nIter+warm)/30)+1)/3),matrix(extract(stanfit)$K, ncol=3), type="l", xlab="Iteration", ylab=expression(sigma^'obs'), lty=1, col=c("#88888888"))

hist(extract(stanfit)$r, breaks = 10, main="", xlab="r", xlim=c(3,3.8)) ; abline(v=3.7, lwd=2, lty=2)
hist(extract(stanfit)$K, breaks = 10, main="", xlab="K", xlim=c(0.8,1.2)) ; abline(v=1, lwd=2, lty=2)
hist(extract(stanfit)$obsSD, breaks = 10, main="", xlab=expression(sigma^obs), xlim=c(0.1,0.5)) ; abline(v=0.2, lwd=2, lty=2)

par(mfrow=(c(1,1)))

betterPairs(extract(stanfit,  pars = c('r', "K", "obsSD")))

v2 <- kernel_density_aprox(stanModel2, pars=c("r", "K"))


lyapunov(r=v2[[1]][1],K=v2[[2]][1])


sm2_extr <- extract(stanModel2, params2)

plot(density(sm2_extr$r))
plot(density(sm2_extr$K))
plot(density(sm2_extr$`pop[1]`))

#--------------------------------------------------------
# Synthetic Likelihood model for logistic map
#--------------------------------------------------------

#simmulator for logistic map
logisticSimul <- function(param, nsim, extraArgs, ...){
  if (!all(c("nObs", "burn") %in% names(extraArgs))) 
    stop("extraArgs should contain burn and nObs")
  if(length(param)<2)stop("no r and K submitted")
  burn <- extraArgs$burn
  nObs <- extraArgs$nObs
  if (is.null(extraArgs$p0) && length(param < 3)) 
    p0 = runif(1)
  else if (is.null(extraArgs$p0)) p0 = runif(1)
  else {p0 = param[4]; burn=0}
  if (is.null(extraArgs$trueprocerror)) trueprocerror = 0.005
  else trueprocerror <- extraArgs$trueprocerror
  simuls <- numeric()
  if (!is.vector(param)) 
    param <- simplify2array(param)

  simuls <- replicate(nsim, rlnorm(n = nObs, meanlog = log(pmax(0.0000001,population(r=param[1], K=param[2], p0=p0, burn=burn, procSD = trueprocerror,nObs = nObs, lag = 1, FUN= "logistic"))), sdlog = param[3])) 
  return(t(simuls))
}


#summary statistics simmilar to Wood (2010)
logisticStats <- function (x, extraArgs, ...) 
{
  obsData <- as.vector(extraArgs$obsData)
  stopifnot(length(obsData) != 0)
  if (!is.matrix(x)) 
    x <- matrix(x, 1, length(x))
  tx <- t(x)
  X0 <- t(orderDist(tx, obsData, np = 3, diff = 1))
  X0 <- cbind(X0, t(nlar(tx^0.3, lag = c(1, 1), power = c(1,2))))
  X0 <- cbind(X0,rowMeans(x))
  X0
}


#Model 
logistic_sl <- synlik(simulator = logisticSimul,
                      summaries = logisticStats,
                      param = c( r = 3.7, K = 1, obsSD=obsSD ),
                      extraArgs = list("nObs" = 50, "burn" = 50)
)

logistic_sl@data <- trueObs

logistic_sl@extraArgs$obsData <- logistic_sl@data

checkNorm(logistic_sl)


cl <- registerCluster(5)

ptm <- proc.time()
logistic_smcmc <- smcmc(logistic_sl, 
                        initPar = c("r"=3.5, "K"=0.8, "ObsSD"=0.01), 
                        niter = 10000, 
                        burn = 1000,
                        priorFun = function(input, ...) {
                          if(input[1]<0) input[1] <- Inf
                          if(input[1]>4) input[1] <- Inf
                          if(input[2]<0) input[2] <- Inf
                          if(input[3]<=0) input[3] <- Inf
                          if(input[3]>=0.5) input[3] <- Inf
                          sum(input)}, 
                        propCov = diag(c(0.1, 0.1, 0.1))^2, 
                        nsim = 500,
                        cluster=cl,
                        multicore=T)
(SynlikTime <- proc.time() - ptm)
stopCluster(cl)


addline1 <- function(parNam, ...) abline(h = logistic_smcmc@param[parNam], lwd = 2, lty = 2, col = 3) 
addline2 <- function(parNam, ...) abline(v = logistic_smcmc@param[parNam], lwd = 2, lty = 2, col = 3)

plot(logistic_smcmc, addPlot1 = "addline1", addPlot2 = "addline2")

# ABC rejection sampling

logistic_prior=list(c("unif",1,4),c("unif",0,1.2),c("unif",0.1,0.5))

logistic_abc_sum <- function(x) logisticStats(x, extraArgs = list(obsData=trueObs))
logistic_abc_model <- function(y){
  as.vector(logisticStats(logisticSimul(param=y, nsim=1, extraArgs= list("nObs" = 50, "burn" = 50)), extraArgs = list(obsData=trueObs)))
}
#logistic_abc_model(c(4,0.6,0.5))
target <- as.vector(logistic_abc_sum(trueObs))

n <- 1000000
p <- 0.1
logistic_abc <- ABC_rejection(model=logistic_abc_model, prior=logistic_prior, nb_simul=n, summary_stat_target=target,tol=p, progress_bar = T, n_cluster=1, use_seed = T)

hist(logistic_abc$param[,1], main="", xlab="r")
hist(logistic_abc$param[,2], main="", xlab="K")
hist(logistic_abc$param[,3], main="", xlab=expression(sigma^{obs}))

rej<-abc(target, logistic_abc$param, logistic_abc$stats, tol=0.1, method="rejection")
boxplot(rej$unadj.values)


hist(rej$unadj.values[,1], main="", xlab="r")
hist(rej$unadj.values[,2], main="", xlab="K")
hist(rej$unadj.values[,3], main="", xlab=expression(sigma^{obs}))

##-------------------------------
## Model Validation logistic map
##-------------------------------


## create Data

set.seed(1337)
nTests <- 100
r <- 3.7
K <- 1
burn <- 200
nObs <- 50
p0 <- expression(runif(1))
procSD <- 0.005
fun="logistic"
obsSD  <- sqrt(1/25.5)

test_states <- replicate(nTests, population(r=r, K=K, p0=eval(p0), procSD=procSD, burn=burn, nObs=nObs,lag=1,FUN=fun, rf=1))
testset <- matrix(rlnorm(length(test_states), meanlog=log(test_states), obsSD), nrow= nObs)

matplot(testset, type="l", xlab="t", ylab="n", col="#333300a0")



##------------------
## fit models
##------------------
cl <- registerCluster(n = 16)

stanlist <- foreach(i = 1:ncol(testset), .packages="rstan") %dopar% { 
  
  rstan_options(auto_write = TRUE)
  options(mc.cores = 1)
  
  obs <- testset[,i]
  data = list(observed = matrix(obs, nrow=segmentSize), steps=nObs/segmentSize, substeps = segmentSize)
  stanfit = sampling(sM, data=data, iter=nIter, init=inits,chains = 3, seed=1337, thin= 30)
  stanfit
}


stopCluster(cl)


## analyze results


parsNames <- c("r", "K", "obsSD")
pars <- sapply(stanlist, median_stan, parameters=parsNames)

par(mfrow=c(1,3))
for(i in 1:length(parsNames)) hist(as.numeric(pars[i,]), xlab=parsNames[i])
par(mfrow=c(1,1))


cl <- registerCluster(n = detectCores())

synliklist <- foreach(i = 1:ncol(testset), .packages="synlik") %dopar% { 
  
  logistic_sl@data <- logistic_sl@extraArgs$obsData <- testset[,i]
  
  logistic_smcmc <- smcmc(logistic_sl, 
                          initPar = c("r"=3.5, "K"=0.8, "ObsSD"=0.01), #
                          niter = 2500, 
                          burn = 100,
                          priorFun = function(input, ...) {
                            if(input[1]<0) input[1] <- Inf
                            if(input[1]>4) input[1] <- Inf
                            if(input[2]<0) input[2] <- Inf
                            if(input[3]<=0) input[3] <- Inf
                            if(input[3]>=0.5) input[3] <- Inf
                            sum(input)}, 
                          propCov = diag(c(0.1, 0.1, 0.1))^2, 
                          nsim = 250
  )
  logistic_smcmc
}

stopCluster(cl)


## analyze synlik results

parsSynlik <- sapply(synliklist, function(y) apply(y@chains, 2, median))

par(mfrow=c(1,3))
for(i in 1:length(parsNames)) hist(as.numeric(parsSynlik[i,]), xlab=parsNames[i], main = NULL)
par(mfrow=c(1,1))


## create Plot

lims <- list(r=c(2,4), K=c(0.5,1.5), obsSD=c(0.1,0.5))
par(mfrow=c(2,3))

for(i in 1:length(parsNames)){ 
  hist(as.numeric(pars[i,]), xlab=parsNames[i], main=NULL, xlim=lims[[i]], breaks = seq(lims[[i]][1], lims[[i]][2], length.out = 20))
  abline(v=c(r,K, obsSD)[i], col=2, lty=2)
}
for(i in 1:length(parsNames)){
  hist(as.numeric(parsSynlik[i,]), xlab=parsNames[i], main = NULL, xlim=lims[[i]], breaks = seq(lims[[i]][1], lims[[i]][2], length.out = 20))
  abline(v=c(r,K, obsSD)[i], col=2, lty=2)
}

##--------------------------------
## Non-parametric
##--------------------------------

# fit s-map model

for(i in 1:10){
  if(i==1) k <- i
  opt <-  optimize(f= function(x) NltsFn(trueObs, PredInterval = 1, Nembed=i,  Theta=x, Method = "Smap")$Corr *-1, interval= c(0,50))
  if(i==1) {opt2 <- opt} 
  if(opt2$objective>opt$objective) {opt2 <- opt; k <- i}
  
}  

logistic_smap <- NltsFn(trueObs, PredInterval = 1, Nembed=k,  Theta=opt2$minimum, Method = "Smap")
matplot(logistic_smap$Table, type="l")
lines(trueObs)

#fit gam
logistic_gam_df <- data.frame(obs=trueObs, par=c(NA,head(trueObs,-1)))
for(i in 2:k){logistic_gam_df <- cbind(logistic_gam_df, c(rep(NA,i),head(trueObs,-i)))}
names(logistic_gam_df) <- c("obs", paste0("par", 1:k))

logistic_gam <- gam(obs ~ s(par1, bs="cr") + s(par2, bs="cr") + s(par3, bs="cr") + s(par4, bs="cr") + s(par5, bs="cr"), data=logistic_gam_df)

logistic_gam_pred <- c(rep(NA,k-1),predict(logistic_gam))
lines(logistic_gam_pred, col=3)


# compare gam ans S-map
cor(tail(trueObs,-k), tail(logistic_smap$Table[,2],-k))
cor(tail(trueObs,-k), tail(as.numeric(logistic_gam_pred)[-51],-k))

nrmse(tail(truePop,-k), tail(logistic_smap$Table[,2],-k))
nrmse(tail(truePop,-k), tail(as.numeric(logistic_gam_pred)[-51],-k))

nrmse(tail(trueObs,-k), tail(logistic_smap$Table[,2],-k))
nrmse(tail(trueObs,-k), tail(as.numeric(logistic_gam_pred)[-51],-k))

p <- trueObs[1:k]
for(i in 1:50) p <- c(p, predict(logistic_gam, newdata = data.frame(par1=p[length(p)-(k-1)],par2=p[length(p)-(k-2)],par3=p[length(p)-(k-3)],par4=p[length(p)-(k-4)],par5=p[length(p)])))
plot(trueObs, type="l")
lines(p, col=4)

#fit gam with time lag = 1 
logistic_gam_df1 <- data.frame(obs=trueObs, par=c(NA,head(trueObs,-1)))
logistic_gam1 <- gam(obs ~ s(par, bs="cr"), data=logistic_gam_df1)

logistic_gam_pred1 <- c(NA,predict(logistic_gam1))

# compare gam ans S-map
cor(tail(trueObs,-k), tail(logistic_smap$Table[,2],-k))
cor(tail(trueObs,-k), tail(as.numeric(logistic_gam_pred)[-51],-k))
cor(tail(trueObs,-k), tail(as.numeric(logistic_gam_pred1)[-51],-k))

nrmse(tail(truePop,-k), tail(logistic_smap$Table[,2],-k))
nrmse(tail(truePop,-k), tail(as.numeric(logistic_gam_pred)[-51],-k))
nrmse(tail(truePop,-k), tail(as.numeric(logistic_gam_pred1)[-51],-k))

nrmse(tail(trueObs,-k), tail(logistic_smap$Table[,2],-k))
nrmse(tail(trueObs,-k), tail(as.numeric(logistic_gam_pred)[-51],-k))
nrmse(tail(trueObs,-k), tail(as.numeric(logistic_gam_pred1)[-51],-k))
