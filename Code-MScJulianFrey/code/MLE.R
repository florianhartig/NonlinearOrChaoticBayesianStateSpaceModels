source("Packages_Functions_DataCreation.R")

#--------------------------------------------------------
# Fitting MLE without segmentation
#--------------------------------------------------------

# grid approximation for R 

ptm <- proc.time()
R <- seq(2,3.9,10^-4)
llh_par_r <- sapply(X=R, FUN=likelihood, p0=truePop[1])
proc.time() -ptm


plot(R,llh_par_r, type="l", ylab="Log-Likelihood")
points(R[which.max(llh_par_r)], max(llh_par_r))

llh <- likelihood
likelihood <- function(x) llh(r=x[1], K=x[2], p0 = x[3], sd=x[4], procERR=x[5:length(x)])
optimizer_results <- data.frame(Method=numeric(), r=numeric(), K=numeric(), p0=numeric(), sd=numeric())

# approximate r using samplers
optimizer_results[1,] <- c("optim_SANN", optim(par=c(r=3, K=0.8, p0=0.5, sd=0.1), fn=likelihood, method="SANN" , control= list(fnscale=-1, maxit=100000, temp=100, tmax=10))$par)
optimizer_results[2,] <- c("optim_Nelder_Mead",optim(par=c(r=3, K=0.8, p0=0.5, sd=0.1), fn=likelihood, method="Nelder-Mead" , control= list(fnscale=-1, maxit=50000))$par)
optimizer_results[3,] <- c("optim_BFGS",optim(par=c(r=3, K=0.8, p0=0.5, sd=0.1), fn=likelihood, method="BFGS" , control= list(fnscale=-1, maxit=50000))$par)
optimizer_results[4,] <- c("optim_CG",optim(par=c(r=3, K=0.8, p0=0.5, sd=0.1), fn=likelihood, method="CG" , control= list(fnscale=-1, maxit=50000))$par)
optimizer_results[5,] <- c("optim_L-BFGS-B",optim(par=c(r=3, K=0.8, p0=0.5, sd=0.1), fn=likelihood, method="L-BFGS-B" , control= list(fnscale=-1, maxit=50000))$par)
                                                                                                                                       

ptm <- proc.time()
gen_r <- GenSA(fn=function(x){-1*likelihood(x)}, lower=c(r=0,K=0,p0=0, sd=0), upper=c(r=4,K=1.5,p0=1.5,sd=1), control = list(temperature=1000, nb.stop.improvement=10000))
str(gen_r)
proc.time() -ptm
optimizer_results[6,] <- c("GenSA",gen_r$par)

optimizer_results1 <- optimizer_results

cl <- registerCluster(n=2)


optimizer_results2<- foreach(i=1:50, .packages = "GenSA",.combine="rbind") %dopar%{
  obsSD <- sqrt(1/25.5)
  r <- 3.7
  K <- 1
  nObs = 50
  
  truePop <- population(r=r, K= K, p0 = 0.5, procSD = 0.005, burn = 100, nObs=nObs)
  trueObs <- rlnorm(length(truePop), mean = log(truePop) , sd = obsSD)
  gen_r <- GenSA(fn=function(x){-1*likelihood(x)}, lower=c(r=0,K=0,p0=0, sd=0, procERR=rep(0.9, nObs-1)), upper=c(r=4,K=1.5,p0=1.5,sd=1, procERR=rep(1.1, nObs-1)), control = list(temperature=1000, nb.stop.improvement=10000))
  gen_r$par
}

optimizer_results_NM2 <- foreach(i=1:50, .packages = "GenSA",.combine="rbind") %dopar%{
  obsSD <- sqrt(1/25.5)
  r <- 3.7
  K <- 1
  nObs = 50
  
  truePop <- population(r=r, K= K, p0 = 0.5, procSD = 0.005, burn = 100, nObs=nObs)
  trueObs <- rlnorm(length(truePop), mean = log(truePop) , sd = obsSD)
  optim(par=c(r=3, K=0.8, p0=0.5, sd=0.1, procERR=rep(1, nObs-1)), fn=likelihood, method="Nelder-Mead" , control= list(fnscale=-1, maxit=50000))$par
}
stopCluster(cl)

xlim=list(c(1,4),c(0,8),c(0,1.6),c(0,0.5))
par(mfrow=c(2,4))
for(i in 1:4) {hist(optimizer_results_NM[,i],breaks=seq(xlim[[i]][1],xlim[[i]][2], length.out = 20) ,xlab=c("r","K",expression(sigma^{obs}),expression("p"[1]))[i], main="", xlim=xlim[[i]]) 
  abline(v=c(3.7,1,obsSD,0.5)[i], lty
         =2, lwd=2, col=2)}

for(i in 1:4) {hist(optimizer_results[,i],breaks=seq(xlim[[i]][1],xlim[[i]][2], length.out = 20), xlab=c("r","K",expression(sigma^{obs}),expression("p"[1]))[i], main="", xlim=xlim[[i]]) 
  abline(v=c(3.7,1,obsSD,0.5)[i], lty
         =2, lwd=2, col=2)}

par(mfrow=c(1,1))


#--------------------------------------------------------
# Fitting MLE with segmentation
#--------------------------------------------------------

# grid approximation for R 

ptm <- proc.time()
llh_par_r <- sapply(X=R, FUN=llh, p0=truePop[seq(1,nObs,5)], segmentSize=5) #truePop[seq(1,nObs,5)]
proc.time() -ptm

plot(R,llh_par_r, type="l", ylab="Log-Likelihood")
points(R[which.max(llh_par_r)], max(llh_par_r))



likelihood <- function(x) {llh(r=x[1], K=x[2], sd=x[3], p0=x[4:13], segmentSize=5)}

optimizer_results <- data.frame(Method=numeric(), r=numeric(), K=numeric(), sd=numeric(), p01=numeric(), p02=numeric(), p03=numeric(), p04=numeric(), p05=numeric(), p06=numeric(), p07=numeric(), p08=numeric(), p09=numeric(), p10=numeric())
#optimize(likelihood,interval=c(0,4), maximum=T, tol=10^-12)
optimizer_results[1,] <- c("optim_SANN", optim(par=c(r=3, K=1.2, sd=0.1, p0=rep(0.5,10)), fn=likelihood, method="SANN" , control= list(fnscale=-1, maxit=100000, temp=100, tmax=10))$par)
#optim(par=c(r=3),lower=0, upper=4, fn=likelihood, method="Brent" , control= list(fnscale=-1, maxit=50000), p0=truePop[seq(1,nObs,5)], segmentSize=5)
optimizer_results[2,] <- c("optim_Nelder_Mead",optim(par=c(r=3, K=1.2, sd=0.1, p0=rep(0.5,10)), fn=likelihood, method="Nelder-Mead" , control= list(fnscale=-1, maxit=50000))$par)
optimizer_results[3,] <- c("optim_BFGS",optim(par=c(r=3, K=1.2, sd=0.1, p0=rep(0.5,10)), fn=likelihood, method="BFGS" , control= list(fnscale=-1, maxit=50000))$par)
optimizer_results[4,] <- c("optim_CG",optim(par=c(r=3, K=1.2, sd=0.1, p0=rep(0.5,10)), fn=likelihood, method="CG" , control= list(fnscale=-1, maxit=50000))$par)
optimizer_results[5,] <- c("optim_L-BFGS-B",optim(par=c(r=3, K=1.2, sd=0.1, p0=rep(0.5,10)), fn=likelihood, method="L-BFGS-B" , control= list(fnscale=-1, maxit=50000))$par)

likelihood <- function(x) {-1*llh(r=x[1], K=x[2], sd=x[3], p0=x[4:13], procERR=x[14:length(x)], segmentSize=5)}

ptm <- proc.time()
gen_r <- GenSA(par=c(r=3, K=1.2, sd=0.1, p0=rep(0.5,10)), fn=likelihood, lower=c(r=0.1,K=0.1, sd=0.1 ,p0= rep(0.1, 10)), upper=c(r=4,K=1.5, sd=1,p0= rep(1.5, 10)), control = list(temperature=1000, nb.stop.improvement=10000))
str(gen_r$par)
proc.time() -ptm
optimizer_results[6,] <- c("GenSA",gen_r$par)

optimizer_results


cl <- registerCluster(n=2)
optimizer_results<- foreach(i=1:50, .packages = "GenSA",.combine="rbind") %dopar%{
  obsSD <- sqrt(1/25.5)
  r <- 3.7
  K <- 1
  nObs = 50
  
  truePop <- population(r=r, K= K, p0 = 0.5, procSD = 0.005, burn = 100, nObs=nObs)
  trueObs <- rlnorm(length(truePop), mean = log(truePop) , sd = obsSD)
  gen_r <- GenSA(par=c(r=3, K=1.2, sd=0.1, p0=rep(0.3,10), procERR=rep(1,40)), fn=likelihood, lower=c(r=0.1,K=0.1, sd=0.1 ,p0= rep(0.1, 10), procERR=rep(0.9,40)), upper=c(r=4,K=1.5, sd=1,p0= rep(1.5, 10), procERR=rep(1.1,40)), control = list(temperature=1000, nb.stop.improvement=10000))
  gen_r$par
}
optimizer_results_NM <- foreach(i=1:50, .packages = "GenSA",.combine="rbind") %dopar%{
  obsSD <- sqrt(1/25.5)
  r <- 3.7
  K <- 1
  nObs = 50
  
  truePop <- population(r=r, K= K, p0 = 0.5, procSD = 0.005, burn = 100, nObs=nObs)
  trueObs <- rlnorm(length(truePop), mean = log(truePop) , sd = obsSD)
  optim(par=c(r=3, K=1.2, sd=0.1, p0=rep(0.5,10), procERR=rep(1,40)), fn=likelihood, method="Nelder-Mead" , control= list(maxit=50000))$par
}
stopCluster(cl)

xlim=list(c(2,4),c(0.6,1.4),c(0.1,0.5),c(0,1.2))
par(mfrow=c(2,4))
for(i in 1:4) {hist(optimizer_results_NM[,i],breaks=seq(xlim[[i]][1],xlim[[i]][2], length.out = 20) ,xlab=c("r","K",expression(sigma^{obs}),expression("p"[1]))[i], main="", xlim=xlim[[i]]) 
abline(v=c(3.7,1,obsSD,0.5)[i], lty
       =2, lwd=2, col=2)}

for(i in 1:4) {hist(optimizer_results[i,],breaks=seq(xlim[[i]][1],xlim[[i]][2], length.out = 20), xlab=c("r","K",expression(sigma^{obs}),expression("p"[1]))[i], main="", xlim=xlim[[i]]) 
  abline(v=c(3.7,1,obsSD,0.5)[i], lty
         =2, lwd=2, col=2)}

par(mfrow=c(1,1))

# Bifurcation Plot

n <- 1
R <- seq(0,4,length=500)
f <- expression(a*x*(1-x))
data <- matrix(0,200,501)

X <- seq(0,1,length.out = 501)

for(a in R){
  x <- X[n]
  #x <- runif(1) # random initial condition
  ## first converge to attractor
  for(i in 1:200){
    x <- eval(f)
  } # collect points on attractor
  for(i in 1:200){
    x <- eval(f)
    data[i,n] <- x
  }
  n <- n+1
}

data <- data[,1:500]
plot(R,data[1,], pch=14,  ylab="", xlab="", axes=F, col="#888888")
for(i in 2:200) points(R,data[i,],pch=24, col="#888888")
