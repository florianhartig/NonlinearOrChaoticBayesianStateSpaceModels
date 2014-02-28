################################################################################
# Code for reply to Perretti et al
#
# Author: Florian Hartig, http://florianhartig.wordpress.com
#
# This script draws the plots for the paper
# 
# Requires the results of run inference for all data to be run before and data
# being stored in results files
#
################################################################################

rm(list=ls(all=TRUE))

set.seed(3)

library(compiler)


# collected input data from different model runs 

setwd("D:/home/Projekte/Papers_in_Progress/13-ResponseToPerrettiEtAl/temp")

load(file = "results21.RData")
load(file = "results22.RData")

# runs were split across two processors
load(file = "results111.RData")
model2resultsfixed1 = model2resultsfixed
load(file = "results112.RData")
for (i in 1:3) model2resultsfixed[[i]] = rbind(model2resultsfixed1[[i]], model2resultsfixed[[i]])
model2resultsfixed[[4]] = c(model2resultsfixed1[[4]], model2resultsfixed[[4]])

load(file = "results121.RData")
model2resultsvary1 = model2resultsvary
load(file = "results122.RData")
for (i in 1:3) model2resultsvary[[i]] = rbind(model2resultsvary1[[i]], model2resultsvary[[i]])
model2resultsvary[[4]] = c(model2resultsvary1[[4]], model2resultsvary[[4]])




setwd("D:/home/Projekte/Papers_in_Progress/13-ResponseToPerrettiEtAl/article")



conv = 1.2  # mpsrf convergence criterion

pop.sd=function(x)sqrt(mean(((mean(x) - x))^2))




# function to calculate and plot different metrics for error
# error 3 is the version used by Perretti et al.
ploterrors <- function(name, parameters, use){
  
  numpars <- dim(parameters)[1]
  error <- matrix(nrow = numpars, ncol = sampleSize) 
  names <- gsub(" ", "", name)
  names <- gsub("\\.", "", names) 

  
  # calculates forecast at quadradic error
  pdf(paste("pred-uncertainty-",names,".pdf", sep = ""),  width=5, height=5) 
  #cex = 1.3
  plot(pop[51:100],type="n",ylim =c(0,1), xlab = "Time interval", ylab = "Population size", main = name)
  for (i in 1:numpars){
    observed <- truepops[i,51:100]
    predicted <- logistic(parameters[i,1],parameters[i,4],truepops[i,51],sampleSize)  
    error[i,] <- (predicted - observed)^2
    if (i < 20) lines(logistic(trueK,trueR,truepops[1,51],sampleSize), col="#00000020") 
  }
  lines(logistic(trueK,trueR,truepops[1,51],sampleSize))
  dev.off()
  
  
  # calculates forecast at quadradic error
  png(paste("pred-uncertainty-",names,".png", sep = ""))
  cex = 1.3
  plot(pop[51:100],type="n",ylim =c(0,1), xlab = "Time interval", ylab = "Population size", main = name)
  for (i in 1:numpars){
    observed <- truepops[i,51:100]
    predicted <- logistic(parameters[i,1],parameters[i,4],truepops[i,51],sampleSize)  
    error[i,] <- (predicted - observed)^2
    if (i < 20) lines(logistic(trueK,trueR,truepops[1,51],sampleSize), col="#00000020") 
  }
  lines(logistic(trueK,trueR,truepops[1,51],sampleSize))
  dev.off()
  
  
  
  prederror <- sqrt(apply(error[use,], 2, mean)) / sd(as.vector(truepops[use,51:100]))
  erroruncertainty <- apply(error[use,], 2, sd) / sd(as.vector(truepops[use,51:100]))
  
  NRMSE <- rep(NA, 49)
  NRMSESD <- rep(NA, 49)
  
  # calculates forecastig error for different intervals, replicates Perretti et al.
  error2 <- matrix(nrow = numpars, ncol = (sampleSize-2))
  for (i in 1:numpars){
    for (pred in 1:48){
      store = rep(NA, 50-pred)
      k = 1
      for (j in (51+pred):100){
        observed <- truepops[i,j]
        predicted <- logistic(parameters[i,1],parameters[i,4],truepops[i,j-pred],pred+1)[pred+1]
        store[k] <- (observed - predicted)^2
        k = k + 1
      }
      error2[i,pred]= sqrt(mean(store)) / pop.sd(truepops[i,(51+pred):100])
    }
  }
  NRMSE2 <- apply(error2[use,], 2, mean)
  NRMSESD2 <- apply(error2[use,], 2, pop.sd)
  
  modelerrors <- list("prederror" = prederror, "erroruncertainty" = erroruncertainty, "NRMSE" = NRMSE, "NRMSESD" = NRMSESD, "NRMSE2" = NRMSE2, "NRMSESD2" = NRMSESD2, "length" = length(use))
  
  # alternative calculation of RMSE
  NRMSEarray <- sqrt(error[,2]) / apply(truepops[,51:100], 1, sd) #
  modelerrors$NRMSE[1] <- mean(NRMSEarray)
  modelerrors$NRMSESD[1] <- sd(NRMSEarray)
  for (i in 2:49){
    NRMSEarray <- sqrt(apply(error[,2:(i+1)], 1, mean)) / sd(as.vector(truepops[,51:100])) #apply(truepops[,52:(50+i)], 1, sd) #
    modelerrors$NRMSE[i] <- mean(NRMSEarray)
    modelerrors$NRMSESD[i] <- sd(NRMSEarray)
    
  }

  
  ###################
  # Plots
  
  xval = 1:49
  
  #   pdf(paste("pred-error-", names, ".pdf", sep = ""), width=5, height=5)
  #   plot(modelerrors$prederror,type = "l", xlab = "Time interval", ylab = "Point-wise RMSE / Time series SD", main = name, ylim = c(0,2.5))
  #   lines(modelerrors$prederror+modelerrors$erroruncertainty/sqrt(50))
  #   lines(modelerrors$prederror-modelerrors$erroruncertainty/sqrt(50))
  #   lines(modelerrors$prederror+modelerrors$erroruncertainty, lty = 2)
  #   lines(modelerrors$prederror-modelerrors$erroruncertainty, lty = 2)
  #   abline(h=1)
  #   dev.off()
#   
#   pdf(paste("pred-error2-", names, ".pdf", sep = ""), width=5, height=5)
#   plot(xval,modelerrors$NRMSE,type = "l", xlab = "Forecasting interval", ylab = "NRMSE / Forecasting SD", main = name, ylim = c(0,1.6))
#   lines(xval,modelerrors$NRMSE+modelerrors$NRMSESD/sqrt(length(use)))
#   lines(xval,modelerrors$NRMSE-modelerrors$NRMSESD/sqrt(length(use)))
#   lines(xval,modelerrors$NRMSE+modelerrors$NRMSESD, lty = 2)
#   lines(xval,modelerrors$NRMSE-modelerrors$NRMSESD, lty = 2)
#   abline(h=1)
#   dev.off()
#   
  
  xval = 1:30
  
  pdf(paste("pred-error3-", names, ".pdf", sep = ""), width=5, height=5)
  plot(xval,modelerrors$NRMSE2[xval],type = "l", xlab = "Forecasting interval", ylab = "NRMSE / Forecasting SD", main = name, ylim = c(0,1.6))
  lines(xval,modelerrors$NRMSE2[xval]+modelerrors$NRMSESD2[xval]/sqrt(length(use)))
  lines(xval,modelerrors$NRMSE2[xval]-modelerrors$NRMSESD2[xval]/sqrt(length(use)))
  lines(xval,modelerrors$NRMSE2[xval]+modelerrors$NRMSESD2[xval], lty = 2)
  lines(xval,modelerrors$NRMSE2[xval]-modelerrors$NRMSESD2[xval], lty = 2)
  abline(h=1)
  dev.off()
  
  return (modelerrors)
}

ploterrors <- cmpfun(ploterrors)


name = "True parameters"
idpars <- matrix(rep(c(trueK, NA, NA, trueR), times = 50), ncol = 4, byrow = T)
errors0 <- ploterrors(name, parameters = idpars, use = rep(T,50))

name="Perretti et al."
errors2 <- ploterrors(name, parameters = model2resultsfixed$medians, use = model2resultsfixed$mpsrf<conv)

name = "Alternative method"
errors3 <- ploterrors(name, parameters = model3resultsfixed$medians, use = model3resultsfixed$mpsrf<conv)


# Plots for alternative error metrics
# 
# xval = 1:49
# 
# 
# pdf("errorcomparison.pdf", width=5, height=5)
#   oldpar <-par(cex.lab = 1.2)
#   plot(xval, errors2$NRMSE,type = "n", xlab = "Forecasting interval", ylab = "RMSE / Forecasting time-series SD", main = "Predictive error", ylim = c(0,1.5))
#   
#   polygon(c(xval,rev(xval)), c(errors2$NRMSE+errors2$NRMSESD/sqrt(errors2$length), rev(errors2$NRMSE-errors2$NRMSESD/sqrt(errors2$length))), col = "#8B000022", border = NA )
#   lines(xval, errors2$NRMSE, col = "darkred")
#   points(xval, errors2$NRMSE, col = "darkred", cex = .7)
# 
#   polygon(c(xval,rev(xval)), c(errors3$NRMSE+errors3$NRMSESD/sqrt(errors3$length), rev(errors3$NRMSE-errors3$NRMSESD/sqrt(errors3$length))), col = "#00640022", border = NA )
#   lines(xval, errors3$NRMSE, col = "darkgreen")
#   points(xval, errors3$NRMSE, col = "darkgreen", pch = 2, cex = .7)
#   
#   polygon(c(xval,rev(xval)), c(errors0$NRMSE+errors0$NRMSESD/sqrt(errors0$length), rev(errors0$NRMSE-errors0$NRMSESD/sqrt(errors0$length))), col = "#00008B22", border = NA )
#   lines(xval, errors0$NRMSE, col = "darkblue")
#   points(xval, errors0$NRMSE, col = "darkblue", pch = 4, cex = .7)
#   
#   abline(h=1)
#   
#   legend("bottomright", legend = c("Perretti et al.", "Alternative", "Identical"), col = c("darkred","darkgreen","darkblue"),lty = 1, pch = c(1,2,4), bg = "white" )
# dev.off()



xval = 1:30

pdf("errorcomparison2.pdf", width=5, height=5)
  oldpar <-par(cex.lab = 1.2)
  plot(xval, errors2$NRMSE2[xval],type = "p", xlab = "Forecasting interval", ylab = "RMSE / Forecasting time-series SD", main = "Predictive error", ylim = c(0,1.5), col = "darkred", cex = .8 )
  
  polygon(c(xval,rev(xval)), c(errors2$NRMSE2[xval]+errors2$NRMSESD2[xval]/sqrt(errors2$length)*1.96, rev(errors2$NRMSE2[xval]-errors2$NRMSESD2[xval]/sqrt(errors2$length)*1.96)), col = "#8B000022", border = NA )
  lines(xval, errors2$NRMSE2[xval], col = "darkred")
  points(xval, errors2$NRMSE2[xval], col = "darkred", cex = .8)

  polygon(c(xval,rev(xval)), c(errors3$NRMSE2[xval]+errors3$NRMSESD2[xval]/sqrt(errors3$length)*1.96, rev(errors3$NRMSE2[xval]-errors3$NRMSESD2[xval]/sqrt(errors3$length)*1.96)), col = "#00640022", border = NA )
  lines(xval, errors3$NRMSE2[xval], col = "darkgreen")
  points(xval, errors3$NRMSE2[xval], col = "darkgreen", pch = 2, cex = .8)
  
  polygon(c(xval,rev(xval)), c(errors0$NRMSE2[xval]+errors0$NRMSESD2[xval]/sqrt(errors0$length)*1.96, rev(errors0$NRMSE2[xval]-errors0$NRMSESD2[xval]/sqrt(errors0$length)*1.96)), col = "#00008B22", border = NA )
  lines(xval, errors0$NRMSE2[xval], col = "darkblue")
  points(xval, errors0$NRMSE2[xval], col = "darkblue", pch = 4, cex = .8)
  
  abline(h=1)
  mtext("B", side = 3, cex = 2.5, line=2.2,at=-5.5)
  

  legend("bottomright", legend = c("Perretti et al.", "Alternative method", "True parameters"), col = c("darkred","darkgreen","darkblue"),lty = 1, pch = c(1,2,4), bg = "white" )
dev.off()


# Histogram of the fitted R

pdf("hist-r-model0.pdf", width=5, height=5)
hist(rep(3.7,50), freq = T, xlab = "Growth rates r", ylab = "Frequency", main = "True parameters", xlim=c(2,4.5), 
     breaks = seq(from=2,to=4.5, length.out=100))
dev.off()

pdf("hist-r-model2.pdf", width=5, height=5)
hist(model2resultsfixed$medians[model2resultsfixed$mpsrf<conv,4], freq = T, xlab = "Posterior median of estimated growth rate r", ylab = "Frequency", main = "Perretti et al.", xlim=c(2,4.5), breaks = seq(from=2,to=4.5, length.out=100))
dev.off()

pdf("hist-r-model3.pdf", width=5, height=5)
hist(model3resultsfixed$medians[model3resultsfixed$mpsrf<conv,4], freq = T, xlab = "Posterior median of estimated growth rate r", ylab = "Frequency", main = "Alternative method", xlim=c(2,4.5), breaks = seq(from=2,to=4.5, length.out=100))
dev.off()
      

# Comparison r fitted true

pdf("true-estimated-r-model2.pdf", width=5, height=5)
par(cex.lab = 1.4, cex.axis = 1.4)
plot(trueRArray[model2resultsvary$mpsrf<conv], model2resultsvary$medians[model2resultsvary$mpsrf<conv,4], ylab = "Estimated growt rates r", xlab = "True growth rate r", main = "Perretti et al.", ylim = c(2,4))
abline(a=0,b=1)
dev.off()


pdf("true-estimated-r-model3.pdf", width=5, height=5)
par(cex.lab = 1.4, cex.axis = 1.4)
plot(trueRArray[model3resultsvary$mpsrf<conv], model3resultsvary$medians[model3resultsvary$mpsrf<conv,4], ylab = "Estimated growt rates r", xlab = "True growth rate r", main = "Alternative method", ylim = c(2,4))
abline(a=0,b=1)
dev.off()



pdf("true-estimated.pdf", width=5, height=5)

oldpar <-par(cex.lab = 1.2)

plot(trueRArray[model2resultsvary$mpsrf<conv], model2resultsvary$medians[model2resultsvary$mpsrf<conv,4], ylab = "Estimated growt rates r", xlab = "True growth rate r", main = "True vs. estimated", pch = 7, ylim = c(2,4))
points(trueRArray[model3resultsvary$mpsrf<conv], model3resultsvary$medians[model3resultsvary$mpsrf<conv,4], pch=4)
abline(a=0,b=1)

legend("bottomleft", legend = c("Perretti et al.", "Alternative method"), pch = c(7,4), bg = "white" )

mtext("A", side = 3, cex = 2.5, line=2.2,at=3.04)

dev.off()



# function used to calculate lyapunov exp, slightly modified from Peretti et al (personal communication)
lyapunov <- function(r){
  
  # This script calculate the lyapunov exponent for each value of r for the logistic
  # model.
  
  burn      <- 100000
  T         <- 1000
  x.burn    <- vector(mode="numeric", length=burn)
  x.burn[1] <- 0.5
  x.a       <- vector(mode="numeric", length=T)
  x.b       <- vector(mode="numeric", length=T)
  lyap      <- vector(mode="numeric", length=T-1)
  avlyap    <- vector(mode="numeric", length=T-1)
  K         <- 1.0
  tau.proc  <- 40000
  sd.proc   <- sqrt(1/tau.proc)
  d0        <- 10^-12 # Distance to separate orbits
  
  ## Process model
  # Burn-in
  for(t in 1:(burn-1)){
    x.burn[t+1] <- r*x.burn[t]*(1-x.burn[t]/K)
  }
  
  # Calculate lyapunov exponent
  x.a[1] <- x.burn[length(x.burn)] # Initiate on attractor
  x.b[1] <- x.burn[length(x.burn)] + d0
  
  for(t in 1:(T-1)){
    
    x.a[t+1] <- r*x.a[t]*(1-x.a[t]/K)
    x.b[t+1] <- r*x.b[t]*(1-x.b[t]/K)
    
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
lyap <- sapply(trueRArray, lyapunov)


bias2med <- trueRArray[model2resultsvary$mpsrf<conv] - model2resultsvary$medians[model2resultsvary$mpsrf<conv,4]
bias2mod <- trueRArray[model2resultsvary$mpsrf<conv] - model2resultsvary$modes[model2resultsvary$mpsrf<conv,4]

bias3med <- trueRArray[model3resultsvary$mpsrf<conv] - model3resultsvary$medians[model3resultsvary$mpsrf<conv,4]
bias3mod <- trueRArray[model3resultsvary$mpsrf<conv] - model3resultsvary$modes[model3resultsvary$mpsrf<conv,4]



pdf("Lyapunov-bias-r-model2.pdf", width=5, height=5)
par(cex.lab = 1.4, cex.axis = 1.4)
plot(lyap[model2resultsvary$mpsrf<conv], bias2med, xlab = "Lyapunov Exponent true model", ylab = "Bias of estimated r", main = "Perretti et al.", ylim = c(-0.2, 1.5))
abline(h=0)
dev.off()


pdf("Lyapunov-bias-r-model3.pdf", width=5, height=5)
par(cex.lab = 1.4, cex.axis = 1.4)
plot(lyap[model3resultsvary$mpsrf<conv], bias3med, xlab = "Lyapunov Exponent true model", ylab = "Bias of estimated r", main = "Alternative method", ylim = c(-0.2, 1.5))
abline(h=0)
dev.off()

