################################################################################
# Code for reply to Perretti et al
#
# Author: Florian Hartig, http://florianhartig.wordpress.com
#
# This scripts runs a single inference with Jags for the two model specifications
# and creates summary plots of the MCMC chains
#
# Requires data creation and modelSpecification to be run before
#
################################################################################


rm(list=ls(all=TRUE))

library(rjags) 
library(modeest)
library(IDPmisc)

# assumes data was stored in c:/temp
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

# betterpairs specifies nicer pair density plots
# good for looking at correlations in mcmc chains
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

betterpairs<-function(data)pairs(as.matrix(lumpchains(data)), lower.panel=function(...){par(new=TRUE);ipanel.smooth(..., nrpoints=0)}, diag.panel=panel.hist, upper.panel=panel.cor)


# fitting the original model

inits <- list(
  list("r" = 2.2, "K"=1, "N0" = 0.5, "ObservationPrecision"=4),
  list("r" = 2.5, "K"=1, "N0" = 0.5, "ObservationPrecision"=4),
  list("r" = 2.8, "K"=1, "N0" = 0.5, "ObservationPrecision"=4)
)

data = list(observed = observed[1:(sampleSize+1)], nObs=sampleSize)
params = c("r", "K", "obsSD", "pop[1]") # Parameters to be monitored

ptm <- proc.time()

jagsModel = jags.model( file= "model1.txt", data=data, n.chains = 3, n.adapt= 500000, inits = inits)
update(jagsModel, 500000)
results = coda.samples( jagsModel , variable.names=params,n.iter=5000000, thin = 500)

proc.time() - ptm

png("traceplots-model1.png", width = 1500, height = 1500, pointsize = 30)
plot(results)
dev.off()

png("correlation-model1.png", width = 1000, height = 1000, pointsize = 30)
betterpairs(results)
dev.off()


png("gelmanplots-model1.png", width = 1000, height = 1000, pointsize = 30)
gelman.plot(results)
mtext(paste("m-psrf",round(gelman.diag(results)$mpsrf, digits=3)), side = 3, outer = F, padj=-3)
dev.off()


# fitting the alternative model
data = list(observed = matrix(observed[1:sampleSize], nrow=5), steps=10, substeps = 5)
params = c("r", "K", "obsSD", "pop[1,1]") # Parameters to be monitored

ptm <- proc.time()

jagsModel = jags.model( file= "model2.txt", data=data, n.chains = 3 , n.adapt= 20000, inits = inits)
update(jagsModel, 100000)
results = coda.samples( jagsModel , variable.names=params,n.iter=300000, thin = 30)

proc.time() - ptm


png("traceplots-model2.png", width = 1500, height = 1500, pointsize = 30)
plot(results)
dev.off()

png("correlation-model2.png", width = 1000, height = 1000, pointsize = 30)
betterpairs(results)
dev.off()


png("gelmanplots-model2.png", width = 1000, height = 1000, pointsize = 30)
gelman.plot(results)
mtext(paste("m-psrf",round(gelman.diag(results)$mpsrf, digits=3)), side = 3, outer = F, padj=-3)
dev.off()

