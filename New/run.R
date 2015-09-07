################################################################################
# 
#
# Author: Florian Hartig, http://florianhartig.wordpress.com
#
# 
#
################################################################################

source("New/1-settings.R")
source("New/models/logisticModel.r")
source("New/models/observation.R")


truevalues <- logisticModel(K=1,r=3.7,N0=0.7, samplesize=20, processerror=0.005)
observedValues <- observationError(data = truevalues)


plot(observedValues, type = "b", col = "darkred", pch = 2, cex = 1.2, lwd = 2, ylim = c(0.2,1.2), xlab = "", ylab = "N/K", lty = 2)
legend("topright", legend = c("Observed population size"), lty = 2, pch = 2, lwd = 2, col = "darkred")

par(cex =1.2)

plot(truevalues, type = "b", lty = 2, lwd =2 , ylim = c(0.2,1.2), xlab = "", ylab = "N/K")
legend("topright", legend = c("True population size"), lty = 2, pch = 1, lwd = 2)


for (i in 1:200){
  truevaluesTemp <- logisticModel(K=1,r=3.7,N0 = rnorm(1,0.7,0.02), samplesize=20, processerror=0.005)
  lines(truevaluesTemp,  col = "#11111111", lwd = 0.5)
}
points(truevalues, type = "b", lty = 2, lwd =2 )



for (i in 1:200){
  observedValuesTemp <- observationError(data = truevalues)  
  lines(observedValuesTemp,  col = "#00990011")
}

lines(truevalues, lty = 2, lwd =2 )
points(truevalues, lwd =2)
legend("topright", legend = c("True population size"), lty = 2, pch = 1, lwd = 2)



points(observedValues,  col = "darkred", pch = 2, cex = 1.2, lwd = 3)
legend("topright", legend = c("True population size", "Observed population size"), lty = c(2,0), pch = c(1, 2), lwd = 2, col = c("black", "darkred"))



#################################################


truevalues <- logisticModel(K=1,r=3.7,N0=0.7, samplesize=20, processerror=0.005)
observedValues <- observationError(data = truevalues)
alternative = rep(0.6,20)


par(cex = 1.2, mar = c(2.5,2.5,2.5,2.5))

plot(observedValues,  col = "darkred", pch = 2, cex = 1, lwd = 2, ylim = c(0,1.2))


  for (i in 1:400){
    truevaluesTemp <- logisticModel(K=1,r=3.7,N0 = rnorm(1,observedValues[1],0.005), samplesize=200, processerror=0.005)
    #lines(observationError(truevaluesTemp),  col = "#00990011", lwd = 0.5)
    lines(truevaluesTemp,  col = "#44444411", lwd = 1)
  }

points(observedValues,  col = "darkred", pch = 2, cex = 1, lwd = 2)



plot(observedValues,  col = "darkred", pch = 2, cex = 1, lwd = 2, ylim = c(0,1.2))
for (i in 1:400){
  lines(observationError(rep(0.5,20), param = c(0.3)) , lwd =1, col = "#00000022" )
}
#abline(h= 0.6, col = "blue", lwd = 3)
points(observedValues,  col = "darkred", pch = 2, cex = 1, lwd = 2)




plot(observedValues,  col = "darkred", pch = 2, cex = 1, lwd = 2, ylim = c(0,1.2))
for (j in c(1,6,11,16)){
  
  for (i in 1:400){
    truevaluesTemp <- logisticModel(K=1,r=3.7,N0 = rnorm(1,observedValues[j],0.01), samplesize=5, processerror=0.005)
    lines(j:(j+4), truevaluesTemp,  col = "#11111111", lwd = 0.5)
  }  
}
points(observedValues,  col = "darkred", pch = 2, cex = 1, lwd = 2)


abline(v = c(5.5, 10.5, 15.5))



## chopping up annimation 






