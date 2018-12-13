

# Here a model with truly nonlinear dynamics, including the option to controll process and observation error
logisticModel <- function(K = 1,r = 3,N0 = 0.2,samplesize = 100, processerror=0.1, obsError = 0){
  pop = rep(NA,(samplesize))
  pop[1] = N0
  for (i in 2:(samplesize)){
    pop[i] = pop[i-1] * r * (1 - pop[i-1]/K) * rlnorm(1,0,processerror)
  } 
  out= list()
  out$pop <- pmax(0,pop)
  out$popObs <- rlnorm(samplesize, mean = log(pop) , sd = obsError)
  return(out)
}

x = logisticModel()
plot(ts(x$pop))


# Here simulations of a standard autoregressive model, i.e. this is only uncorrelated noise 

t1 = arima.sim(n = 100, list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)),
          sd = sqrt(0.1796))
plot(t1)

# mildly long-tailed
t2 = arima.sim(n = 100, list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)),
          rand.gen = function(n, ...) sqrt(0.1796) * rt(n, df = 5))
plot(t2)

# Standard time series analysis 

# https://a-little-book-of-r-for-time-series.readthedocs.io/en/latest/src/timeseries.html

# See also https://cran.r-project.org/web/views/TimeSeries.html


# Forecasting

# https://cran.r-project.org/web/packages/tsfknn/vignettes/tsfknn.html
# https://rpubs.com/mattBrown88/TimeSeriesMachineLearning

# Nonlinearity analysis

# http://math.furman.edu/~dcs/courses/math47/R/library/tseries/html/terasvirta.test.html
# https://www.rdocumentation.org/packages/fNonlinear/versions/3042.79/topics/NonLinTests


# s-map
# I think everything Christian has explained in the last meeting is here 
# https://cran.r-project.org/web/packages/rEDM/vignettes/rEDM-tutorial.html
# install.packages("rEDM")
# There is another package with the s-map, not sure if we need this
# https://github.com/James-Thorson/NLTS
# 
