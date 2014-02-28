################################################################################
# Code for reply to Perretti et al
#
# Author: Florian Hartig, http://florianhartig.wordpress.com
#
# This scripts writes the JAGS model specifications to the working directory
#
# Run after data creation
#
################################################################################

setwd("C:/temp")
#setwd("D:/home/Projekte/Papers_in_Progress/13-ResponseToPerrettiEtAl/temp")

# original model

sink("model1.txt")

cat(
  "model{
    
    r ~ dunif(2,4.5)
    K ~ dunif(0.01,10)
  
    N0 ~ dunif(0.2,1.5)
    pop[1] <- N0
    
    ProcessPrecision <- 1/0.005/0.005
    ObservationPrecision ~ dunif(0,40)
    obsSD <- sqrt( 1 / ObservationPrecision)
    
    observed[1] ~ dlnorm(log(pop[1]), ObservationPrecision)
    
    for(Time in 1:(nObs-1))
    {   

      logpop[Time+1] <- log(max(0, r * pop[Time] * (1 - pop[Time] / K) ))

      pop[Time+1] ~ dlnorm(logpop[Time+1], ProcessPrecision) 
      
      observed[Time+1] ~ dlnorm(log(pop[Time+1]), ObservationPrecision)
      
    }
  }
  ",fill = TRUE)

sink()



# Modifield model with "resetting" of the state space to avoid biasing
# estimations by wrong state estimates due to the chaotic behavior
# behavior of the model


# Save BUGS description of the model to working directory
sink("model2.txt")

cat(
  "model{
  
    # Priors
    
    r ~ dunif(2,4.5)
    K ~ dunif(0.01,10)
    
    ProcessPrecision <- 1/0.005/0.005
    ObservationPrecision ~ dunif(0.1,40)
    obsSD <- sqrt( 1 / ObservationPrecision)

    # Likelihood
    
    for(i in 1:steps)
    {
      pop[1,i] ~ dunif(0.2,1.5)
      observed[1,i] ~ dlnorm(log(pop[1,i]), ObservationPrecision)
      
      for(j in 1:(substeps-1))
      {   
        procerr[j,i] ~ dlnorm(0, ProcessPrecision)
        
        pop[j+1,i] <- max(0, r * pop[j,i] * (1 - pop[j,i] / K) ) * procerr[j,i]
        
        observed[j+1,i] ~ dlnorm(log(pop[1+j,i]), ObservationPrecision)
      }
    }
  }
  ",fill = TRUE)

sink()
