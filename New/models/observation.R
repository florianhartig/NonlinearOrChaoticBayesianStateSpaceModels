



observationError <- function(data, type = "lognormal", param = c(sqrt(1/25.5))){
  
  if (type == "lognormal") 
      return(rlnorm(length(data), mean = log(data) , sd = param[1]))
  
}
                             
                             
                             
                             