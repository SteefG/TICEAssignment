#install.packages("fGarch")
#install.packages('bootUR')

library(fGarch)
library(bootUR)
rm(list = ls())

seed <- 123

#generates p ""a" coefficients and q "b" coefficients, s.t their sum < 1
GARCHcoeff <- function(p,q,upbound){
  a <- rep(0,p)
  b <- rep(0,q)
  a[1] <- runif(1, 0, upbound) #initalises the first "a" coefficient
  #upper bound so the first coefficient isn't excessively large compared to the others
  for (i in 2:p){
    a[i] <- runif(1, 0, (1-sum(a))) #generated such that the sum of the coefficient stay < 1
  }
  for (i in 1:q){
    b[i] <- runif(1, 0, (1-(sum(a))+sum(b))) 
    #generated such that the sum of the coefficients of a and b stay < 1
  }
return(list(a,b))
}

GARCHcoeff(2,2,0.8)


#generates n datapoints following a GARCH(p,q) with a list of coefficients coeff
#omega is by default
simulationGARCH <- function(coeff,n){
  spec = garchSpec(model = coeff)
  data <- garchSim(spec, n)
  plot(data, type = "l", col = "blue", main = "Generated GARCH")
  return(data)
}

#example with the generated coefficients
sim_series <- simulationGARCH(GARCHcoeff(2,2,0.8),1000)
#no trend, looks stationary
stationarity <- function(series){
  adftest <- adf(series, deterministics = "intercept")
  if(adftest$p.value < 0.05){
    print("The series is stationary. Hurray!")
  }
  else{
    ("The series is not stationary.")
  }
}

stationarity(sim_series)

#example with own coefficents
#simulationGARCH(list(alpha = c(0.2, 0.05), beta = c(0.09,0.02)),20)

