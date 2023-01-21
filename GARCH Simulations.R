#install.packages("fGarch")
library(fGarch)

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
  return(garchSim(spec, n))
}

#example with the generated coefficients
simulationGARCH(GARCHcoeff(2,2,0.8),20)

#example with own coefficents
simulationGARCH(list(alpha = c(0.2, 0.05), beta = c(0.09,0.02)),20)



