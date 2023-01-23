#install.packages("fGarch")
#install.packages('bootUR')
#install.packages("rugarch")

library(fGarch)
library(bootUR)
library(rugarch)
rm(list = ls())

seed <- 123

#generates p ""a" coefficients and q "b" coefficients, s.t their sum < 1
GARCHcoeff <- function(p,q,upbound){
  a <- rep(0,p)
  b <- rep(0,q)
  a[1] <- runif(1, 0, upbound) #initalises the first "a" coefficient
  #upper bound so the first coefficient isn't excessively large compared to the others
  if (p > 1){
    for (i in 2:p){
      a[i] <- runif(1, 0, (1-sum(a))) #generated such that the sum of the coefficient stay < 1
    }
  }
  for (i in 1:q){
    b[i] <- runif(1, 0, (1-(sum(a))+sum(b))) 
    #generated such that the sum of the coefficients of a and b stay < 1
  }
return(list(a,b))
}

#GARCHcoeff(2,2,0.8)
trueCoef_11 <- GARCHcoeff(1,1,1)

#generates n datapoints following a GARCH(p,q) with a list of coefficients coeff
#omega is by default
simulationGARCH <- function(coeff,n){
  spec = garchSpec(model = coeff)
  data <- garchSim(spec, n)
  plot(data, type = "l", col = "blue", main = "Generated GARCH")
  return(data)
}

#example with the generated coefficients
#simSeries <- simulationGARCH(trueCoef_11,1000)
simSeries <- simulationGARCH(list(alpha=0.5,beta=0.1),1000)

#no trend, looks stationary
stationarity <- function(series){
  adftest <- adf(series, deterministics = "intercept")
  if(adftest$p.value < 0.05){
    print("The series is stationary.")
  }
  else{
    ("The series is not stationary.")
  }
}

stationarity(simSeries)

#example with own coefficents
#simulationGARCH(list(alpha = c(0.2, 0.05), beta = c(0.09,0.02)),20)


#find the best lag order for a GARCH model (to do later, first GARCH(1,1))

##### 1. Estimate the parameters #####
#estimate the parameters of a simple GARCH(1,1)
garchSpec_11 <- ugarchspec(variance.model = list(garchOrder = c(1, 1)), 
                           mean.model = list(armaOrder = c(0, 0)))
#armaOrder = c(0, 0) because the mean equation doesn't have ARMA components
garchFit_11 <- ugarchfit(data = simSeries, spec = garchSpec_11)
#can have problems converging (??)
coef(garchFit_11) #estimates are horrendous

#theta_hat vector, containing the estimated (omega, alpha, beta)
estCoeff_11 <- c(as.numeric(coef(garchFit_11)[2]), 
                 as.numeric(coef(garchFit_11)[3]), as.numeric(coef(garchFit_11)[4]))

##### 2. Compute the volatility & residuals #####
#estimated conditional variances
estCondVar_11 <- function(y, w, a, b){
  n = length(y)
  sigmaHat2 <- rep(0,n)
  sigmaHat2[1] <- w/(1-a-b) #marginal variance
  for (i in 2:n){
    sigmaHat2[i] <- w + a*(y[i-1])^2 + b*(sigmaHat2[i-1])^2
  }
  return(sigmaHat2)
}

s <- estCondVar_11(simSeries,estCoeff_11[1],estCoeff_11[2],estCoeff_11[3])

#residuals
residual_11 <- function(sigmaHat2,y){
  n = length(y)
  e <- rep(0,n)
  for (i in 1:n){
    e[i] <- y[i]/sqrt(sigmaHat2[i])
  }
  return(e)
}

r <- residual_11(s,simSeries)

##### 3. Obtain bootstrap replicates #####
bootRep_11 <- function(y,w,a,b,sigmaHat2,eSample){
  n = length(y)
  sigmaStar2 <- rep(0,n)
  sigmaStar2[1] <- sigmaHat2[1]
  yStar <- rep(0,n)
  eStar <- sample_edf(y,n) #e sampled from edf with replacement
  yStar[1] <- eStar[1]*sqrt(sigmaStar2)[1]
  for (i in 2:n){
    sigmaStar2[i] <- w + a*(yStar[i-1])^2 + b*(sigmaStar2[i-1])^2
    yStar[i] <- eStar[i]*sqrt(sigmaStar2)[i]
  }
  return(yStar)
}

seriesStar <- bootRep_11(simSeries,estCoeff_11[1],estCoeff_11[2],estCoeff_11[3],s,r)


##### 4. Estimate the parameters of the bootstrap series #####
garchStar_11 <- ugarchfit(data = seriesStar, spec = garchSpec_11)
coef(garchStar_11)
coeffStar <- c(as.numeric(coef(garchStar_11)[2]), 
               as.numeric(coef(garchStar_11)[3]), as.numeric(coef(garchStar_11)[4]))


##### 5. Bootstrap forecasts of future values #####
bootForecast_11 <- function(y,w,a,b,K,eSample){ #K is the forecasting horizon
  n = lenght(y)
  sigmaStar2 <- rep(0,K)
  sum = 0
  for (j in 0:(n-2)){
    sum <- sum + b^j 
  }
  sigmaStar2[1] <- (w/(1-a-b))+a
  yStar <- rep(0,K)
  eStar <- sample_edf(y,K) #e sampled from edf with replacement
  yStar[1] <- eStar[1]*sqrt(sigmaStar2)[1]
  for (k in 2:K){
    sigmaStar2[i] <- w + a*(yStar[i-1])^2 + b*(sigmaStar2[i-1])^2
    yStar[i] <- eStar[i]*sqrt(sigmaStar2)[i]
  }
  return(yStar)
}

bootForecast_11(seriesStar,coeffStar[1],coeffStar[2],coeffStar[3],5)


