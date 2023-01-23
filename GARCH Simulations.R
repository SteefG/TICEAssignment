#install.packages("fGarch")
#install.packages('bootUR')
#install.packages("rugarch")

rm(list = ls())

library(fGarch)
library(bootUR)
library(rugarch)

set.seed <- 123


##### GARCH(1,1) simulation #####
GARCH11 <- function(omega, alpha, beta, n) {
  sigma2 <- rep(0, n) #Initalises the sigma^2 vector
  sigma2[1] <- 1 #Randomly set the first value to 1
  
  eps <- rnorm(n, mean=0, sd = 1) #Generates white noise process with unit variance
  
  y <- rep(0, n) ##Initalises y vector
  y[1] <- sqrt(sigma2[1])*eps[1]
  
  for (i in 2:n) {
    sigma2[i] <- omega + alpha * eps[i-1]^2 + beta * sigma2[i-1]
    y[i] <- sqrt(sigma2[i])*eps[i]
  }
  
  return(y)
}

simSeries <- GARCH11(omega = 0.01, alpha = 0.1, beta = 0.3, n = 1000)


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


