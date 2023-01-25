#install.packages("fGarch")
#install.packages('bootUR')
#install.packages("rugarch")

rm(list = ls())

library(fGarch)
library(bootUR)
library(rugarch)
#Runs *entire* Start.R, we need to only have function definitions in this file
source("Start.R")

set.seed(123)


##### GARCH(1,1) simulation #####
GARCH11 <- function(omega, alpha, beta, n) {
  sigma2 <- rep(0, n) #Initialises the sigma^2 vector
  sigma2[1] <- 0.0001 #Randomly set the first value to 1
  
  eps <- rnorm(n, mean=0, sd = 1) #Generates white noise process with unit variance
  
  y <- rep(0, n) ##Initialises y vector
  y[1] <- sqrt(sigma2[1])*eps[1]
  
  for (i in 2:n) {
    sigma2[i] <- omega + alpha * y[i-1]^2 + beta * sigma2[i-1]
    y[i] <- sqrt(sigma2[i])*eps[i]
  }
  
  y_vec <- matrix(y, nrow = n, ncol = 1)
  sigma2_vec <- matrix(y, nrow = n, ncol = 1)
  garch11 <- cbind(y_vec, sigma2_vec)
  colnames(garch11) <- c("y","simga^2")
  
  return(garch11)
}


simSeries <- GARCH11(omega = 0.01, alpha = 0.1, beta = 0.3, n = 1000)[,1]
simVolatility <- GARCH11(omega = 0.01, alpha = 0.1, beta = 0.3, n = 1000)[,2]
#plot(simSeries)
plot.ts(simSeries)
plot.ts(simVolatility)

##### 1. Estimate the parameters #####
#estimate the parameters of a simple GARCH(1,1)
garchSpec_11 <- ugarchspec(variance.model = list(garchOrder = c(1, 1)), 
                           mean.model = list(armaOrder = c(0, 0)))
garchSpec_11

#armaOrder = c(0, 0) because the mean equation doesn't have ARMA components
garchFit_11 <- ugarchfit(data = simSeries, spec = garchSpec_11)
garchFit_11

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

bootRep_11 <- function(y,w,a,b,sigmaHat2,eSample, k){
  n = length(y)
  sigmaStar2 <- rep(0,n)
  sigmaStar2[1] <- sigmaHat2[1]
  yStar <- rep(0,n)
  eStar <- block_sampler(r,n,k) #e sampled from EDF of residuals with replacement
  yStar[1] <- eStar[1]*sqrt(sigmaStar2)[1]
  for (i in 2:n){
    sigmaStar2[i] <- w + a*(yStar[i-1])^2 + b*(sigmaStar2[i-1])^2
    yStar[i] <- eStar[i]*sqrt(sigmaStar2)[i]
  }
  return(yStar)
}

#seriesStar <- bootRep_11(simSeries,estCoeff_11[1],estCoeff_11[2],estCoeff_11[3],s,r, 10) #block length of 10 for testing


##### 4. Estimate the parameters of the bootstrap series #####
#garchStar_11 <- ugarchfit(data = seriesStar, spec = garchSpec_11)
#coef(garchStar_11)
#coeffStar <- c(as.numeric(coef(garchStar_11)[2]), 
#               as.numeric(coef(garchStar_11)[3]), as.numeric(coef(garchStar_11)[4])) #Theta star hat


##### 5. Bootstrap forecasts of future values #####
bootForecast_11 <- function(y,w,a,b,K,eSample,blockSize){ #K is the forecasting horizon
  n = length(y)
  sigmaStar2 <- rep(0,K)
  sigmaStar2[1] <- (w/(1-a-b)) + (a*sum((y[1:(n-2)])^2 - (w/(1-a-b))))
  yStar <- rep(0,K)
  eStar <- block_sampler(r,K,blockSize) #e sampled from EDF of residuals with replacement
  yStar[1] <- eStar[1]*sqrt(sigmaStar2)[1]
  for (k in 2:K){
    sigmaStar2[k] <- w + a*(yStar[k-1])^2 + b*(sigmaStar2[k-1])^2
    yStar[k] <- eStar[k]*sqrt(sigmaStar2)[k]
  }
  return(c(yStar, sigmaStar2)) #First K values are for y, others are for sigma
}

#bootForecast_11(seriesStar,coeffStar[1],coeffStar[2],coeffStar[3],5,r,10) #First forecast of y seems to be very off

#Makes bootstrap forecasts for y as a B by K matrix
make_y_forecast <- function(y, K, B, blockSize){
  n <- length(y)
  out <- matrix(0, nrow = B, ncol = K)
  s <- estCondVar_11(y,estCoeff_11[1],estCoeff_11[2],estCoeff_11[3])
  r <- residual_11(s,y)

  for (b in 1:B){
    seriesStar <- bootRep_11(simSeries,estCoeff_11[1],estCoeff_11[2],estCoeff_11[3],s,r, blockSize)
    garchStar_11 <- ugarchfit(data = seriesStar, spec = garchSpec_11)
    coeffStar <- c(as.numeric(coef(garchStar_11)[2]), 
                   as.numeric(coef(garchStar_11)[3]), as.numeric(coef(garchStar_11)[4]))
    out[b,] <- bootForecast_11(seriesStar,coeffStar[1],coeffStar[2],coeffStar[3],K,r,blockSize)[1:K]
  }
  return(out)
}

#Makes bootstrap forecasts for sigma as a B by K matrix
make_sigma_forecast <- function(y, K, B, blockSize){
  n <- length(y)
  out <- matrix(0, nrow = B, ncol = K)
  #The estCoeff are constants we could put the calculation in this 
  #function depends on what you want
  s <- estCondVar_11(y,estCoeff_11[1],estCoeff_11[2],estCoeff_11[3])
  r <- residual_11(s,y)
  
  for (b in 1:B){ #Repeats bootstrap B times
    seriesStar <- bootRep_11(simSeries,estCoeff_11[1],estCoeff_11[2],estCoeff_11[3],s,r, blockSize)
    garchStar_11 <- ugarchfit(data = seriesStar, spec = garchSpec_11)
    coeffStar <- c(as.numeric(coef(garchStar_11)[2]), 
                   as.numeric(coef(garchStar_11)[3]), as.numeric(coef(garchStar_11)[4]))
    out[b,] <- bootForecast_11(seriesStar,coeffStar[1],coeffStar[2],coeffStar[3],K,r,blockSize)[(K+1):(2*K)]
  }
  return(out)
}

#Definitely have to play around with the block length. 100 seems to be ok but still not great. 
#Sometimes may fail to converge not sure why
make_y_forecast(simSeries, 5, 10, 100)
make_sigma_forecast(simSeries, 5, 10, 100)
