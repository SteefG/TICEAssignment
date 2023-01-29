#install.packages("fGarch")
#install.packages('bootUR')
#install.packages("rugarch")

rm(list = ls())

# Load Packages 
library(fGarch)
library(bootUR)
library(rugarch)

#Runs *entire* Start.R, we need to only have function definitions in this file
source("Start.R")

set.seed(123)


##### GARCH(1,1) simulation #####
GARCH11 <- function(w, a, b, n) {
  ### Takes as input parameters, w, a, b and n signifying the parameters of a GARCH(1,1) and simulates a GARCH series 
  ### of length n.
  sigma2 <- rep(0, n) #Initialises the sigma^2 vector
  sigma2[1] <- 1 #Randomly set the first value to 1
  
  epsilon <- rnorm(n, mean=0, sd = 1) #Generates white noise process with unit variance
  
  y <- rep(0, n) ##Initialises y vector
  y[1] <- sqrt(sigma2[1])*epsilon[1]
  
  for (i in 2:n) {
    sigma2[i] <- w + a * y[i-1]^2 + b * sigma2[i-1]
    y[i] <- sqrt(sigma2[i])*epsilon[i]
  }
  
  y_vec <- matrix(y, nrow = n, ncol = 1)
  sigma2_vec <- matrix(sigma2, nrow = n, ncol = 1)
  garch11 <- cbind(y_vec, sigma2_vec)
  colnames(garch11) <- c("y","sigma^2")
  
  return(garch11)
}

# Generating Simulated Series 
sim <- GARCH11(w = 0.1, a = 0.5, b = 0.1, n = 10000)
simSeries <- sim[,1]
simVolatility <- sim[,2]

# PLotting Simulated Series
#plot.ts(simSeries)
#plot.ts(simVolatility)


##### 1. Estimate the parameters #####

#estimate the parameters of a simple GARCH(1,1)
garchSpec_11 <- ugarchspec(variance.model = list(garchOrder = c(1, 1)), 
                           mean.model = list(armaOrder = c(0, 0)))


#armaOrder = c(0, 0) because the mean equation doesn't have ARMA components
garchFit_11 <- ugarchfit(spec = garchSpec_11, data = simSeries)


coef(garchFit_11)

#theta_hat vector, containing the estimated (omega, alpha, beta)
estCoeff_11 <- c(coef(garchFit_11)[2], coef(garchFit_11)[3], coef(garchFit_11)[4])


##### 2. Compute the volatility & residuals #####

#estimated conditional variances
estimate_conditional_variances <- function(y){
  ### Takes as input, series y that is a GARCH(1,1) and outputs a vector of estimated conditional variances
  
  garchFit_11 <- ugarchfit(spec = garchSpec_11, data = y)
  estCoeff_11 <- c(coef(garchFit_11)[2], coef(garchFit_11)[3], coef(garchFit_11)[4])
  
  sigmaHat2 <- rep(0,length(y))
  sigmaHat2[1] <- estCoeff_11[1]/(1-estCoeff_11[2]-estCoeff_11[3]) #marginal variance
  
  
  for (i in 2:length(y)){
    sigmaHat2[i] <- estCoeff_11[1] + estCoeff_11[2]*(y[i-1])^2 + estCoeff_11[3]*(sigmaHat2[i-1])
  }
  
  return(sigmaHat2)
}

estimated_conditional_variances <- estimate_conditional_variances(simSeries) ##### doesn't need to be run

#residuals
residual_11 <- function(y){
  ### Takes as input a series y and outputs the residuals after fitting a GARCH(1,1)
  
  sigmaHat2 <- estimate_conditional_variances(y)
  
  e <- rep(0, length(y))
  for (i in 1:length(y)){
    e[i] <- y[i]/sqrt(sigmaHat2[i])
  }
  resids <- y[1:length(y)]/sqrt(sigmaHat2[1:length(y)])
  # print(resids)
  # print("This is e")
  # print(e)
  return(e)
}

resids <- residual_11(simSeries)

##### 3. Obtain bootstrap replicates #####

bootRep_11 <- function(y, block_size){
  ### Takes as input the series y and block size block_size and outputs a bootstrapped ystar vector
  
  garchFit_11 <- ugarchfit(spec = garchSpec_11, data = y)
  estCoeff_11 <- c(coef(garchFit_11)[2], coef(garchFit_11)[3], coef(garchFit_11)[4])
  ##### Why are those in the function? 
  sigmaHat2 <- estimate_conditional_variances(y)
  residuals <- residual_11(y)
  
  sigmaStar2 <- rep(0, length(y)) ##### Change variable name cuz a bit confusing
  sigmaStar2[1] <- sigmaHat2[1]
  
  yStar <- rep(0, length(y))
  
  eStar <- block_sampler(residuals, length(y), block_size) #e sampled from EDF of residuals with replacement
  yStar[1] <- eStar[1]*sqrt(sigmaStar2)[1]

  
  for (i in 2:length(y)){
    sigmaStar2[i] <- estCoeff_11[1] + estCoeff_11[2]*(yStar[i-1])^2 + estCoeff_11[3]*(sigmaStar2[i-1])
    yStar[i] <- eStar[i]*sqrt(sigmaStar2)[i]
  }

  return(yStar)
}

seriesTest <- bootRep_11(simSeries, 10) #block length of 10 for testing


##### 4. Estimate the parameters of the bootstrap series #####
#garchStar_11 <- ugarchfit(spec = garchSpec_11, data = seriesStar)
#coef(garchStar_11)
#coeffStar <- c(as.numeric(coef(garchStar_11)[2]), 
#               as.numeric(coef(garchStar_11)[3]), as.numeric(coef(garchStar_11)[4])) #Theta star hat



##### 5. Bootstrap forecasts of future values #####
bootForecast_11 <- function(y, forecast_length, block_size){ #forecast_length is the forecasting horizon
  
  garchFit_11 <- ugarchfit(spec = garchSpec_11, data = y)
  estCoeff_11 <- c(coef(garchFit_11)[2], coef(garchFit_11)[3], coef(garchFit_11)[4])
  residuals <- residual_11(y)
  
  sigmaStar2 <- rep(0, forecast_length)
  sum <- 0
  
  for (i in 0:(length(y)-2)){
    sum <- sum + (estCoeff_11[3]^i)*((y[length(y)-i-1]^2) - (estCoeff_11[1]/(1-estCoeff_11[2]-estCoeff_11[3])))
  }
  
  sigmaStar2[1] <- (estCoeff_11[1]/(1-estCoeff_11[2]-estCoeff_11[3])) + estCoeff_11[2]*sum
  yStar <- rep(0, forecast_length)
  eStar <- block_sampler(residuals, forecast_length, block_size) #e sampled from EDF of residuals with replacement
  yStar[1] <- eStar[1]*sqrt(sigmaStar2)[1]
  
  for (k in 2:forecast_length){
    sigmaStar2[k] <- estCoeff_11[1] + estCoeff_11[2]*(yStar[k-1])^2 + estCoeff_11[3]*(sigmaStar2[k-1])
    yStar[k] <- eStar[k]*sqrt(sigmaStar2)[k]
  }
  return(c(yStar, sigmaStar2)) #First K values are for y, others are for sigma
}


#bootForecast_11(seriesStar,coeffStar[1],coeffStar[2],coeffStar[3],5,r,10) #First forecast of y seems to be very off

#Makes bootstrap forecasts for y as a B by forecast_length matrix
make_y_forecast <- function(y, forecast_length, B, block_size){
  
  out <- matrix(0, nrow = B, ncol = forecast_length)
  garchFit_11 <- ugarchfit(spec = garchSpec_11, data = y)
  estCoeff_11 <- c(coef(garchFit_11)[2], coef(garchFit_11)[3], coef(garchFit_11)[4])
  estimated_conditional_variances<- estimate_conditional_variances(y)
  residuals <- residual_11(y)

  for (b in 1:B){
    seriesStar <- bootRep_11(y, block_size)
    garchStar_11 <- ugarchfit(spec = garchSpec_11, data = seriesStar)
    coeffStar <- c(as.numeric(coef(garchStar_11)[2]), 
                   as.numeric(coef(garchStar_11)[3]), as.numeric(coef(garchStar_11)[4]))
    out[b,] <- bootForecast_11(seriesStar, forecast_length, block_size)[1:forecast_length]
  }
  return(out)
}

#Makes bootstrap forecasts for sigma as a B by forecast_length matrix
make_sigma_forecast <- function(y, forecast_length, B, block_size){
  
  
  out <- matrix(0, nrow = B, ncol = forecast_length)
  #The estCoeff are constants we could put the calculation in this 
  #function depends on what you want
  
  estimated_conditional_variances<- estimate_conditional_variances(y)
  
  for (b in 1:B){ #Repeats bootstrap B times
    seriesStar <- bootRep_11(simSeries, block_size)
    garchStar_11 <- ugarchfit(spec = garchSpec_11, data = seriesStar)
    coeffStar <- c(as.numeric(coef(garchStar_11)[2]), 
                   as.numeric(coef(garchStar_11)[3]), as.numeric(coef(garchStar_11)[4]))
    out[b,] <- bootForecast_11(seriesStar, forecast_length, block_size)[(forecast_length+1):(2*forecast_length)]
  }
  return(out)
}

#Definitely have to play around with the block length. 100 seems to be ok but still not great. 
#Sometimes may fail to converge not sure why
#yK <- make_y_forecast(simSeries, 5, 100, 100)
#sigma2K <- make_sigma_forecast(simSeries, 5, 100, 100)



get_y_CI <- function(y, forecast_length, B, block_size, alpha){ #gets Kth forecast CI for y
  
  forecast <- make_y_forecast(y, forecast_length, B, block_size)
  

  forecast[,forecast_length] = sort(forecast[,forecast_length])
  
  lower <- quantile(forecast[,forecast_length], (alpha/2))
  upper <- quantile(forecast[,forecast_length], (1-(alpha/2)))

  return(list(lb = lower, ub = upper))
}



get_sigma_CI <- function(y, forecast_length, B, block_size, alpha){ #gets Kth forecast CI for simga^2

  forecast <- make_sigma_forecast(y, forecast_length, B, block_size)
  

  forecast[,forecast_length] = sort(forecast[,forecast_length])
  
  lower <- quantile(forecast[,forecast_length], (alpha/2))
  upper <- quantile(forecast[,forecast_length], (1-(alpha/2)))

  return(c(lb=lower, ub=upper))
}

get_y_CI(simSeries, 5, 499, 10, 0.05)
get_sigma_CI(simSeries, 5, 10, 10, 0.05)


#Comparing asymptotically
get_asymptotic_CI <- function(y, alpha){
  
  lower <- (mean(y)-(qnorm((1-(alpha/2)))*(sqrt(var(y))/sqrt(length(y)))))
  upper <- (mean(y)+(qnorm((1-(alpha/2)))*(sqrt(var(y))/sqrt(length(y)))))
  
  return(list(lb=lower, ub=upper))
  
}

get_asymptotic_CI(simSeries,0.05)
get_asymptotic_CI(simVolatility,0.05)



