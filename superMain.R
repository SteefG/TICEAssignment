rm(list = ls())

# Load Packages 
library(bootUR)
library(rugarch)
library(fGarch)

set.seed(123)	

# Function take checks if two arrays have same elements
SameElements <- function(a, b) return(identical(sort(a), sort(b)))


# Performs a Block Bootstrap Sample, returns a n length vector
block_sampler <- function(x, n, k) {
  ### Takes parameters, x - series, n - bootstrap sample size and block length k.
  
  n_b <- ceiling(n/k) #Number of blocks in sample
  out <- list()
  start_indices <- floor(runif(n_b, min = 1, max = length(x)-k+2)) 
  
  for(i in 1:n_b){
    start <- start_indices[i]
    length_list <- length(out)
    if (length_list + k > n){ #Cuts off the last block
      out <- append(out, x[start:(start+n-length_list-1)])
      return(unlist(out))
    }
    
    else if (length_list == n){ #As a fail-safe
      return(unlist(out))
    }
    
    else { #Adds block to list
      test <- x[start:(start+k-1)]
      out <- append(out, x[start:(start+k-1)])
    }
  }
  
  return(unlist(out))
}

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

#estimated conditional variances
estimate_conditional_variances <- function(y){
  ### Takes as input, series y that is a GARCH(1,1) and outputs a vector of estimated conditional variances
  
  garchSpec_11 <- ugarchspec(variance.model = list(garchOrder = c(1, 1)), 
                             mean.model = list(armaOrder = c(0, 0)))
  garchFit_11 <- ugarchfit(spec = garchSpec_11, data = y)
  estCoeff_11 <- c(coef(garchFit_11)[2], coef(garchFit_11)[3], coef(garchFit_11)[4])
  
  sigmaHat2 <- rep(0,length(y))
  sigmaHat2[1] <- estCoeff_11[1]/(1-estCoeff_11[2]-estCoeff_11[3]) #marginal variance
  
  
  for (i in 2:length(y)){
    sigmaHat2[i] <- estCoeff_11[1] + estCoeff_11[2]*(y[i-1])^2 + estCoeff_11[3]*(sigmaHat2[i-1])
  }
  
  return(sigmaHat2)
  
}

residual_11 <- function(y){
  ### Takes as input a series y and outputs the residuals after fitting a GARCH(1,1)
  
  sigmaHat2 <- estimate_conditional_variances(y)
  resids <- y[1:length(y)]/sqrt(sigmaHat2[1:length(y)])
  return(resids)
}

bootRep_11 <- function(y, block_size){
  ### Takes as input the series y and block size block_size and outputs a bootstrapped ystar vector
  
  
  garchSpec_11 <- ugarchspec(variance.model = list(garchOrder = c(1, 1)), 
                             mean.model = list(armaOrder = c(0, 0)))
  
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

bootForecast_11 <- function(y, forecast_length, block_size){ #forecast_length is the forecasting horizon
  
  garchSpec_11 <- ugarchspec(variance.model = list(garchOrder = c(1, 1)), 
                             mean.model = list(armaOrder = c(0, 0)))
  
  
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

#Makes bootstrap forecasts for y as a B by forecast_length matrix
# make_y_forecast <- function(y, forecast_length, B, block_size){
#   
#   out <- matrix(0, nrow = B, ncol = forecast_length)
#   garchSpec_11 <- ugarchspec(variance.model = list(garchOrder = c(1, 1)), 
#                              mean.model = list(armaOrder = c(0, 0)))
#   
#   garchFit_11 <- ugarchfit(spec = garchSpec_11, data = y)
#   estCoeff_11 <- c(coef(garchFit_11)[2], coef(garchFit_11)[3], coef(garchFit_11)[4])
#   estimated_conditional_variances<- estimate_conditional_variances(y)
#   residuals <- residual_11(y)
#   
#   for (b in 1:B){
#     seriesStar <- bootRep_11(y, block_size)
#     garchStar_11 <- ugarchfit(spec = garchSpec_11, data = seriesStar)
#     coeffStar <- c(as.numeric(coef(garchStar_11)[2]), 
#                    as.numeric(coef(garchStar_11)[3]), as.numeric(coef(garchStar_11)[4]))
#     out[b,] <- bootForecast_11(seriesStar, forecast_length, block_size)[1:forecast_length]
#   }
#   return(out)
# }

#Makes bootstrap forecasts for sigma as a B by forecast_length matrix
make_forecast <- function(y, forecast_length, B, block_size){
  
  garchSpec_11 <- ugarchspec(variance.model = list(garchOrder = c(1, 1)), 
                             mean.model = list(armaOrder = c(0, 0)))
  
  
  out_y <- matrix(0, nrow = B, ncol = forecast_length)
  out_sigma <- matrix(0, nrow = B, ncol = forecast_length)
  
  
  for (b in 1:B){ #Repeats bootstrap B times
    seriesStar <- bootRep_11(y, block_size)
    # garchStar_11 <- ugarchfit(spec = garchSpec_11, data = seriesStar)
    # coeffStar <- c(as.numeric(coef(garchStar_11)[2]),
    #                as.numeric(coef(garchStar_11)[3]), as.numeric(coef(garchStar_11)[4]))
    out_sigma[b,] <- bootForecast_11(seriesStar, forecast_length, block_size)[(forecast_length+1):(2*forecast_length)]
    out_y[b,] <- bootForecast_11(seriesStar, forecast_length, block_size)[1:forecast_length]

  }
  # seriesStar <- replicate(B, bootRep_11(y, block_size))
  
  # out_sigma[1:B,] <- apply(X = seriesStar, MARGIN = c(1), FUN = bootForecast_11, forecast_length = forecast_length, block_size = block_size)[(forecast_length+1):(2*forecast_length)]
  # out_sigma[1:B,] <- bootForecast_11(seriesStar, forecast_length, block_size)[(forecast_length+1):(2*forecast_length)]
  # out_y[1:B,] <- bootForecast_11(seriesStar, forecast_length, block_size)[1:forecast_length]
  # out_y[1:B,] <- apply(X = seriesStar, MARGIN = c(1), FUN = bootForecast_11, forecast_length = forecast_length, block_size = block_size)[1:forecast_length]
  
  # print(out_y)
  # print("BREAK")
  # print(out_sigma)
  return(cbind(out_y, out_sigma))
}

get_CI <- function(y, forecast_length, B, block_size, alpha){ #gets Kth forecast CI for y
  
  forecast <- make_forecast(y, forecast_length, B, block_size)
  # print(forecast)
  forecast_y <- forecast[,1:forecast_length]
  forecast_sigma <- forecast[,(forecast_length + 1):(2*forecast_length)]
  
  # print(forecast_y)
  forecast_y[,forecast_length] = sort(forecast_y[,forecast_length])
  forecast_sigma[,forecast_length] = sort(forecast_sigma[,forecast_length])
  
  # print("forecast y is ")
  # print(forecast_y)
  
  lower_y <- quantile(forecast_y[,forecast_length], (alpha/2))
  lower_sigma <- quantile(forecast_sigma[,forecast_length], (alpha/2))
  upper_y <- quantile(forecast_y[,forecast_length], (1-(alpha/2)))
  upper_sigma <- quantile(forecast_sigma[,forecast_length], (1-(alpha/2)))
  
  # print(lower_y)
  # print(upper_y)
  # print(lower_sigma)
  # print(upper_sigma)
  
  return(cbind(lower_y, upper_y, lower_sigma, upper_sigma))
}


# get_sigma_CI <- function(y, forecast_length, B, block_size, alpha){ #gets Kth forecast CI for simga^2
#   
#   forecast <- make_sigma_forecast(y, forecast_length, B, block_size)
#   
#   
#   forecast[,forecast_length] = sort(forecast[,forecast_length])
#   
#   lower <- quantile(forecast[,forecast_length], (alpha/2))
#   upper <- quantile(forecast[,forecast_length], (1-(alpha/2)))
#   
#   return(list(lb=lower, ub=upper))
# }

##### 3. Obtain the confidence intervals #####


n = 1000
nr.sim = 1
forecast_length = 5
block_size = 5
B = 10
alpha = 0.05
i=1

#Function to check if the mean the kth forecast is contained in the interval
check_interval <- function(interval, forecast_k) {
  in_interval <- (interval$lb <= forecast_k) & (forecast_k <= interval$ub)	# Here I use a small R "shortcut". The output is a so-called Boolean TRUE or FALSE. However in calculations TRUE gets treated automatically as 1, and FALSE as 0.
  return(in_interval)
}



monte_carlo_simulation <- function(n, nr.sim, forecast_length, block_size, B, alpha){
  
  accept_y <- rep(0, times = nr.sim)  # Vector to store acceptances
  accept_sigma2 <- rep(0, times = nr.sim)
  
  for (i in 1:nr.sim){
    print(i)
    ## Step 1: Simulate ##
    sim <- GARCH11(w = 0.05, a = 0.1, b = 0.85, (n+forecast_length))
    simSeries_k <- sim[n+forecast_length,1] #Obtain the quantity of interest (n+kth oberservation)
    simVolatility_k <- sim[n+forecast_length,2]
    simSeries <- head(sim[,1],-forecast_length)
    # print(c(simSeries_k,simVolatility_k))
    # print("This was in the MONTE")
    
    
    ## Step 2: Apply ##
    ##### 1. Estimate the parameters #####
    garchSpec_11 <- ugarchspec(variance.model = list(garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0)))
    garchFit_11 <- ugarchfit(spec = garchSpec_11, data = simSeries)
    estCoeff_11 <- c(coef(garchFit_11)[2], coef(garchFit_11)[3], coef(garchFit_11)[4])
    
    
    
    CI <- get_CI(simSeries, forecast_length, B, block_size, alpha)
    # print("CI is ")
    # print(CI)
    series_CI <- CI[1:2]
    volatility_CI <- CI[3:4]
    # print(series_CI)
    # print(volatility_CI)
    
    
    ## Step 3: Evaluate ##
    #check_interval(simSeries_k,block_CI$series_CI)
    #check_interval(simVolatility_k,block_CI$volatility_CI)
    if (simSeries_k < series_CI[2] && simSeries_k > series_CI[1]) {accept_y[i] <- 1}
    if (simVolatility_k < volatility_CI[2] && simVolatility_k > volatility_CI[1]) {accept_sigma2[i] <- 1}
  }
  
  ## Step 4: Summarize ##
  coverage_probability_y <- mean(accept_y)
  coverage_probability_sigma2 <- mean(accept_sigma2)
  
  return(c(coverage_probability_y, coverage_probability_sigma2))
}


monte_carlo_simulation(100,3,2,5,100,0.05)


get_asymptotic_CI <- function(y, alpha){
  
  lower <- (mean(y)-(qnorm((1-(alpha/2)))*(sqrt(var(y))/sqrt(length(y)))))
  upper <- (mean(y)+(qnorm((1-(alpha/2)))*(sqrt(var(y))/sqrt(length(y)))))
  
  return(list(lb=lower, ub=upper))
  
}



#Coverage plot
n.vec <- 10*1:10														# Choose the sample sizes to run the simulation over
simulation_results <- array(0, c(length(n.vec), 3))						# Initialize a matrix to store the simulation results
for (i in 1:length(n.vec)) {											# Start the loop over the sample sizes stored in "n.vec"
  n <- n.vec[i]														# Set the sample size n
  coverage <- array(0, c(N, 3))										# Initialize a matrix to store for each simulation (1,...,N) - and all 3 intervals - whether the interval contains the mean or not
  for (s in 1:N) {													# Start the loop over the Monte Carlo simulations 1,...,N
    sim <- GARCH11(w = 0.1, a = 0.5, b = 0.1, (n+forecast_length))
    simSeries_k <- sim[n+forecast_length,1] 
    simVolatility_k <- sim[n+forecast_length,2]
    simSeries <- head(sim[,1],-forecast_length)
    
    
    asymptotic_volatility_CI <- get_asymptotic_CI(simSeries,alpha)					# Construct an asymptotic (1-alpha) confidence interval for the mean
    asymptotic_volatility_CI <- get_asymptotic_CI(simVolatility,alpha)
    #iid_boot__CI <- iid_bootstrap.interval(x, alpha, B)		# Construct an equal-tailed (1-alpha) percentile-t confidence interval for the mean using the iid bootstrap
    #iid_boot_CI <- parametric.bootstrap.interval(x, alpha, B)	# Construct an equal-tailed (1-alpha) percentile-t confidence interval for the mean using the parametric bootstrap
    
    
    coverage[s,1] <- check.interval(asymptotic_volatility_CI,alpha)				# Check if the mean beta is contained in the asymptotic confidence interval (1 if yes, 0 if no)
    #coverage[s,2] <- check.interval(iid.boot.conf.int, beta)		# Check if the mean beta is contained in the iid bootstrap confidence interval (1 if yes, 0 if no)
    #coverage[s,3] <- check.interval(par.boot.conf.int, beta)		# Check if the mean beta is contained in the parametric bootstrap confidence interval (1 if yes, 0 if no)
  }
  simulation.results[i,] <- colMeans(coverage)						# Calculate the coverage probabilities for the three intervals
}


# Plot the coverage probabilities
matplot(n.vec, cbind(1 - alpha, simulation.results), type = "l", lty = 1:6, col = c(1, 3, 2, 4), lwd = 2, xlab = expression(italic(n)), ylab = "Coverage")
legend.text <- c("Desired Coverage", "Asymptotic", "Iid bootstrap", "Blocl bootstrap")
legend("bottomright", legend.text, lty = 1:6, col = c(1, 3, 2, 4))

