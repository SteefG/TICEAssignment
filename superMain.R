rm(list = ls())

# Load Packages 
library(bootUR)
library(rugarch)
library(stargazer)

set.seed(123)	

# Function take checks if two arrays have same elements
SameElements <- function(a, b) return(identical(sort(a), sort(b)))


# Performs a Block Bootstrap Sample, returns a n length vector
block_sampler <- function(x, n, block_size) {
  ### Takes parameters, x - series, n - bootstrap sample size and block length block_size and outputs a vector Bootstrap samples of length n
  
  number_blocks <- ceiling(n/block_size) #Number of blocks in sample
  block_bootsrap_sample <- list()
  start_indices <- floor(runif(number_blocks, min = 1, max = length(x)-block_size+2)) 
  
  for(i in 1:number_blocks){
    start <- start_indices[i]
    length_list <- length(block_bootsrap_sample)
    if (length_list + block_size > n){ #Cuts off the last block
      block_bootsrap_sample <- append(block_bootsrap_sample, x[start:(start+n-length_list-1)])
      return(unlist(block_bootsrap_sample))
    }
    
    else if (length_list == n){ #As a fail-safe
      return(unlist(block_bootsrap_sample))
    }
    
    else { #Adds block to list
      test <- x[start:(start+block_size-1)]
      block_bootsrap_sample <- append(block_bootsrap_sample, x[start:(start+block_size-1)])
    }
  }
  
  return(unlist(block_bootsrap_sample))
}


##### GARCH(1,1) simulation #####
GARCH11 <- function(w, a, b, n, starting_sigma2) {
  ### Takes as input parameters, w, a, b and n signifying the parameters of a GARCH(1,1) and simulates a GARCH series 
  ### of length n.
  
  sigma2 <- rep(0, n) #Initialises the sigma^2 vector
  sigma2[1] <- starting_sigma2 #Randomly set the first value to 1
  
  epsilon <- rnorm(n, mean=0, sd = 1) #Generates white noise process with unit variance
  
  y <- rep(0, n) ##Initialises y vector
  y[1] <- sqrt(sigma2[1])*epsilon[1]
  
  for (i in 2:n) {
    sigma2[i] <- w + a * y[i-1]^2 + b * sigma2[i-1]
    y[i] <- sqrt(sigma2[i])*epsilon[i]
  }
  
  garch11 <- cbind(y, sigma2)
  colnames(garch11) <- c("y","sigma^2")
  
  return(garch11)
}

#estimated conditional variances
estimate_conditional_variances <- function(y, theta){
  ### Takes as input, series y that is a GARCH(1,1) and outputs a vector of estimated conditional variances
  
  sigmaHat2 <- rep(0,length(y))
  sigmaHat2[1] <- theta[1]/(1-theta[2]-theta[3])  #Estimated marginal variance
  
  
  for (i in 2:length(y)){
    sigmaHat2[i] <- theta[1] + theta[2]*(y[i-1])^2 + theta[3]*(sigmaHat2[i-1])
  }
  
  return(sigmaHat2)
  
}

residual_11 <- function(y, theta){
  ### Takes as input a series y and outputs the residuals after fitting a GARCH(1,1)
  
  sigmaHat2 <- estimate_conditional_variances(y, theta)
  resids <- y[1:length(y)]/sqrt(sigmaHat2[1:length(y)])
  return(resids)
}


bootRep_11 <- function(y, block_size, theta){
  ### Takes as input the series y and block size block_size and outputs a bootstrapped ystar vector (bootstrap replicates)
  
  sigmaHat2 <- estimate_conditional_variances(y, theta)
  residuals <- residual_11(y, theta)
  
  sigmaStar2 <- rep(0, length(y)) 
  sigmaStar2[1] <- sigmaHat2[1]
  
  yStar <- rep(0, length(y))
  
  eStar <- block_sampler(residuals, length(y), block_size) #e sampled from EDF of residuals with replacement
  yStar[1] <- eStar[1]*sqrt(sigmaStar2)[1]
  
  
  for (i in 2:length(y)){
    sigmaStar2[i] <- theta[1] + theta[2]*(yStar[i-1])^2 + theta[3]*(sigmaStar2[i-1])
    yStar[i] <- eStar[i]*sqrt(sigmaStar2[i])
  }
  
  return(yStar)
}



bootForecast_11 <- function(y, forecast_length, block_size, theta){ #forecast_length is the forecasting horizon
  
  garchSpec_11 <- ugarchspec(variance.model = list(garchOrder = c(1, 1)), 
                             mean.model = list(armaOrder = c(0, 0)))
  
  
  garchFit_11 <- ugarchfit(spec = garchSpec_11, data = y)
  estCoeff_11 <- c(coef(garchFit_11)[2], coef(garchFit_11)[3], coef(garchFit_11)[4])
  residuals <- residual_11(y, theta)
  
  sigmaStar2 <- rep(0, (forecast_length+1))
  sum <- 0
  
  for (i in 0:(length(y)-2)){
    sum <- sum + (estCoeff_11[3]^i)*((y[length(y)-i-1]^2) - (estCoeff_11[1]/(1-estCoeff_11[2]-estCoeff_11[3])))
  }
  
  sigmaStar2[1] <- (estCoeff_11[1]/(1-estCoeff_11[2]-estCoeff_11[3])) + estCoeff_11[2]*sum
  yStar <- rep(0, (forecast_length+1))
  eStar <- block_sampler(residuals, (forecast_length+1), block_size) #e sampled from EDF of residuals with replacement
  yStar[1] <- y[length(y)]
  
  for (k in 2:(forecast_length+1)){
    sigmaStar2[k] <- estCoeff_11[1] + estCoeff_11[2]*(yStar[k-1])^2 + estCoeff_11[3]*(sigmaStar2[k-1])
    yStar[k] <- eStar[k]*sqrt(sigmaStar2[k])
  }
  return(c(yStar[2:(forecast_length+1)], sigmaStar2[2:(forecast_length+1)])) #First K values are for y, others are for sigma
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
make_forecast <- function(y, forecast_length, B, block_size, theta){
  
  #garchSpec_11 <- ugarchspec(variance.model = list(garchOrder = c(1, 1)), 
                             #mean.model = list(armaOrder = c(0, 0)))
  
  
  out_y <- matrix(0, nrow = B, ncol = forecast_length)
  out_sigma <- matrix(0, nrow = B, ncol = forecast_length)
  
  
  for (b in 1:B){ #Repeats bootstrap B times
    seriesStar <- bootRep_11(y, block_size, theta)
    # garchStar_11 <- ugarchfit(spec = garchSpec_11, data = seriesStar)
    # coeffStar <- c(as.numeric(coef(garchStar_11)[2]),
    #                as.numeric(coef(garchStar_11)[3]), as.numeric(coef(garchStar_11)[4]))
    forecast_variable <-  bootForecast_11(seriesStar, forecast_length, block_size, theta)
    out_sigma[b,] <- forecast_variable[(forecast_length+1):(2*forecast_length)]
    out_y[b,] <- forecast_variable[1:forecast_length]
  }
  
  # seriesStar <- replicate(B, bootRep_11(y, block_size))
  # # print(seriesStar)
  # 
  # out_sigma[1:B,] <- apply(X = seriesStar, MARGIN = c(2), FUN = bootForecast_11, forecast_length = forecast_length, block_size = block_size)[(forecast_length+1):(2*forecast_length)]
  # out_y[1:B,] <- apply(X = seriesStar, MARGIN = c(2), FUN = bootForecast_11, forecast_length = forecast_length, block_size = block_size)[1:forecast_length]
  
  # print(out_y)
  # print("BREAK")
  # print(out_sigma)
  return(cbind(out_y, out_sigma))
}

get_CI <- function(y, forecast_length, B, block_size, alpha, theta){ #gets Kth forecast CI for y
  
  forecast <- make_forecast(y, forecast_length, B, block_size, theta)
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


get_asymptotic_CI <- function(y, alpha){
  
  lower <- (mean(y)-(qnorm((1-(alpha/2)))*(sqrt(var(y))/sqrt(length(y)))))
  upper <- (mean(y)+(qnorm((1-(alpha/2)))*(sqrt(var(y))/sqrt(length(y)))))
  
  return(c(lower,upper))
  
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


#n = 1000
#nr.sim = 1
#forecast_length = 5
#block_size = 10
#B = 10
#alpha = 0.05
#i=1
#R = 10


monte_carlo_simulation <- function(n, nr.sim, forecast_length, block_size, B, alpha,R){
  coverage_y <- matrix(0, nrow=nr.sim, ncol = 2)  # Matrix to store the average coverage for the returns of each simulation, 1st column is for bootstrap, 2nd is for asymptotic
  coverage_sigma2 <- matrix(0, nrow=nr.sim, ncol = 2)  # Matrix to store the average coverage for the volatility of each simulation
  average_above_y_CI <- matrix(0, nrow=nr.sim, ncol = 2)
  average_above_sigma2_CI <- matrix(0, nrow=nr.sim, ncol = 2)
  average_below_y_CI <- matrix(0, nrow=nr.sim, ncol = 2)
  average_below_sigma2_CI <- matrix(0, nrow=nr.sim, ncol = 2)
  length_y_CI <- matrix(0, nrow=nr.sim, ncol = 2)
  length_sigma2_CI <- matrix(0, nrow=nr.sim, ncol = 2)
  
  for (i in 1:nr.sim){
    print(i)
    accept_y <- matrix(0, nrow=R, ncol = 2)  # Vector to store acceptances
    accept_sigma2 <- matrix(0, nrow=R, ncol = 2)
    above_y_CI <- matrix(0, nrow=R, ncol = 2)
    above_sigma2_CI <- matrix(0, nrow=R, ncol = 2)
    below_y_CI <- matrix(0, nrow=R, ncol = 2)
    below_sigma2_CI <- matrix(0, nrow=R, ncol = 2)
    
    ## Step 1: Simulate ##
    #simulates n obervations of a GARCH(1,1)
    sim <- GARCH11(w = 0.05, a = 0.1, b = 0.85, n , 1)
    simSeries <- sim[,1]
    simVolatility <- sim[,2]
    
    
    sim_k <- matrix(0, nrow = R, ncol = 2)  #Initalizes a vector to store R future values, first column is for returns, second for volatility
    
    #generates R futures values (at horizon n+forecast_length) of returns and or volatilities
    for (r in 1:R){
      sim_horizon <- GARCH11(w = 0.05, a = 0.1, b = 0.85, forecast_length, simVolatility[n])
      sim_k[r,1] <- sim_horizon[forecast_length,1] #Obtain the quantity of interest (n+kth oberservation)
      sim_k[r,2] <- sim_horizon[forecast_length,2]
    }
    
    
    ## Step 2: Apply ##
    ##### 1. Estimate the parameters #####
    garchSpec_11 <- ugarchspec(variance.model = list(garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0)))
    garchFit_11 <- ugarchfit(spec = garchSpec_11, data = simSeries)
    estimated_coefficients <- c(coef(garchFit_11)[2], coef(garchFit_11)[3], coef(garchFit_11)[4])
    
    
    
    boot_CI <- get_CI(simSeries, forecast_length, B, block_size, alpha, estimated_coefficients)
    # print("CI is ")
    # print(CI)
    series_CI <- boot_CI[1:2]
    volatility_CI <- boot_CI[3:4]
    asymptotic_series_CI <- get_asymptotic_CI(simSeries, alpha)
    asymptotic_volatility_CI <- get_asymptotic_CI(simVolatility, alpha)
    # print(series_CI)
    # print(volatility_CI)
    
    length_y_CI[i,1] <- series_CI[2]-series_CI[1]
    length_sigma2_CI[i,1] <- volatility_CI[2]-volatility_CI[1]
    length_y_CI[i,2] <- asymptotic_series_CI[2]-asymptotic_series_CI[1]
    length_sigma2_CI[i,2] <- asymptotic_volatility_CI[2]-asymptotic_volatility_CI[1]
    
    ## Step 3: Evaluate ##
    for (r in 1:R){
      #check for bootstrap
      if (sim_k[r,1] < series_CI[2] && sim_k[r,1] > series_CI[1]) {accept_y[r,1] <- 1}
      else if (sim_k[r,1] > series_CI[2]) {above_y_CI[r,1] <-1 }
      else if(sim_k[r,1] < series_CI[1]) {below_y_CI[r,1] <-1 }
      
      if (sim_k[r,2] < volatility_CI[2] && sim_k[r,2] > volatility_CI[1]) {accept_sigma2[r,1] <- 1}
      else if (sim_k[r,2] > volatility_CI[2]) {above_sigma2_CI[r,1] <-1 }
      else if(sim_k[r,2] < volatility_CI[1]) {below_sigma2_CI[r,1] <-1 }
      
      #check for asymptotic
      if (sim_k[r,1] < asymptotic_series_CI[2] && sim_k[r,1] > asymptotic_series_CI[1]) {accept_y[r,2] <- 1}
      else if (sim_k[r,1] > asymptotic_series_CI[2]) {above_y_CI[r,2] <-1 }
      else if(sim_k[r,1] < asymptotic_series_CI[1]) {below_y_CI[r,2] <-1 }
      
      if (sim_k[r,2] < asymptotic_volatility_CI[2] && sim_k[r,2] > asymptotic_volatility_CI[1]) {accept_sigma2[r,2] <- 1}
      else if (sim_k[r,2] > asymptotic_volatility_CI[2]) {above_sigma2_CI[r,2] <-1 }
      else if(sim_k[r,2] < asymptotic_volatility_CI[1]) {below_sigma2_CI[r,2] <-1 }
      
    }
    
    coverage_y[i,] <- c(mean(accept_y[,1]),mean(accept_y[,2]))
    coverage_sigma2[i,] <- c(mean(accept_sigma2[,1]),mean(accept_sigma2[,2]))
    average_above_y_CI[i,] <- c(mean(above_y_CI[,1]),mean(above_y_CI[,2]))
    average_above_sigma2_CI[i,] <- c(mean(above_sigma2_CI[,1]),mean(above_sigma2_CI[,2]))
    average_below_y_CI[i,] <- c(mean(below_y_CI[,1]),mean(below_y_CI[,2]))
    average_below_sigma2_CI[i,] <- c(mean(below_sigma2_CI[,1]),mean(below_sigma2_CI[,2]))
  }
  ## Step 4: Summarize ##
  coverage_probability_y <- mean(coverage_y[,1])
  sd_coverage_y <- sd(coverage_y[,1])
  average_length_y <- mean(length_y_CI[,1])
  sd_length_y <- sd(length_y_CI[,1])
  coverage_probability_sigma2 <- mean(coverage_sigma2[,1])
  sd_coverage_sigma2 <- sd(coverage_sigma2[,1])
  average_length_sigma2 <- mean(length_sigma2_CI[,1])
  sd_length_sigma2 <- sd(length_sigma2_CI[,1])
  above_y <- mean(average_above_y_CI[,1])
  below_y <- mean(average_below_y_CI[,1])
  above_sigma2 <- mean(average_above_sigma2_CI[,1])
  below_sigma2 <- mean(average_below_sigma2_CI[,1])
  
  asymptotic_coverage_probability_y <- mean(coverage_y[,2])
  asymptotic_sd_coverage_y <- sd(coverage_y[,2])
  asymptotic_average_length_y <- mean(length_y_CI[,2])
  asymptotic_sd_length_y <- sd(length_y_CI[,2])
  asymptotic_coverage_probability_sigma2 <- mean(coverage_sigma2[,2])
  asymptotic_sd_coverage_sigma2 <- sd(coverage_sigma2[,2])
  asymptotic_average_length_sigma2 <- mean(length_sigma2_CI[,2])
  asymptotic_sd_length_sigma2 <- sd(length_sigma2_CI[,2])
  asymptotic_above_y <- mean(average_above_y_CI[,2])
  asymptotic_below_y <- mean(average_below_y_CI[,2])
  asymptotic_above_sigma2 <- mean(average_above_sigma2_CI[,2])
  asymptotic_below_sigma2 <- mean(average_below_sigma2_CI[,2])
  
  
  summary_boot_y <- matrix(c(coverage_probability_y, sd_coverage_y, above_y, below_y, 
                        average_length_y, sd_length_y), nrow = 1, ncol = 6)
  colnames(summary_boot_y) <- c("Average coverage", "SD", "Av. coverage above", 
                           "Av. coverage below", "Average length", "SD")
  
  summary_boot_sigma2 <- matrix(c(coverage_probability_sigma2, sd_coverage_sigma2, 
                             above_sigma2, below_sigma2, average_length_sigma2, 
                             sd_length_sigma2), nrow = 1, ncol = 6)
  
  colnames(summary_boot_sigma2) <- c("Average coverage", "SD", "Av. coverage above", 
                                "Av. coverage below", "Average length", "SD")
  
  
  
  summary_asymptotic_y <- matrix(c(asymptotic_coverage_probability_y, 
                                   asymptotic_sd_coverage_y, asymptotic_above_y, 
                                   asymptotic_below_y, 
                                   asymptotic_average_length_y, 
                                   asymptotic_sd_length_y), nrow = 1, ncol = 6)
  colnames(summary_asymptotic_y) <- c("Average coverage", "SD", "Av. coverage above", 
                                "Av. coverage below", "Average length", "SD")
  
  summary_asymptotic_sigma2 <- matrix(c(asymptotic_coverage_probability_sigma2, 
                                        asymptotic_sd_coverage_sigma2, 
                                        asymptotic_above_sigma2, asymptotic_below_sigma2, 
                                        asymptotic_average_length_sigma2, 
                                        asymptotic_sd_length_sigma2), nrow = 1, ncol = 6)
  
  colnames(summary_asymptotic_sigma2) <- c("Average coverage", "SD", "Av. coverage above", 
                                     "Av. coverage below", "Average length", "SD")
  
  return(rbind.data.frame(summary_boot_y, summary_boot_sigma2, 
                          summary_asymptotic_y, summary_asymptotic_sigma2))
}


result <- monte_carlo_simulation(1000,1,5,10,10,0.05,100)  #two first rows are the result for the bootstrap, the last two rows for the asymptotic

#stargazer(result, type = "latex", title = "Summary table", summary = FALSE,
          #column.labels = c("Average coverage", "SD", "Av. coverage above", 
                            #"Av. coverage below", "Average length", "SD"))





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

