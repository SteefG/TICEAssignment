rm(list = ls())

# Load Packages 
library(rugarch)

set.seed(1)	



block_sampler <- function(x, n, block_size) {
  ### Takes parameters, x - series, n - bootstrap sample size and 
  ### block length block_size, and outputs a vector of bootstrap samples 
  ### of length n
  
  number_blocks <- ceiling(n/block_size) # Number of blocks in the sample
  block_bootsrap_sample <- list()
  start_indices <- floor(runif(number_blocks, min = 1, 
                               max = length(x)-block_size+2)) 
  
  for(i in 1:number_blocks){
    start <- start_indices[i]
    length_list <- length(block_bootsrap_sample)
    if (length_list + block_size > n){ # Cuts off the last block
      block_bootsrap_sample <- append(block_bootsrap_sample, 
                                      x[start:(start+n-length_list-1)])
      return(unlist(block_bootsrap_sample))
    }
    
    else if (length_list == n){ # As a fail-safe
      return(unlist(block_bootsrap_sample))
    }
    
    else { # Adds block to list
      test <- x[start:(start+block_size-1)]
      block_bootsrap_sample <- append(block_bootsrap_sample, 
                                      x[start:(start+block_size-1)])
    }
  }
  
  return(unlist(block_bootsrap_sample))
}



GARCH11 <- function(w, a, b, n, starting_sigma2) {
  ### Takes as input parameters, w, a and b (parameters of a simple GARCH(1,1)), 
  ### the sample size n, and starting_sigma2, the first observed volatility
  ### Returns GARCH(1,1) returns and volatilies of length n
  
  sigma2 <- rep(0, n) # Initialises the volatility vector
  sigma2[1] <- starting_sigma2
  
  epsilon <- rnorm(n, mean=0, sd = 1) 
  # Generates innovations (here assumed to be standard normal)
  
  y <- rep(0, n) # Initialises the returns vector
  y[1] <- sqrt(sigma2[1])*epsilon[1]
  
  for (i in 2:n) { # Creates volatilities and returns of a GARCH(1,1)
    sigma2[i] <- w + a * y[i-1]^2 + b * sigma2[i-1]
    y[i] <- sqrt(sigma2[i])*epsilon[i]
  }
  
  garch11 <- cbind(y, sigma2)
  colnames(garch11) <- c("y","sigma^2")
  
  return(garch11)
}



estimate_conditional_variances <- function(y, theta){
  ### Takes as input, series y that is a GARCH(1,1) and estimated coefficients 
  ### theta = (w_hat, a_hat, b_hat) and outputs a vector of estimated 
  ### conditional variances
  
  sigmaHat2 <- rep(0,length(y)) 
  # Initialises the estimated conditional variance vector
  sigmaHat2[1] <- theta[1]/(1-theta[2]-theta[3])  #Estimated marginal variance
  
  
  for (i in 2:length(y)){ 
    # Estimates the GARCH volatility using the estimated coefficients
    sigmaHat2[i] <- theta[1] + theta[2]*(y[i-1])^2 + theta[3]*(sigmaHat2[i-1]) 
  }
  
  return(sigmaHat2)
  
}



residuals <- function(y, theta){
  ### Takes as input a series y and estimated coefficients
  ### theta = (w_hat, a_hat, b_hat), and outputs the residuals after fitting a 
  ### GARCH(1,1)
  
  sigmaHat2 <- estimate_conditional_variances(y, theta) 
  # First estimates the conditional variances
  
  resids <- y[1:length(y)]/sqrt(sigmaHat2[1:length(y)]) 
  # Residuals = observed returns/estimated standard deviation
  
  return(resids)
}



boot_replicates <- function(y, block_size, theta){
  ### Takes as input a series y, block size block_size, and estimated GARCH 
  ### coefficients theta = (w_hat, a_hat, b_hat) and outputs a bootstrapped
  ### ystar vector (bootstrap replicates)
  
  sigmaHat2 <- estimate_conditional_variances(y, theta) 

  e <- residuals(y, theta)
  
  sigmaStar2 <- rep(0, length(y)) 
  # Initialises a vector of volatility bootstrap replicates
  sigmaStar2[1] <- sigmaHat2[1] 
  # Sets the first replicate equal to the first estimate volatility
  
  yStar <- rep(0, length(y)) 
  # Initialises a vector of return bootstrap replicates
  
  eStar <- block_sampler(e, length(y), block_size) 
  # e* sampled with replacement from the EDF of the centered residuals

  yStar[1] <- eStar[1]*sqrt(sigmaStar2)[1]
  
  
  for (i in 2:length(y)){ 
    # Obtains the bootstrap replicates y* of the returns by using the sampled 
    # residuals in the GARCH(1,1) model
    
    sigmaStar2[i] <- theta[1] + theta[2]*(yStar[i-1]^2) + 
      theta[3]*(sigmaStar2[i-1])
    
    yStar[i] <- eStar[i]*sqrt(sigmaStar2[i])
  }
  
  return(yStar)
}



boot_forecast <- function(y, forecast_length, block_size, theta){
  ### Takes as input a series y, a forecasting horizon forecast_length, 
  ### a block size for the block sampling bootstrap, and the vector of estimated
  ### coefficients theta = (w_hat, a_hat, b_hat), and returns a vector of length
  ### 2*forecast_length containing the forecasted returns and then the 
  ### forecasted volatilities

  garch_specifications <- ugarchspec(variance.model = list(garchOrder=c(1, 1)), 
                             mean.model = list(armaOrder = c(0, 0))) 
  # Specifies the GARCH order (1,1), and ARMA order of the mean model (0,0)
  
  garch_fit <- ugarchfit(spec = garch_specifications, data = y) 
  # Fits the simple GARCH(1,1) to the series y
  
  estimated_coef <- c(coef(garch_fit)[2], coef(garch_fit)[3], coef(garch_fit)[4])
  # Extracts the estimated coefficients w,a,b
  
  e <- residuals(y, theta) 
  
  sigmaStar2 <- rep(0, (forecast_length+1)) 
  # Initialises a vector of bootstrap forecasted volatilities, 
  # the first entry is the last "observed" bootstrap volatility replicate: 
  # sigma^2*_n
  
  sum <- 0
  for (i in 0:(length(y)-2)){
    sum <- sum + 
      (estimated_coef[3]^i)*((y[length(y)-i-1]^2) - 
         (estimated_coef[1]/(1-estimated_coef[2]-estimated_coef[3])))
  }
  
  sigmaStar2[1] <- (estimated_coef[1]/(1-estimated_coef[2]-estimated_coef[3])) + 
    estimated_coef[2]*sum
  # sigma^2*_n
  
  eStar <- block_sampler(e, (forecast_length+1), block_size)
  
  yStar <- rep(0, (forecast_length+1))
  # Initialises a vector of bootstrap forecasted returns, 
  # the first entry is the last "observed" bootstrap return replicate: y*_n
  yStar[1] <- eStar[1]*sqrt(sigmaStar2[1])
  
  
  for (k in 2:(forecast_length+1)){ 
    # Computes the bootstrapped returns and volatilities at n+1, ..., n+forecast_length
    sigmaStar2[k] <- estimated_coef[1] + 
      estimated_coef[2]*(yStar[k-1]^2) + 
      estimated_coef[3]*sigmaStar2[k-1]
    yStar[k] <- eStar[k]*sqrt(sigmaStar2[k])
  }
  
  return(c(yStar[2:length(yStar)], sigmaStar2[2:length(sigmaStar2)])) 
  # We remove the first entry as it's the forecast at horizon 0 
  # (last observation)
}



make_forecast <- function(y, forecast_length, B, block_size, theta){
  ### Takes as input a series y, a forecasting horizon forecast_length,  
  ### a number of of bootstrap replicates B, a block size for the block sampling 
  ### bootstrap, and the vector of estimated coefficients 
  ### theta = (w_hat, a_hat, b_hat), and returns a matrices containing B
  ### bootstrap forecasts (from k=1,..., forecast_length) 

  
  forecast_y <- matrix(0, nrow = B, ncol = forecast_length)
  forecast_sigma <- matrix(0, nrow = B, ncol = forecast_length)
  
  
  for (b in 1:B){
    seriesStar <- boot_replicates(y, block_size, theta)

    forecast_sigma[b,] <- boot_forecast(seriesStar, forecast_length, block_size, 
                                        theta)[(forecast_length+1):(2*forecast_length)]
    forecast_y[b,] <- boot_forecast(seriesStar, forecast_length, 
                                    block_size, theta)[1:forecast_length]
    
  }
  
  return(cbind(forecast_y, forecast_sigma))
}



get_CI <- function(y, forecast_length, B, block_size, alpha, theta){ 
  ### Also takes the significance level alpha we want to test, and 
  ### returns the confidence interval for the forcasted return and volatility 
  ### at horizon forecast_length
  
  forecast <- make_forecast(y, forecast_length, B, block_size, theta)
  forecast_y <- forecast[,1:forecast_length]
  forecast_sigma <- forecast[,(forecast_length + 1):(2*forecast_length)]
  
  forecast_y[,forecast_length] = sort(forecast_y[,forecast_length])
  forecast_sigma[,forecast_length] = sort(forecast_sigma[,forecast_length])
  # Sorts the last forecasted bootstrap values in increasing order

  
  # Computes upper and lower bounds of the intervals
  lower_y <- quantile(forecast_y[,forecast_length], (alpha/2))
  upper_y <- quantile(forecast_y[,forecast_length], (1-(alpha/2)))
  
  lower_sigma <- quantile(forecast_sigma[,forecast_length], (alpha/2))
  upper_sigma <- quantile(forecast_sigma[,forecast_length], (1-(alpha/2)))
  
  
  return(cbind(lower_y, upper_y, lower_sigma, upper_sigma))
}



get_asymptotic_CI <- function(y, alpha, theta, forecast_length, sigma2){
  ### Returns the asymptotic confidence interval for a GARCH(1,1) forecasted 
  ### return at horizon forecast_length. The method used comes from the paper.
  
  a <- theta[1]
  b <- theta[2]
  w <- theta[3]
  
  n <- length(y)
  
  expect_sigma2T1 <- (w/(1-a-b)) + 
    ((a + b)^(forecast_length) * (w+a*y[n]^2+b*sigma2[n]))
  
  lower <- qnorm(alpha/2)*expect_sigma2T1
  upper <- -1*qnorm(alpha/2)*expect_sigma2T1
  
  return(c(lower,upper))
  
}



monte_carlo <- function(n, nr.sim, forecast_length, block_size, B, alpha, R){
  ### Runs a Monte Carlo simulation for the bootstrap CI for forecasted returns 
  ### and volatilities of a GARCH(1,1), and an asymptotic CI for the return.
  ### The simulations are run nr.sim times, and the function returns the average 
  ### coverage and its standard deviation, the percentage of forecasts falling 
  ### above and bellow the CI, and the average length of the CI and its 
  ### standard deviation
  
  coverage <- list(y = matrix(0, nrow=nr.sim, ncol = 2), 
                   sigma2 = matrix(0, nrow=nr.sim, ncol = 1))
  average_above <- list(y = matrix(0, nrow=nr.sim, ncol = 2), 
                           sigma2 = matrix(0, nrow=nr.sim, ncol = 1))
  average_below <- list(y = matrix(0, nrow=nr.sim, ncol = 2), 
                           sigma2 = matrix(0, nrow=nr.sim, ncol = 1))
  length_CI <- list(y = matrix(0, nrow=nr.sim, ncol = 2), 
                    sigma2 = matrix(0, nrow=nr.sim, ncol = 1))
  # First column of a matrix for the bootstrap CI, second for the asymptotic CI
  
  
  for (i in 1:nr.sim){
    print(i)
    accept <- list(y = matrix(0, nrow=R, ncol = 2), 
                   sigma2 = matrix(0, nrow=R, ncol = 1))
    above <- list(y = matrix(0, nrow=R, ncol = 2), 
                  sigma2 = matrix(0, nrow=R, ncol = 1))
    below <- list(y = matrix(0, nrow=R, ncol = 2), 
                  sigma2 = matrix(0, nrow=R, ncol = 1))

    
    # Simulates n observations of a GARCH(1,1)
    sim <- GARCH11(w = 0.05, a = 0.1, b = 0.85, n , 1) 
    simSeries <- sim[,1]
    simVolatility <- sim[,2]
    
    
    sim_future <- matrix(0, nrow = R, ncol = 2)  
    
    
    for (r in 1:R){
      # Generates R futures values (at horizon forecast_length) of returns and 
      # volatilities
      sim_horizon <- GARCH11(w = 0.05, a = 0.1, b = 0.85, 
                             (forecast_length+1), simVolatility[n])
      # The starting volatility value is the last observed volatility 
      
      sim_future[r,1] <- sim_horizon[(forecast_length+1),1] 
      # Contains R future returns
      sim_future[r,2] <- sim_horizon[(forecast_length+1),2] 
      # Contains R future volatilities
    }
    
    
    # Estimates the GARCH parameters from the simulated data
    garch_specifications <- ugarchspec(variance.model = list(garchOrder=c(1,1)), 
                                       mean.model = list(armaOrder = c(0, 0)))
    garch_fit_11 <- ugarchfit(spec = garch_specifications, data = simSeries)
    estimated_coefficients_11 <- c(coef(garch_fit_11)[2], coef(garch_fit_11)[3], 
                                   coef(garch_fit_11)[4])
    
    
    # Obtains the bootstrap CIs
    boot_CI <- get_CI(simSeries, forecast_length, B, block_size, alpha, 
                      estimated_coefficients_11)
    series_CI <- boot_CI[1:2]
    volatility_CI <- boot_CI[3:4]
  
      
    # Obtains the asymptotic CI
    asymptotic_series_CI <- get_asymptotic_CI(simSeries, alpha, 
                                              estimated_coefficients_11, 
                                              forecast_length, simVolatility)
    
  
    #Computes the length of the CIs
    length_CI$y[i,1] <- series_CI[2]-series_CI[1]
    length_CI$sigma2[i] <- volatility_CI[2]-volatility_CI[1]
    length_CI$y[i,2] <- asymptotic_series_CI[2]-asymptotic_series_CI[1]
    
    
  
    for (r in 1:R){
      # Checks the bootstrap CI for each future value generated
      # If in the interval, accept = 1, if above the upper bound, above = 1, 
      # if below the lower bound, below = 1
      if (sim_future[r,1] < series_CI[2] && sim_future[r,1] > series_CI[1]) {
        accept$y[r,1] <- 1
        }
      else if (sim_future[r,1] > series_CI[2]) {above$y[r,1] <-1 } 
      else if(sim_future[r,1] < series_CI[1]) {below$y[r,1] <-1 }
      
      if (sim_future[r,2] < volatility_CI[2] && sim_future[r,2] > volatility_CI[1]) {
        accept$sigma2[r] <- 1
        }
      else if (sim_future[r,2] > volatility_CI[2]) {above$sigma2[r] <-1 }
      else if(sim_future[r,2] < volatility_CI[1]) {below$sigma2[r] <-1 }
      
      # Checks the asymptotic CI
      if (sim_future[r,1] < asymptotic_series_CI[2] && 
          sim_future[r,1] > asymptotic_series_CI[1]) {
        accept$y[r,2] <- 1
        }
      else if (sim_future[r,1] > asymptotic_series_CI[2]) {above$y[r,2] <-1 }
      else if(sim_future[r,1] < asymptotic_series_CI[1]) {below$y[r,2] <-1 }
      
    }
    
    
    coverage$y[i,] <- colMeans(accept$y) # Return coverage of this simulation
    # Bootstrap coverage is the first column, asymptotic the second
    coverage$sigma2[i] <- mean(accept$sigma2) # Computes the volatility coverage 
    #of this simulation
    
    average_above$y[i,] <- colMeans(above$y) # Computes the return proportion
    # above the CI in this simulation. Bootstrap first column, asymptotic second
    average_above$sigma2[i] <- mean(above$sigma2) # For volatility
    
    average_below$y[i,] <- colMeans(below$y) # Same for below
    average_below$sigma2[i] <- mean(below$sigma2)
    
  }
  
  
  # Summarises the results for future returns using the bootstrap CI
  boot_coverage_y <- mean(coverage$y[,1]) # Average of the coverages found in 
  #each simulation
  boot_sd_coverage_y <- sd(coverage$y[,1]) # and its standard deviation
  
  boot_above_y <- mean(average_above$y[,1]) # Average proportion of futures 
  # values above the CI
  boot_below_y <- mean(average_below$y[,1]) # Same for below
  
  boot_length_y <- mean(length_CI$y[,1]) # Average length of the CI
  boot_sd_length_y <- sd(length_CI$y[,1])
  
  
  summary_boot_y <- matrix(c(boot_coverage_y, boot_sd_coverage_y, boot_above_y, 
                             boot_below_y, 
                             boot_length_y, boot_sd_length_y), nrow = 1, 
                           ncol = 6)
  colnames(summary_boot_y) <- c("Average coverage", "SD", "Av. coverage above", 
                                "Av. coverage below", "Average length", "SD")
  
  
  # Summarises the results for future volatilities using the bootstrap CI
  boot_coverage_sigma2 <- mean(coverage$sigma2)
  boot_sd_coverage_sigma2 <- sd(coverage$sigma2)
  
  above_sigma2 <- mean(average_above$sigma2)
  below_sigma2 <- mean(average_below$sigma2)
  
  boot_length_sigma2 <- mean(length_CI$sigma2)
  boot_sd_length_sigma2 <- sd(length_CI$sigma2)
  
  
  summary_boot_sigma2 <- matrix(c(boot_coverage_sigma2, boot_sd_coverage_sigma2, 
                             above_sigma2, below_sigma2, 
                             boot_length_sigma2, boot_sd_length_sigma2), 
                             nrow = 1, ncol = 6)
  colnames(summary_boot_sigma2) <- c("Average coverage", "SD", 
                                     "Av. coverage above", 
                                "Av. coverage below", "Average length", "SD")
  
  
  # Summarises the results for future returns using the asymptotic CI
  asymptotic_coverage_y <- mean(coverage$y[,2]) 
  asymptotic_sd_coverage_y <- sd(coverage$y[,2]) 
  
  asymptotic_above_y <- mean(average_above$y[,2]) 
  asymptotic_below_y <- mean(average_below$y[,2])
  
  asymptotic_length_y <- mean(length_CI$y[,2])
  asymptotic_sd_length_y <- sd(length_CI$y[,2])

  
  summary_asymptotic_y <- matrix(c(asymptotic_coverage_y, 
                                   asymptotic_sd_coverage_y, asymptotic_above_y, 
                                   asymptotic_below_y, 
                                   asymptotic_length_y, 
                                   asymptotic_sd_length_y), nrow = 1, ncol = 6)
  colnames(summary_asymptotic_y) <- c("Average coverage", "SD", 
                                      "Av. coverage above", 
                                      "Av. coverage below", 
                                      "Average length", "SD")
  
  
  return(rbind.data.frame(summary_boot_y, summary_boot_sigma2, 
                          summary_asymptotic_y))
}


monte_carlo(500, 25, 50, 10, 199, 0.05, 100)

monte_carlo(500, 25, 50, 20, 199, 0.05, 100)
