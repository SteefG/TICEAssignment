rm(list = ls())

# Load Packages 
library(bootUR)
library(rugarch)
library(fGarch)

set.seed(123)																								# Set the true value of the mean


source("GARCH Simulations.R")

##### Block Sampler #####

#Performs a Block Bootstrap Sample, returns a n length vector
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

#Makes the block sample as a B by n matrix
make_block_sample <- function(x, n, k, B){
    out <- matrix(0, nrow = B, ncol = n)
    for (b in 1:B){
        out[b,] <- block_sampler(x, n, k)
    }
    return(out)
}




n = 100
nr.sim = 1
forecast_length = 5
block_size = 5
B = 10
alpha = 0.05

#Function to check if the mean the kth forecast is contained in the interval
check_interval <- function(interval, forecast_k) {
    in_interval <- (as.numeric(interval$lb) <= forecast_k) & (forecast_k <= as.numeric(interval$ub))	# Here I use a small R "shortcut". The output is a so-called Boolean TRUE or FALSE. However in calculations TRUE gets treated automatically as 1, and FALSE as 0.
    return(in_interval)
}



monte_carlo_simulation <- function(n, nr.sim, forecast_length, block_size, B, alpha){
    accept_y <- rep(0, times = nr.sim)  # Vector to store acceptances
    accept_sigma2 <- rep(0, times = nr.sim)
    
    for (i in 1:nr.sim){
        ## Step 1: Simulate ##
        sim <- GARCH11(w = 0.1, a = 0.5, b = 0.1, (n+forecast_length))
        simSeries_k <- sim[n+forecast_length,1] #Obtain the quantity of interest (n+kth oberservation)
        simVolatility_k <- sim[n+forecast_length,2]
        simSeries <- head(sim[,1],-forecast_length)
        print(c(simSeries_k,simVolatility_k))
        
        
        
        ## Step 2: Apply ##
        ##### 1. Estimate the parameters #####
        garchSpec_11 <- ugarchspec(variance.model = list(garchOrder = c(1, 1)), 
                                   mean.model = list(armaOrder = c(0, 0)))
        garchFit_11 <- ugarchfit(spec = garchSpec_11, data = simSeries)
        estCoeff_11 <- c(coef(garchFit_11)[2], coef(garchFit_11)[3], 
                         coef(garchFit_11)[4])
        
        
        ##### 2. Compute the volatility & residuals #####
        estimated_conditional_variances <- estimate_conditional_variances(simSeries)
        
        resids <- residual_11(simSeries)
        
        
        ##### 3. Obtain the confidence intervals #####
        series_CI <- get_y_CI(simSeries, forecast_length, B, block_size, alpha)
        print(series_CI)
        
        volatility_CI <- get_sigma_CI(simSeries, forecast_length, B, block_size, alpha)
        print(volatility_CI)
        
        
        ## Step 3: Evaluate ##
        #check_interval(simSeries_k,series_CI)
        #check_interval(simVolatility_k,volatility_CI)
        if (simSeries_k < series_CI[2] && simSeries_k > series_CI[1]) {accept_y[i] <- 1}
        if (simVolatility_k < volatility_CI[2] && simVolatility_k > volatility_CI[1]) {accept_sigma2[i] <- 1}
    }
    
    ## Step 4: Summarize ##
    coverage_probability_y <- mean(accept_y)
    coverage_probability_sigma2 <- mean(accept_sigma2)
    
    return(c(coverage_probability_y,coverage_probability_sigma2))
}



monte_carlo_simulation(100,3,2,5,10,0.05)


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
