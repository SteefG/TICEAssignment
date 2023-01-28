set.seed(123)																								# Set the true value of the mean

source("Start.R")
source("GARCH Simulations.R")

n = 100
nr.sim = 1
forecast_length = 5
block_size = 5
B = 10
alpha = 0.05

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
        garchFit_11 <- ugarchfit(spec = garchSpec_11, data = simSeries)
        estCoeff_11 <- c(coef(garchFit_11)[2], coef(garchFit_11)[3], coef(garchFit_11)[4])
        
        
        ##### 2. Compute the volatility & residuals #####
        estimated_conditional_variances <- estimate_conditional_variances(simSeries)
        
        resids <- residual_11(simSeries)
        
        
        ##### 3. Obtain the confidence intervals #####
        series_CI <- get_y_CI(simSeries, forecast_length, B, block_size, alpha)
        print(series_CI)
    
        volatility_CI <- get_sigma_CI(simSeries, forecast_length, B, block_size, alpha)
        print(volatility_CI)
        
        
        ## Step 3: Evaluate ##
        if (simSeries_k < series_CI[2] && simSeries_k > series_CI[1]) {accept_y[i] <- 1}
        if (simVolatility_k < volatility_CI[2] && simVolatility_k > volatility_CI[1]) {accept_sigma2[i] <- 1}
    }
    
    ## Step 4: Summarize ##
    coverage_probability_y <- mean(accept_y)
    coverage_probability_sigma2 <- mean(accept_sigma2)
    
    return(c(coverage_probability_y,coverage_probability_sigma2))
}


monte_carlo_simulation(100,1,2,5,10,0.05)
