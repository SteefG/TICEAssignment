#clear global environment
# rm(list = ls())



monte_carlo_normal <- function(N, sample_size, alpha, mu, sigma){
  count <- 0
  for (i in 1:N ) {
    rand_norm <- rnorm(n = sample_size, mean = mu, sd = sigma)
    CI_lower <- (mean(rand_norm)-(qnorm((1-(alpha/2)))*(sqrt(variance(rand_norm))/sqrt(sample_size))))
    CI_upper <- (mean(rand_norm)+(qnorm((1-(alpha/2)))*(sqrt(variance(rand_norm))/sqrt(sample_size))))
    if (mu >= CI_lower && mu <= CI_upper) {
      count <- count + 1
    }
  }
  coverage <- count/N
}
# As N increases then coverage should converge to 1-alpha by LLN
# Plot the coverage as N increases and as sample size increases

#Empirical Distribution Function implementation
edf <- function(x, h) {
  out<-numeric(length(h))
  length <- length(h) #To stop it from recalculating this in the loop
  for(i in 1:length) {
    indicator <- as.numeric(x<=h[i])
    out[i] <- sum(indicator)/length(h)
  }
  out
}


#Sampling from data with replacement, returns an n length vector
sample_edf <- function(x, n){ 
  out <- rep(0, n)
  #+1 added to upper-bound since we floor the random number
  values <- floor(runif(n, min = 1, max = length(x)+1))
  for (i in 1:n){
    out[i] <- x[values[i]]
  }
  return(out)
}

#Makes bootstrap sample as a B by n matrix
make_boot_sample <- function(x, n, B){ 
  out <- matrix(0, nrow = B, ncol = n)
  for (b in 1:B){
    out[b,] <- sample_edf(x, n)
  }
  return(out)
}

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



