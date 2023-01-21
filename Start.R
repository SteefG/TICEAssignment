#clear global environment
rm(list = ls())

rand_norm <- rnorm(n = 10, mean = 0, sd = 1)

hist(rand_norm, xlab = "Random value (X)", col = "grey",
     main = "", cex.lab = 1.5, cex.axis = 1.5)

qnorm(0.05)
qnorm(0.95)
qnorm(0.975)

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


edf <- function(x, increment ) {
  out<-numeric(length(increment))
  for(i in 1:length(increment)) {
    indicator <- as.numeric(x<=increment[i])
    out[i] <- sum(indicator)/length(increment)
  }
  out
}

# test to see if the edf function matches the base R one
x <- rnorm(1000)
plot(ecdf(x))
lines(seq(-4, 4, length=1000), 
      edf(x, seq(-4, 4, length=1000)), 
      col='red')


