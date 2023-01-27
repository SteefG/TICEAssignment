set.seed(515)													# Set the seed
nr.sim <- 2000													# Number of simulations
B <- 499														# Number of bootstrap replications
n <- 100														# Sample size
alpha <- 0.05													# Nominal level of the test
mu <- 0															# Set the true value of the mean
reject <- rep(0, times = nr.sim)								# Vector to store rejections

for (i in 1:nr.sim){											# Start the simulations
    ## Step 1: Simulate ##
    X <- rnorm(n, mean = mu)									# Draw X
    ## Step 2: Apply ##
    X.bar <- mean(X)											# Sample mean of X
    St.Dev <- sd(X)												# Standard deviation of X
    Q <- sqrt(n)*X.bar / St.Dev									# Test statistic

    Q.star <- rep(NA, times = B)								# Vector for bootstrap quantities
    for (b in 1:B) {
        J <- sample.int(n, size = n, replace = TRUE)			# Draw the indices J
        X.star <- X[J]											# Draw the bootstrap sample
        X.bar.star <- mean(X.star)								# Bootstrap sample mean
        St.Dev.star <- sd(X.star)								# Bootstrap standard deviation
        Q.star[b] <- sqrt(n)*(X.bar.star - X.bar) / St.Dev.star	# Bootstrap test statistic
    }
    cv.star <- quantile(Q.star, probs = 1-alpha)				# Bootstrap critical value

    ## Step 3: Evaluate ##
    if (Q > cv.star) {reject[i] <- 1}							# Check if the null hypothesis is rejected
}

## Step 4: Summarize ##
ERF <- mean(reject)												# Empirical rejection frequency

## Give output on screen ##
if (mu == 0) {
    print(paste("Size using bootstrap:", ERF))
} else if (mu > 0) {
    print(paste("Power using bootstrap:", ERF))
}
