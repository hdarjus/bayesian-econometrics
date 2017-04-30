library("metRology")

# a) and b)
mu.true <- 2
sigma.true <- 1.5
alpha <- 0.05

#exc.true <- rep(0, 10000)
exc.meas <- rep(0, 10000)

for (j in 1:10000) {
  series <- rnorm(111, mu.true, sigma.true)
  
  count <- 0
  for (i in seq(10, 110)) {
    mu.hat <- mean(series[1:i])
    sigma2.hat <- var(series[1:i])
    var.hat <- qnorm(1-alpha, mu.hat, sqrt(sigma2.hat), lower.tail = T)
    count <- count + (series[i+1] > var.hat)
    
    #cat(i, ":\n mu:    \t", mu.hat, "\n sigma2:\t", sigma2.hat, "\n var:   \t", var.hat, "\n next:  \t", series[i+1], "\n")
  }
  #exc.true[j] <- sum(series > qnorm(1-alpha, mu.true, sigma.true, lower.tail = T))
  exc.meas[j] <- count
  
  #cat("Measured exceedences: ", count, "\nReal exceedenes:     ", sum(series > qnorm(1-alpha, mu.true, sigma.true, lower.tail = T)), "\n")
}

# The measured VaR using the plug-in estimate is smaller than the real one.

hist(exc.meas)
mean(exc.meas) # ~ 5.6 > 5

# c)
alpha <- 0.05
# priors
m <- 2
M <- 0.92  # just to avoid using 1, it would be boring
a <- 2
b <- 0.42

exc.meas <- rep(0, 1000)

for (j in 1:1000) {
  sigma2.true <- 1/rgamma(1, a, scale = b)
  mu.true <- rnorm(1, m, sqrt(M*sigma2.true))
  series <- rnorm(111, mu.true, sqrt(sigma2.true))
  
  count <- 0
  for (i in seq(10, 110)) {
    m.prime <- (m + M*i*mean(series[1:i]))/(M*i + 1)
    M.prime <- M/(M*i + 1)
    a.prime <- a + i/2
    b.prime <- b + 1/2 * (i-1) * var(series[1:i]) + (i * (mean(series[1:i]) - m)^2) / (2*(M*i + 1))
    var.hat <- qt.scaled(1-alpha, df = 2*a.prime, mean = m.prime, sd = sqrt(b.prime * (M.prime + 1) / a.prime), lower.tail = T)
    count <- count + (series[i+1] > var.hat)
    
    #cat(i, "\n VaR:   \t", var.hat, "\n next:  \t", series[i+1], "\n")
  }
  exc.meas[j] <- count
  
  #cat("Measured exceedences: ", count, "\nReal exceedenes:     ", sum(series > qnorm(1-alpha, mu.true, sigma.true, lower.tail = T)), "\n")
}

hist(exc.meas)
mean(exc.meas)

# d) and e)
K <- 100000
alpha <- 0.05
# prior constants
m <- 2
M <- 0.92
a <- 2
b <- 0.42

mu.true <- 2
sigma.true <- 1.5
series <- rnorm(111, mu.true, sqrt(sigma2.true))

var.values <- matrix(rep(0, K), K, 3)
  
count <- 0
for (i in seq(10, 110)) {
  # draw K times from the posterior, slide 22
  mean.y <- mean(series[1:i])
  sigma2.draws <- 1/rgamma(K, a+i/2, rate=b+.5*(sum((series[1:i]-mean.y)^2)+i/(M*i+1)*(mean.y-m)^2))
  mu.draws <- rnorm(K, (m+M*i*mean.y)/(M*i+1), sqrt(M*sigma2.draws/(M*i+1)))
  # now sample y
  y.post.pred <- rnorm(K, mu.draws, sqrt(sigma2.draws))
  # VaR
  var.hat <- quantile(y.post.pred, probs = 1-alpha)  # empirical VaR for d)
  count <- count + (series[i+1] > var.hat)
  
  if (i %in% c(10, 60, 110)) {  # for e)
    var.values[, (i+40)/50] <- qnorm(1-alpha, mu.draws, sqrt(sigma2.draws), lower.tail = T)
  }
  
  #cat(i, "\n VaR:   \t", var.hat, "\n next:  \t", series[i+1], "\n")
}

# credible intervals for e)
quantile(var.values[, 1], probs = c(2.5, 97.5)/100)
quantile(var.values[, 2], probs = c(2.5, 97.5)/100)
quantile(var.values[, 3], probs = c(2.5, 97.5)/100)
