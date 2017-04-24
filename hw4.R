# I only implement the Metropolis samplers, the Gibbs sampler is basically implemented in hw2.R

metropolis <- function (n, B, rate, delta, init = min(B/2, rate), proposal = c("uniform", "normal")) {
  if (missing(proposal) || !(proposal %in% c("uniform", "normal")))
    stop("parameter 'proposal' should be either 'uniform' or 'normal'")
  if (init >= B || init <= 0)
    stop("init has to be positive and smaller than B")
  if (proposal == "uniform") {
    #f.prop.log <- function (x, theta, delta) dunif(x, min = theta-delta, max = theta+delta, log = T)
    gen.prop <- function (n, theta, delta) runif(n, min = theta-delta, max = theta+delta)
  } else {
    #f.prop.log <- function (x, theta, delta) dnorm(x, mean = theta, sd = delta, log = T)
    gen.prop <- function (n, theta, delta) rnorm(n, mean = theta, sd = delta)
  }
  f.true.log <- function (x, B, rate) ifelse(x > 0 && x < B, dexp(x, rate = rate, log = T), -Inf)
  
  x.vec <- rep(0, n)
  x.vec[1] <- init
  for (i in seq(2, n)) {
    x.prev <- x.vec[i-1]
    x.new <- gen.prop(1, x.prev, delta)
    if (log(runif(1)) < f.true.log(x.new, B, rate) - f.true.log(x.prev, B, rate)) {  # accept
      x.vec[i] <- x.new
    } else {
      x.vec[i] <- x.prev
    }
  }
  x.vec
}

# autocorrelation as a function of delta
n <- 40000
B <- 12
delta.n <- 25
delta.values <- exp(seq(log(B/100), log(10*B), length.out = delta.n))
metro.ac.norm1 <- rep(0, delta.n)
metro.ac.unif1 <- rep(0, delta.n)
metro.ac.norm5 <- rep(0, delta.n)
metro.ac.unif5 <- rep(0, delta.n)
metro.accept.norm <- rep(0, delta.n)
metro.accept.unif <- rep(0, delta.n)
for (i in seq_len(delta.n)) {
  metro.norm <- metropolis(n, B, .1, delta.values[i], proposal = "normal")
  metro.unif <- metropolis(n, B, .1, delta.values[i], proposal = "uniform")
  metro.ac.norm1[i] <- acf(metro.norm, plot = F)[[1]][2]  # extract the lag 1 autocorrelation
  metro.ac.unif1[i] <- acf(metro.unif, plot = F)[[1]][2]
  metro.ac.norm5[i] <- acf(metro.norm, plot = F)[[1]][6]  # extract the lag 5 autocorrelation
  metro.ac.unif5[i] <- acf(metro.unif, plot = F)[[1]][6]
  metro.accept.norm[i] <- 1 - (sum(metro.norm[-1] == metro.norm[-n]) / (n-1))  # accepting but staying has probability 0
  metro.accept.unif[i] <- 1 - (sum(metro.unif[-1] == metro.unif[-n]) / (n-1))
}

opar <- par(mfcol = c(2, 1), mar = c(4, 4, 2, 1))
plot(delta.values, metro.ac.norm1, type = "l",
     ylim = c(0, 1), col = "dark red",
     xlab = "", ylab = "autocorrelation", main = "Autocorrelation and acceptance as a function of delta",
     log = "x")
lines(delta.values, metro.ac.unif1, col = "dark blue")
lines(delta.values, metro.ac.norm5, col = "red")
lines(delta.values, metro.ac.unif5, col = "light blue")
legend("bottomleft", lty = 1, col = c("dark red", "red", "dark blue", "light blue"),
       legend = c("normal lag 1", "normal lag 5", "uniform lag 1", "uniform lag 5"))

plot(delta.values, metro.accept.norm, type = "l",
     ylim = c(0, max(metro.accept.norm, metro.accept.unif)), col = "dark red",
     xlab = "", ylab = "acceptance ratio", main = NULL,
     log = "x")
lines(delta.values, metro.accept.unif, col = "dark blue")
legend("bottomleft", lty = 1, col = c("dark red", "dark blue"),
       legend = c("normal", "uniform"))
par(opar)

# acceptance ratio is decreasing with increasing delta because with larger delta
# we look further, and there's a high chance that we propose a value that has low
# probability of acceptance.
# autocorrelation looks a bit differently: it's almost 1 at the extrema but low in the middle.
# ac is influenced by both acceptance and the distance to which we jump on average when we accept,
# acceptance is high on the left but the jumps are too small, and acceptance is too low on the right.
