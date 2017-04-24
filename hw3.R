# I only implement the two rejection samplers and not the whole Gibbs sampling procedure.
# It's basically done in hw2.R

# a)
sample.reject.exp <- function (n, B, rate) {
  x.vec <- rep(0, n)
  count <- 0
  for (i in seq_len(n)) {
    repeat {
      count <- count+1
      x <- rexp(1, rate = rate)
      if (x < B)
        break
    }
    x.vec[i] <- x
  }
  list(samples = x.vec, count = count)
}

# b)
sample.reject.unif <- function (n, B, rate) {
  x.vec <- rep(0, n)
  count <- 0
  C <- 1 / (1-exp(-B*rate))
  M <- B * rate * C
  for (i in seq_len(n)) {
    repeat {
      count <- count+1
      x <- runif(1, min = 0, max = B)
      if (runif(1) < C * rate * exp(-x*rate) / (M * (1 / B)))
        break
    }
    x.vec[i] <- x
  }
  list(samples = x.vec, count = count)
}

# visual test, the 3 histograms should look similar
# the first one has less elements than the other 2
B <- 2
rate <- .5
n <- 10000
x.plain <- rexp(n, rate)
x.exp <- sample.reject.exp(n, B, rate)$samples
x.unif <- sample.reject.unif(n, B, rate)$samples
opar <- par(mfcol = c(1,3))
hist(x.plain[x.plain < B])
hist(x.exp)
hist(x.unif)
par(opar)

# time using several parameters
n <- 20000
B <- 1

rate <- 0.1
system.time({  # slower
  sample.reject.exp(n, B, rate)
})
system.time({  # faster
  sample.reject.unif(n, B, rate)
})

rate <- 1
system.time({  # faster
  sample.reject.exp(n, B, rate)
})
system.time({  # slower
  sample.reject.unif(n, B, rate)
})


B <- 10

rate <- 0.1
system.time({  # faster
  sample.reject.exp(n, B, rate)
})
system.time({  # slower
  sample.reject.unif(n, B, rate)
})

rate <- 1
system.time({  # fastest among all
  sample.reject.exp(n, B, rate)
})
system.time({  # slowest among all
  sample.reject.unif(n, B, rate)
})

# rejections
n <- 1000
rate <- 1
B.n <- 25
B.values <- 10^seq(-1.5, .5, length.out = B.n)
rejects.exp <- rep(0, B.n)
rejects.unif <- rep(0, B.n)
for (i in seq_len(B.n)) {
  rejects.exp[i] <- sample.reject.exp(n, B.values[i], rate)$count / n - 1
  rejects.unif[i] <- sample.reject.unif(n, B.values[i], rate)$count / n - 1
}

# warning: the plot doesn't work if rejects.exp or rejects.unif have zeros
plot(B.values, rejects.exp, col = "blue", type = "l",
     xlim = c(-.5, B.values[B.n]+.5), ylim = c(0.005, max(c(rejects.exp, rejects.unif))),
     main = "Ratio of rejections over successes as the function of B on a log scale",
     xlab = "B", ylab = "#(rejects) / #(successes)",
     log = "y")
points(B.values, rejects.exp, col = "blue")
points(B.values, rejects.unif, col = "red")
lines(B.values, rejects.unif, col = "red")
legend("topright", col = c("blue", "red"), legend = c("exp", "unif"),
       lty = 1)
