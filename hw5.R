# I don't exploit analytic simplifications of the formulae

f.prop.unif <- function (x, y, B) 0
gen.prop.unif <- function (y, B) runif(1, min = 0, max = B)
f.prop.exp <- function (x, y, B) dexp(x, rate = y, log = T)
gen.prop.exp <- function (y, B) rexp(1, rate = y)
f.true <- function (x, y, B) ifelse(x <= 0 || x >= B, -Inf, dexp(x, rate = y, log = T))

gibbs.sampler <- function (n, B, init, proposal = c("exponential", "uniform")) {
  if (proposal == "uniform") {
    gen.prop <- gen.prop.unif
    f.prop <- f.prop.unif
  } else if (proposal == "exponential") {
    gen.prop <- gen.prop.exp
    f.prop <- f.prop.exp
  } else {
    stop("proposal is invalid")
  }
  x.vec <- rep(0, n)
  y.vec <- rep(0, n)
  x.vec[1] <- init[1]
  y.vec[1] <- init[2]
  for (i in seq(2, n)) {
    x <- gen.prop(y.vec[i-1], B)
    if (log(runif(1)) < f.true(x, y.vec[i-1], B) - f.true(x.vec[i-1], y.vec[i-1], B) - f.prop(x, y.vec[i-1], B) + f.prop(x.vec[i-1], y.vec[i-1], B))
      x.vec[i] <- x
    else
      x.vec[i] <- x.vec[i-1]
    
    y <- gen.prop(x.vec[i], B)
    if (log(runif(1)) < f.true(y, x.vec[i], B) - f.true(y.vec[i-1], x.vec[i], B) - f.prop(y, x.vec[i], B) + f.prop(y.vec[i-1], x.vec[i], B))
      y.vec[i] <- y
    else
      y.vec[i] <- y.vec[i-1]
  }
  cbind(x = x.vec, y = y.vec)
}

n <- 10000
B <- 10
B <- 0.1

chain.exp <- gibbs.sampler(n, B, c(1,1), "exponential")
chain.unif <- gibbs.sampler(n, B, c(1,1), "uniform")

# equivalence: the qqplot is a straight line
qqplot(chain.exp[, 1], chain.unif[, 1])

# autocorrelation: exp is much better
opar <- par(mfcol = c(1, 2), mar = c(4,4,1,1))
acf(chain.exp[, 1], lag.max = 15)
acf(chain.unif[, 1], lag.max = 15)
par(opar)

# acceptance rates
cat("Acceptance rates:\n exp: ",
    round(sum(chain.exp[2:n, 1] != chain.exp[1:(n-1), 1])/(n-1), 5),
    " %\n unif: ",
    round(sum(chain.unif[2:n, 1] != chain.unif[1:(n-1), 1])/(n-1), 5), " %\n")

# For large B the exp method is better, for small B the unif method is better.
