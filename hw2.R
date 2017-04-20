B <- Inf
#B <- 5  # try commenting and uncommenting this
n <- 20000

x.vec <- rep(0,n)
y.vec <- rep(0,n)

x.vec[1] <- 0.01
y.vec[1] <- 0.1

for (i in seq(2, n)) {
  repeat {
    x <- rexp(1, y.vec[i-1])
    if (x < B)
      break
  }
  x.vec[i] <- x
  repeat {
    y <- rexp(1, x)
    if (y < B)
      break
  }
  y.vec[i] <- y
}

opar <- par(mfcol = c(3, 1), mar = c(2,3,1,1))
plot(log(x.vec), log(y.vec), type = "l")
plot(log(x.vec), type = "l")
plot(log(y.vec), type = "l")
par(opar)
