hw1 <- function (m, n, seed) {
  # generate joint distribution
  f.joint <- matrix(runif(m*n), m, n)  # rows are Y values, columns are X values,
  # just like the bottom of page 169,
  # because in this paper they multiply with transition matrices from the right
  f.joint <- f.joint / sum(f.joint)  # normalize to probabilities
  #print(f.joint)
  
  # conditional distributions
  f.x.given.y <- f.joint / rowSums(f.joint)
  f.y.given.x <- t(f.joint) / rowSums(t(f.joint))
  
  # rename them to transition matrices :)
  f.x.to.y <- f.y.given.x
  f.y.to.x <- f.x.given.y
  f.x.to.x <- f.x.to.y %*% f.y.to.x
  f.y.to.y <- f.y.to.x %*% f.x.to.y
  
  # initial values as marginal distributions
  f.x <- runif(n)
  f.x <- t(f.x / sum(f.x))  # also transform to a row matrix
  f.y <- runif(m)
  f.y <- t(f.y / sum(f.y))
  
  f.x.init <- f.x
  f.y.init <- f.y
  
  # the main functionality
  max.iterations <- 100000
  converged <- FALSE
  iter <- 0
  while (iter < max.iterations && !converged) {
    f.x.previous <- f.x
    f.y.previous <- f.y
    f.y <- f.x %*% f.x.to.y
    f.x <- f.y %*% f.y.to.x
    if (sum(abs(f.x - f.x.previous)) + sum(abs(f.y - f.y.previous)) < 1e-5)
      converged <- TRUE
    iter <- iter+1
  }
  
  # print results
  cat("Number of iterations: ", iter, "\nResulting distributions:\n initial f.x:\n")
  print(f.x.init)
  cat(" previous f.x:\n")
  print(f.x.previous)
  cat(" last f.x:\n")
  print(f.x)
}
