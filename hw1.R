hw1 <- function (m, n, seed, ruin = FALSE) {
  # generate joint distribution
  f.joint <- matrix(runif(m*n), m, n)  # rows are Y values, columns are X values,
  # just like the bottom of page 169,
  # because in this paper they multiply with transition matrices from the right
  
  if (ruin) {
    f.joint[, ] <- 0  # set almost everything to zero to ruin the sampler
    f.joint[m, n] <- 1  # only almost everything :)
  }
        
  f.joint <- f.joint / sum(f.joint)  # normalize to probabilities
  #print(f.joint)
  
  # conditional distributions
  f.x.given.y <- f.joint / ifelse(rowSums(f.joint) == 0, rep(1, m), rowSums(f.joint))
  f.y.given.x <- t(f.joint) / ifelse(rowSums(t(f.joint)) == 0, rep(1, n), rowSums(t(f.joint)))
  
  # rename them to transition matrices :)
  f.x.to.y <- f.y.given.x
  f.y.to.x <- f.x.given.y
  f.x.to.x <- f.x.to.y %*% f.y.to.x
  f.y.to.y <- f.y.to.x %*% f.x.to.y
  
  # initial values as marginal distributions
  f.x <- runif(n)
  if (ruin) f.x[-1] <- 0
  f.x <- t(f.x / sum(f.x))  # also transform to a row matrix
  f.y <- runif(m)
  if (ruin) f.y[-1] <- 0
  f.y <- t(f.y / sum(f.y))
  
  f.x.init <- f.x
  f.y.init <- f.y
  
  # the main functionality
  max.iterations <- 10000
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

# General run
hw1(10, 11, 42)

# Ruined sampler (results in 0 distribution)
hw1(10, 11, 42, TRUE)
