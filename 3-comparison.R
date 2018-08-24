## FirstGreader (Auxiliary Function)
# Given two matrix A and B of size n x m, this function returns 
# a vector of length n with element i equal to
# min(A[i, j] : A[i, j] >= B[i, j], j = 1,...,m), i = 1,...,n.
FirstGreater <- function(A, B) {
  
  if (is.vector(drop(B))){
    
    return( A[cbind(1:nrow(A), apply( A , 1, function(x) min(which(x >= B)) )) ] )
    
  }
  
  return( A[cbind(1:nrow(A), apply( A >= B, 1, function(x) min(which(x)) )) ] )
  
}

## Test
# (A <- rbind(c(2,4,8,16), c(1,3,9,27)))
# (B <- rbind(c(3, 5, 7, 12), c(5, 13, 17, 23)))
# FirstGreater(A, B)

##
# This function generates the data needed to compare the three 
# possible stopping strategies (estimated boundary and confidence functions)
# against each other and the one associated to the true boundary 
# b(t) = B * sigma * sqrt(T - t) + y.
# INPUTS:
# boundary computation setup: tol, a, y, sigma
# confidence infunctions setup: alpha, eps
# where to compare: t, T, N, q
# sample size: n
# OUTPUT:
# four arrays V.true, V.est, V.low and V.up,
# of dimensions n, length(t), length(q).
# each array stands for one stopping strategy, and its (i, j, k) entry
# indicates the i-th observation of the payoff associated 
# to the current stopping strategy, with initial condition (t(j), q(k))
# and having estimated the volatility with the path from t(0) to t(j - 1)

VS <- function(q = c(0.2, 0.4, 0.6, 0.8), tol = 1e-3, n = 1e3, 
                     a = 0, y = 0, sigma = 1, T = 1, N = 2e2, 
                     t = NULL, t.boundary = NULL, alpha = 0.05, eps = 1e-2) {
  
  set.seed(23456)
  normals <- matrix(rnorm(n * (N - 1), mean = 0, sd = 1), 
                    nrow = n, ncol = N - 1, byrow = TRUE)
  
  if (is.null(t)) {
    
    t <- seq(0, T, by = T / N)
  
  } else {
    
    N <- length(t) - 1
    T <- t[N + 1]
    
  }
  
  K <- length(q)
  
  # Preallocating arrays
  # The second index goes only to N - 2 because it has no sense to compare
  # at the first point and at the two second last points
  true <- array(dim = c(n, N, K))
  est <- array(dim = c(n, N, K))
  up <- array(dim = c(n, N, K))
  low <- array(dim = c(n, N, K))
  
  bBB.true <- 0.8399 * sigma * sqrt(T - t) + y
  
  for (j in 2:N) {
    
    print(paste("j =",as.character(j)))
    
    for (k in 1:K) {
      
      print(paste("k =",as.character(k)))
      
      # Computing the q[k] percentile of the process at time t[j]
      Xq <- qnorm( q[k], mean = (a + (y - a) * t[j]) / T, 
                    sd = sigma * sqrt( t[j] * (T - t[j]) / T ) )
      
      # Generating a sample of n Brownian bridges forced to step by (t[j], Xq)
      samp = rBBB(n = n, a = a, b = Xq, c = y,
                  sigma = sigma, T = T, N1 = j, N = N, normals = normals)
      
      # Estimating the volatility from the Brownian bridges
      Sigma <- drop( SigMl(samp = samp, T = T, N1 = j) )
      
      # Taking only the elements of t.boundary that satisfies
      # t.boundary[i] >= t[j], along with its greater element that satisfies
      # t.boundary[i] >= t[j]. This step improve dramatically the speed of the function
      # compared with taking the hole t.boundary for all j
      t.boundary.efective <- t.boundary[c(max(which(t.boundary <= t[j])), which(t.boundary > t[j]))] 
      
      # Computing the boundaries at t.boundary.efective
      boundary <- bBBCpp(tol = tol, y = y, sigma = Sigma, t = t.boundary.efective)
      
      # Computing the confidence curves
      boundary <- bBBConf(boundary, tol = tol, y = y, sigma = Sigma, 
                          t = t.boundary.efective, nsamp = j, alpha = alpha, eps = eps)
      
      # Allocating tmporary matrices
      temp1 <- matrix(nrow = n, ncol = N + 2 - j) 
      temp2 <- temp1
      temp3 <- temp1
      
      # valuing the boundaries and the confidence curves at t[j:(N + 1)] via splinefun
      for (l in 1:n){
        
        temp1[l, ] <- splinefun(t.boundary.efective, boundary$bBB[l, ])(t[j:(N + 1)])
        temp2[l, ] <- splinefun(t.boundary.efective, boundary$bBB.low[l, ])(t[j:(N + 1)])
        temp3[l, ] <- splinefun(t.boundary.efective, boundary$bBB.up[l, ])(t[j:(N + 1)])
      
      }
      
      # Updating the boundaries and the confidences curves
      boundary$bBB <- temp1
      boundary$bBB.low <- temp2
      boundary$bBB.up <- temp3
      
      # Computing the payoff associated to each stopping strategy
      true[, j - 1, k] <- FirstGreater(samp[, j:(N + 1)], bBB.true[j:(N + 1)])
      est[, j - 1, k] <- FirstGreater(samp[, j:(N + 1)], boundary$bBB)
      up[, j - 1, k] <- FirstGreater(samp[, j:(N + 1)], boundary$bBB.up)
      low[, j - 1, k] <- FirstGreater(samp[, j:(N + 1)], boundary$bBB.low)
      
    }
    
    # Csating into matrices
    true[, N, ] <- matrix(rep(y, n * K), nrow = n)
    est[, N, ] <- matrix(rep(y, n * K), nrow = n)
    up[, N, ] <- matrix(rep(y, n * K), nrow = n)
    low[, N, ] <- matrix(rep(y, n * K), nrow = n)
    
  }
  
  return(list(true = true, est = est, up = up, low = low))

}

## Generating data
 N <- 200
 tlog <- log(seq(exp(0), exp(1), l = N + 1)) 
 system.time(payoff_alpha05_sigma1 <- VS(sigma = 1, t.boundary = tlog))

 save(payoff_alpha05_sigma1, file = "payoff_alpha05_sigma1.RData")

