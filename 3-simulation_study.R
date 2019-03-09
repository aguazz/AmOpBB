## FirstGreater (Auxiliary Function)
# Given two matrices A and B of size n x m, this function returns 
# a vector of length n with element i equal to
# min(A[i, j] : A[i, j] >= B[i, j], j = 1,...,m), i = 1,...,n.
FirstLower <- function(A, B) {
  
  if (is.vector(drop(B))){
    
    return( A[cbind(1:nrow(A), apply( A , 1, function(x) min(which(x <= B)) )) ] )
    
  }
  
  return( A[cbind(1:nrow(A), apply( A <= B, 1, function(x) min(which(x)) )) ] )
  
}

## Test
# (A <- rbind(c(4, 6, 4, 16), c(3, 2, 21, 27)))
# (B <- rbind(c(3, 5, 7, 12), c(1, 13, 17, 23)))
# FirstLower(A, B)

## Data generator ##
# This function generates the data needed to compare the three 
# possible stopping strategies (estimated boundary and confidence functions)
# against each other and the one associated to the true boundary 
# b(t) = S - B * sigma * sqrt(T - t).
# INPUTS:
# * tol: tolerance for the fixed point algorithm inside the function b_BB_AmPut
# * sigma: volatility of the BB's paths
# * r: discount rate (NOTE: THe true boundary is only available for r = 0)
# confidence functions setup: alpha, eps
# * n: Number of BBs paths
# * a: Initial condition X_0 = a
# * S: Final condition X_T = S
# * q: quantiles levels of the marginal distributions of a BB with sigma volatility
# * T, N, t: Profits will be computed at the discretization of [0, T]: 
#     (1) using N + 1 points corresponding to the times Delta * i, 
#         with Delta = T / N and i = 0, 1, ...., N.
#     (2) or simply acording to t if it is provided.
# * t.boundary: partition where the OSB is going to be originially 
#   (before the spline interpolation) computed over
# * eps: increment length for the function b_BB_AmPu_Conf
# * ratio (natural number): determines how many observations of a BB path are 
#   between two consecutive points of t.
# For more details see Subsection 4.4 of the paper 
# "Optimal exercise of American options under stock pinning".
# OUTPUT:
# four arrays V.true, V.est, V.low and V.up,
# of dimensions n, length(t), length(q).
# each array stands for one stopping strategy, and its (i, j, k) entry
# indicates the i-th observation of the payoff of the 
# corresponding stopping strategy, associeted to initial conditions (t(j), q(k))
# and having estimated the volatility with the path from t(0) to t(j - 1)
Profit <- function(q = c(0.2, 0.4, 0.6, 0.8), tol = 1e-3, n = 1e3, 
                     a = 10, S = 10, sigma = 1, r = 0, T = 1, N = 2e2,
                     t = NULL, ratio = 1, t.boundary = NULL, alpha = 0.05, eps = 1e-2) {
  
  # Check if the discretization is given
  if (is.null(t)) {
    
    Delta <- rep(T / N, N + 1)          # Length step matrix
    t <- seq(0, T, by = T / N)          # Discretization
    t.path <- seq(0, T, by = T / (ratio * N))
    N.path <- length(Delta)
    
  } else {
    
    N <- length(t) - 1                  # Number of subintervals
    T <- t[N + 1]                       # Expiration date
    Delta <- t[2:(N + 1)] - t[1:N]      # Length step matrix
    
    t.path <- t[1]
    for(x in 1:N){
      
      v <- 1:ratio / ratio 
      t.path <- c(t.path, Delta[x] * v + t[x])
      
    }
    
    N.path <- length(Delta)
    
  }
  
  N.path <- ratio * N
  
  set.seed(23456)
  normals <- matrix(rnorm(n * (N.path - 1), mean = 0, sd = 1), 
                    nrow = n, ncol = N.path - 1, byrow = TRUE)
  SigmaM <- sigma^2 * outer(t.path, t.path, pmin)
  
  K <- length(q)
  
  # Preallocating arrays
  true.profit <- array(dim = c(n, N, K))
  est.profit <- array(dim = c(n, N, K))
  up.profit <- array(dim = c(n, N, K))
  low.profit <- array(dim = c(n, N, K))
  
  # true boundary
  bBB.true <- S - 0.8399 * sigma * sqrt(T - t.path)
  
  for (j in 2:N) {
    
    print(paste("j =", as.character(j)))
    
    for (k in 1:K) {
      
      print(paste("k =", as.character(k)))
      
      # Computing the q[k] percentile of the process at time t[j]
      Xq <- qnorm( q[k], mean = (a + (S - a) * t[j]) / T, 
                    sd = sigma * sqrt( t[j] * (T - t[j]) / T ) )
      
      # Generating a sample of n Brownian bridges forced to step by (t[j], Xq)
      samp = rBBB_chol(n = n, a = a, b = Xq, c = S, sigma = SigmaM, T = T, 
                  N1 = ratio * (j - 1) + 1, N = N.path, normals = normals)
      
      # Estimating the volatility from the Brownian bridges
      Sigma <- drop( SigMl(samp = samp, T = T, N1 = ratio * (j - 1) + 1) )
      
      # Taking only the elements of t.boundary that satisfies
      # t.boundary[i] > t[j], along with its greater element that satisfies
      # t.boundary[i] <= t[j]. This step improve dramatically the speed of the function
      # compared with taking the hole t.boundary for all j
      t.boundary.efective <- t.boundary[c(max(which(t.boundary <= t[j])), which(t.boundary > t[j]))] 
      
      # Computing the boundaries at t.boundary.efective
      boundary <- b_BB_AmPut(tol = tol, S = S, sigma = Sigma, t = t.boundary.efective, r = 0)
      
      # Computing the confidence curves
      boundary <- b_BB_AmPut_Conf(boundary, tol = tol, S = S, sigma = Sigma, r = 0, 
                                  t = t.boundary.efective, nsamp = ratio * (j - 1) + 1, 
                                  alpha = alpha, eps = eps)
      
      # Allocating temporary matrices
      temp1 <- matrix(nrow = n, ncol =  N.path - ratio * (j - 1) + 1) 
      temp2 <- temp1
      temp3 <- temp1
      
      # valuing the boundaries and the confidence curves at t[j:(N + 1)] via splinefun
      for (l in 1:n){
        
        temp1[l, ] <- splinefun(t.boundary.efective, boundary$bBB[l, ])(t.path[(ratio * (j - 1) + 1):(N.path + 1)])
        temp2[l, ] <- splinefun(t.boundary.efective, boundary$bBB.low[l, ])(t.path[(ratio * (j - 1) + 1):(N.path + 1)])
        temp3[l, ] <- splinefun(t.boundary.efective, boundary$bBB.up[l, ])(t.path[(ratio * (j - 1) + 1):(N.path + 1)])
      
      }
      
      # Updating the boundaries and the confidences curves
      boundary$bBB <- temp1
      boundary$bBB.low <- temp2
      boundary$bBB.up <- temp3
      
      # Computing the payoff associated to each stopping strategy
      true.profit[, j - 1, k] <- FirstLower(samp[, (ratio * (j - 1) + 1):(N.path + 1)], bBB.true[(ratio * (j - 1) + 1):(ratio * N + 1)])
      est.profit[, j - 1, k] <- FirstLower(samp[, (ratio * (j - 1) + 1):(N.path + 1)], boundary$bBB)
      up.profit[, j - 1, k] <- FirstLower(samp[, (ratio * (j - 1) + 1):(N.path + 1)], boundary$bBB.up)
      low.profit[, j - 1, k] <- FirstLower(samp[, (ratio * (j - 1) + 1):(N.path + 1)], boundary$bBB.low)
      
    }
    
    # Csating into matrices
    true.profit[, N, ] <- matrix(rep(S, n * K), nrow = n)
    est.profit[, N, ] <- matrix(rep(S, n * K), nrow = n)
    up.profit[, N, ] <- matrix(rep(S, n * K), nrow = n)
    low.profit[, N, ] <- matrix(rep(S, n * K), nrow = n)
    
  }
  
  return(list(true = S - true.profit, est = S - est.profit, 
              up = S - up.profit, low = S - low.profit))

}

## Generating data
source(file = "0-simulations.R")
source(file = "1-boundary_computation_BB-AmPut_discount.R")
source(file = "2-inference_BB.R")
N <- 200
T <- 1
ratio <- 1

t <- seq(0, T, by = T / N)
t.boundary <- log(seq(exp(0), exp(1), l = N + 1)) 

system.time(payoff_N200_ratio1 <- Profit(sigma = 1, t = t, t.boundary = t.boundary, 
                                          ratio = ratio))

save(payoff_N200_ratio1, file = "payoff_N200_ratio1.RData")

