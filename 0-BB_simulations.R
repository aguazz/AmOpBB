## Simulates n BBs from X_0 = a to X_T = b based on a BM with sigma volatility
rBB <- function(n, a = 0, b = 0, sigma = 1, T = 1, N = 1e2) {
  require(mvtnorm)
  t <- seq(0, T, l = N + 1)
  mu <- a + t * (b - a) / T
  Sigma <- sigma ^ 2 * 
    outer(t, t, function(u, v) (T - pmax(u, v)) * pmin(u, v)) / T
  rmvnorm(n = n, mean = mu, sigma = Sigma)
}
## Simulates n BBs from X_0 = a to X_t = b to X_T = c, 
## where [0, T] is discretized in N + 1 points 
## corresponding to the times Delta * i, 
## with Delta = T / N and i = 0, 1, ...., N. We encode t = Delta * N1.
rBBB <- function(n, a = 0, b = 1, c = 0, sigma = 1, T = 1,
                 N1 = 1e2, N = 1e2, t = NULL, normals) {
  
  # Avoid degenerate cases
  stopifnot(N1 > 1 & N1 < N + 1)
  
  if (is.null(t)){
    
    t <- seq(0, T, by = T / N)
    
  } else {
    
    N <- length(t) - 1
    T <- t[N + 1]
    
  }
  
  # Variance structure of a BM
  Sigma <- sigma^2 * outer(t, t, pmin)
  
  # Conditioning of the BM to pass through X_t = b and X_T = c. Since the 
  # process is Gaussian, this is the conditioning of a multivariate Gaussian
  # on a block of entries, which results in a Gaussian with mean mu and 
  # covariance Sigma. Check Section 8.1.3 of The Matrix Cookbook. Since the BM 
  # is started at a, the Gaussian has mean a.
  indCond <- c(N1, N + 1)
  Sigma1 <- Sigma[-indCond, -indCond]
  Sigma2 <- Sigma[indCond, indCond]
  Sigma12 <- Sigma[-indCond, indCond]
  Sigma12invSigma2 <- Sigma12 %*% solve(Sigma2)
  mu <- drop(a + Sigma12invSigma2 %*% c(b - a, c - a))
  Sigma <- Sigma1 - Sigma12invSigma2 %*% t(Sigma12)
  
  # Compute sqrt of Sigma for transformation of N(0, I)
  ev <- eigen(Sigma, symmetric = TRUE)
  sqrtSigma <- t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))
  
  # Sample N(0, I) and store such that, for the same seed, n + 1 and n samples 
  # have the first n samples in common
  if (missing(normals)) {
    
    normals <- matrix(rnorm(n * (N - 1), mean = 0, sd = 1), 
                      nrow = n, ncol = N - 1, byrow = TRUE)
    
  }
  
  # Transform to a N(mu, Sigma)
  samp <- t(t(normals %*% sqrtSigma) + mu)
  
  # Add the points in which we have conditioned for the trajectories
  if (N1 == N) {
    
    samp <- cbind(samp[, 1:(N1 - 1)], b, c)
  
  } else { 
    
    samp <- cbind(samp[, 1:(N1 - 1)], b, samp[, N1:(N - 1)], c)
  
  }
  
  
  return(samp)
  
}
###### Test ######
## rBB
# T <- 1
# N <- 2e2
# t_line <- seq(0, T, by = T / N)
# BBs <- rBB(n = 1e2, a = 0, b = 0, T = T, N = N)
# matplot(t_line, t(BBs), type = "l", lty = 1)
## rBBB
# n <- 10
# normals <- matrix(rnorm(n * (N - 1), mean = 0, sd = 1), 
#                   nrow = n, ncol = N - 1, byrow = TRUE)
# rBBB(n = n, a = 0, b = 1, c = 0, T = T, N1 = N1, N = N, normals = normals)
# require(manipulate)
# manipulate({
#   BBs <- rBBB(n = n, a = 0, b = 1, c = 0, T = T, N1 = N1, N = N, normals = normals)
#   matplot(t_line, t(BBs), type = "l", lty = 1)
#   }, N1 = slider(2, N - 1, initial = floor(N / 2))
# )