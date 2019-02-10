################### BB ###################
## Simulates n BBs from X_0 = a to X_T = b based on a BM with sigma volatility
rBB <- function(n, a = 0, b = 0, sigma = 1, T = 1, N = 1e2, t = NULL) {
  
  require(mvtnorm)
  
  # Check if t is given
  if (is.null(t)){
    
    t <- seq(0, T, by = T / N)
    
  } else {
    
    N <- length(t) - 1
    T <- t[N + 1]
    
  }
  
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
                 N1 = 1e2, N = 1e2, t = NULL, normals = NULL) {
  
  # Check if t is given
  if (is.null(t)){
    
    t <- seq(0, T, by = T / N)
    
  } else {
    
    N <- length(t) - 1
    T <- t[N + 1]
    
  }
  
  # Avoid degenerate cases
  stopifnot((N1 > 1) && (N1 < N + 1))
  
  # Variance structure of a BM
  if(is.matrix(sigma)){
    
    Sigma <- sigma
    
  } else {
    
    Sigma <- sigma^2 * outer(t, t, pmin)
    
  }
  
  
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
  if (is.null(normals)) {
    
    normals <- matrix(rnorm(n * (N - 1), mean = 0, sd = 1), 
                      nrow = n, ncol = N - 1, byrow = TRUE)
    
  }
  
  # Transform to a N(mu, Sigma)
  samp <- t(t(normals %*% sqrtSigma) + mu)
  
  # Add the point where the paths are forced to stop by
  if (N1 == N) {
    
    samp <- cbind(samp[, 1:(N1 - 1)], b, c)
  
  } else { 
    
    samp <- cbind(samp[, 1:(N1 - 1)], b, samp[, N1:(N - 1)], c)
  
  }
  
  
  return(samp)
  
}

### Cholesky version ###
rBBB_chol <- function(n, a = 0, b = 1, c = 0, sigma = 1, T = 1,
                      N1 = 1e2, N = 1e2, t = NULL, normals = NULL) {
  
  # Avoid degenerate cases
  stopifnot(N1 > 1 & N1 < N + 1)
  
  # Vheck if t is given
  if (is.null(t)){
    
    t <- seq(0, T, by = T / N)
    
  } else {
    
    N <- length(t) - 1
    T <- t[N + 1]
    
  }
  
  # Variance structure of a BM
  if(is.matrix(sigma)){
    
    Sigma <- sigma
    
  } else {
    
    Sigma <- sigma^2 * outer(t, t, pmin)
    
  }
  
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
  # ev <- eigen(Sigma, symmetric = TRUE)
  # R <- t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))
  Sigma[1, 1] <- .Machine$double.eps
  R <- chol(Sigma)
  
  # Sample N(0, I) and store such that, for the same seed, n + 1 and n samples 
  # have the first n samples in common
  if (is.null(normals)) {
    
    normals <- matrix(rnorm(n * (N - 1), mean = 0, sd = 1), 
                      nrow = n, ncol = N - 1, byrow = TRUE)
    
  }
  
  # Transform to a N(mu, Sigma)
  samp <- t(t(normals %*% R) + mu)
  
  # Add the point where the paths are forced to stop by
  if (N1 == N) {
    
    samp <- cbind(samp[, 1:(N1 - 1)], b, c)
    
  } else { 
    
    samp <- cbind(samp[, 1:(N1 - 1)], b, samp[, N1:(N - 1)], c)
    
  }
  
  
  return(samp)
  
}
#
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

################### GBM ###################

## Simulates n GBMs from X_0 = a to X_T, with mu dfrit and sigma volatility
rGBM <- function(n, x = 1, mu = 0, sigma = 1, T = 1, N = 1e2, t = NULL) {
  
  require(mvtnorm)
  
  # Check if t is given
  if (is.null(t)){
    
    t <- seq(0, T, by = T / N)
    
  } else {
    
    N <- length(t) - 1
    T <- t[N + 1]
    
  }
  
  normals <- rnorm(n = (n * N), mean = rep(0, n * N),
                     sd = rep(sqrt(t[2:(N + 1)] - t[1:N]), n) )
  normals <- matrix(normals, ncol = N )
  normals <- cbind(rep(0, n), normals)
  
  BM <- t(apply(normals, 1, cumsum))
  
  BM <- x * exp(sigma * t(BM) + (mu - sigma^2 / 2) * t)
  
  return(t(BM))
  
}
##### Test ######
## rGBM
# T <- 1
# N <- 2e2
# t_line <- seq(0, T, by = T / N)
# GBMs <- rGBM(n = 10, x = 10, T = T, N = N, mu = 1, sigma = 1)
# matplot(t_line, t(GBMs), type = "l", lty = 1)
## rBBB
