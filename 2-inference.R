## Inputs:
# samp  - a matrix with each raw being a Bb path satisfying the same settings
# N1    - an index between 1 and ncol(samp)
# T     - the expiration date (it can be ommited if t is provided)
# t     - the partition where the Bb paths was sampled (it can be omitted if T is provided)
## Output:
# Sigma - a vector which i-th entry is the MLE for the volatility of the Bb path
#         given at the i-th row of samp and using just its first N1 observations

SigMl <- function(samp, T = 1, t = NULL, N1 = 2) {
  
  N <- ncol(samp) - 1
  
  # Avoid the ending point with null variance.
  stopifnot(all(N1 < N + 1))
  
  # Check if the partition is given and otherwise create it 
  if (is.null(t)){
    
    t <- seq(0, T, by = T / N)
    
  }  else {
    
    T <- t[N + 1]
    
  }
  
  n <- nrow(samp)
  if (n == 1) samp <- as.matrix(samp)
  y <- samp[1, N + 1]
  
  mu <- (t(samp[, 1:(N1 - 1), drop = F]) * (T - t[2:N1]) +  
      y * (t[2:N1] - t[1:(N1 - 1)])) / (T - t[1:(N1 - 1)])
  
  sig <- sqrt( (t[2:N1] - t[1:(N1 - 1)]) * (T - t[2:N1]) / (T - t[1:(N1 - 1)]) )
  
  Sigma <- sqrt(colSums( ( (t(samp[, 2:N1, drop = F]) - mu) / sig ) ^ 2 ) / (N1 - 1))
  
  return( Sigma )   
  
}

###### Test ######
# a <- 0        # strike price
# y <- 0        # expiration price
# T <- 1        # expiration date
# sigma <- 1    # volatility
# N <- 2e2      # number of points (N + 1) used to discretize the interval [0, T]

# set.seed(43)
# n <- 5

# samp <- rBB(n = n, a = a, b = y, N = N, T = T, sigma = 1)
# Sigma <- SigMl(samp = samp, T = T)

# require(latex2exp)
# matplot(1:(N - 1), t(Sigma), type = "l", lty = 1, lwd = 1, col = "black", xlab = TeX("Sample size k"), ylab = "Estimation")
# lines(1:(N - 1),rep(sigma, N-1), lty = 2, lwd = 3, col = "black")
# legend("topright", legend = c(TeX("$\\widehat{\\sigma}_{k, N}$"), TeX("$\\sigma = 1$")), lty = c(1,2), lwd = c(1,3), bty = "n")

## Computes the confidence function (at level alpha) for the boundary...
## ...in case sigma was estimated via maximum likelihood.
bBBConf <- function(bnd, tol = 1e-4, y = 0, sigma = 1, T = 1, N = 1e2, 
                    t = NULL, nsamp, alpha = 0.05, eps = 1e-2,
                    oracle = FALSE, Cpp = TRUE) {
  
  # Check if the discretization is given
  if (is.null(t)) {
    
    t <- seq(0, T, by = T / N)
    
  }
  
  # Uses the true boundary (oraculo = TRUE) or the bBB algorithm (oraculo = FALSE) 
  # to compute the derivative of the boundary with respect to sigma
  if (oracle == TRUE) {
    
    # Boundary's derivative using the true boundary
    dbnd <- matrix (rep (0.8399 * sqrt(t[length(t)] - t), each = nrow(bnd)), 
                    nrow = nrow(bnd))
  
  } else {
    
    if (Cpp == TRUE){
      
      # Boundary using the incremented sigma ratio via the Cpp function bBBCpp
      bnd_eps <- bBBCpp(tol = tol, y = y, sigma = sigma + eps, t = t)
      
    } else {
      
      # Boundary using the incremented sigma via the R function bBB
      bnd_eps <- bBB(tol = tol, y = y, sigma = sigma + eps, t = t)
    
    }
    
    # Boundary derivate via the incremental ratio
    dbnd <- (bnd_eps - bnd) / eps
    
  }
  
  nalph <- length(alpha)
  nsigma <- length(sigma)
  # asymptotic lower confidence function
  bnd.low <- qnorm(rep(alpha / 2, each = (N + 1) * nsigma), mean = bnd, sd = abs(dbnd) * sigma / sqrt(2*nsamp))
  dim(bnd.low) <- c(nsigma, N + 1, nalph) 
  # asymptotic upper confidence function
  bnd.up <- qnorm(rep(1 - alpha / 2, each = (N + 1) * nsigma), mean = bnd, sd = abs(dbnd) * sigma / sqrt(2*nsamp)) 
  dim(bnd.up) <- c(nsigma, N + 1, nalph)
    
  return(list(bBB = bnd, bBB.low = drop(bnd.low), bBB.up = drop(bnd.up)))
    
}
###### Test ######
# set.seed(88)
# n <- 1
# alpha = 0.05
# t <- seq(0, T, by = T / N)

# samp <- rBB(n = n, N = N, a = a, b = y, sigma = 1)
# N1 <- floor(N/3)
# Sigma <- SigMl(samp = samp, T = T, N1 = N1)
# boundary <- bBB(y = y, sigma = drop(Sigma), T = T, N = N)
# boundary <- bBBConf(bnd = boundary, y = y, sigma = drop(Sigma), T = T, N = N, nsamp = N1, alpha = alpha)

# plot(t, samp, type = "n", xlab = "", ylab = "", ylim = c(-0.5,0.75), xlim = c(0,1), bty = "n", yaxt = "n", xaxt = "n")
# lines(t[1:N1], samp[1:N1], lty = 1, lwd = 0.5)
# lines(t[N1:(N + 1)], boundary.true[N1:(N + 1)], lty = 1, lwd = 1, col = "red")
# lines(t[N1:(N + 1)], samp[N1:(N + 1)], lty = 1, lwd = 1)
# lines(t[N1:(N + 1)], boundary$bBB[N1:(N + 1)], lty = 1, lwd = 1, col = "blue")
# lines(t[N1:(N + 1)], boundary$bBB.up[N1:(N + 1)], lty = 1, lwd = 1, col = "orange")
# lines(t[N1:(N + 1)], boundary$bBB.low[N1:(N + 1)], lty = 1, lwd = 1, col = "green")
# lines(c(t[N1],t[N1]), c(-0.4, 0.85), lty = 1)
# axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"), padj = -1.5,  line = -0.5)
# axis(2, at = c(-0.4, 0, 0.4, 0.8), labels = c("-0.4", "0", "0.4", "0.8"), padj = 1)
# legend(x = 0.001, y = 0.88, 
#        legend = c(TeX("$b_{\\sigma}$"), 
#                   TeX("$\\tilde{b}_{\\widehat{\\sigma}_{k,N}}$"), 
#                   TeX("$\\tilde{c}_{1, \\widehat{\\sigma}_{k, N}}$"), 
#                   TeX("$\\tilde{c}_{2, \\widehat{\\sigma}_{k, N}}$")), 
#        lty = c(1, 1, 1, 1), lwd = c(1, 1 ,1, 1), 
#        col = c("red", "blue", "orange", "green"), bty = "n")
