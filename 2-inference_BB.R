#### Parallel MLE for the volatility of Brownian bridges ###
## Inputs:
# * same: matrix whose rows are Brownian bridge's paths, not necessarily 
#   with the same volatility
# * N1: an index between 1 and ncol(samp), such that the first N1 observations
#   will feed the MLE algorithm
# * T: finite horizon
# * t: discretization of the interval [0, T] 
#   where the Brownian bridges were sampled
## Output:
# A vector whose i-th entry is the volatility estimated for the 
# Brownian bridge at the i-th row of samp by using its first N1 observations
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
  
  if (is.vector(samp)) samp <- matrix(samp, nrow = 1)
  y <- samp[1, N + 1]
  
  mu <- (t(samp[, 1:(N1 - 1), drop = F]) * (T - t[2:N1]) +  
      y * (t[2:N1] - t[1:(N1 - 1)])) / (T - t[1:(N1 - 1)])
  
  sig <- sqrt( (t[2:N1] - t[1:(N1 - 1)]) * (T - t[2:N1]) / (T - t[1:(N1 - 1)]) )
  
  Sigma <- sqrt(colSums( ( (t(samp[, 2:N1, drop = F]) - mu) / sig ) ^ 2 ) / (N1 - 1))
  
  return( Sigma )   
  
}

###### Test ######
# T <- 1        # expiration date
# sigma <- 1    # volatility
# N <- 2e2      # number of points (N + 1) used to discretize the interval [0, T]

# set.seed(43)
# n <- 5

# samp <- rBB(n = n, a = a, b = y, N = N, T = T, sigma = sigma)
# Sigma <- mapply(function(i) SigMl(samp = samp, T = T, N1 = i), 2:N)

# require(latex2exp)
# matplot(1:(N - 1), t(Sigma), type = "l", lty = 1, lwd = 1, col = "black", 
#         xlab = TeX("Sample size n"), ylab = "Estimation")
# lines(1:(N - 1),rep(sigma, N-1), lty = 2, lwd = 3, col = "black")
# legend("topright", legend = c(TeX("$\\widehat{\\sigma}_{n}$"), 
#                               TeX("$\\sigma = 1$")),
#        lty = c(1,2), lwd = c(1,3), bty = "n")


#### Pointwise confidence curves for the OSBs ###
## Inputs:
# * bnd: Numerical approximation of OSBs formated as the output 
#   of the function b_BB_AmPut
# * tol, S, sigma, r, T, N, t: settings used to compute bnd
# * nsamp: sample sized used to estimated the volatility of the underlying BBs
# * alpha: 100*(1-alpha) is the confidence level
# * eps: increment used to compute the difference quotient 
#   (bnd(sigma + eps) - bnd(sigma)) / eps explained in Subsection 4.3 from the
#   paper "Optimal exercise of American options under stock pinning".
## Output: A list with 3 elements:
#           * bBB (= bnd): boudnary approximation vía Algorithm 1
#           * bBB.low: lower confidence curve
#           * bBB.up: Upper confidence curve
b_BB_AmPut_Conf <- function(bnd, tol = 1e-3, S = 0, sigma = 1, r = 0.5, 
                            T = 1, N = 1e2, t = NULL, nsamp, 
                            alpha = 0.05, eps = 1e-2) {
  
  # Check if the discretization is given
  if (is.null(t)) {
    
    t <- seq(0, T, by = T / N)
    
  } else {
    
    N <- length(t) - 1
    T <- t [N + 1]
    
  }
  
  # Boundary using the incremented sigma ratio via the Cpp function bBBCpp
  bnd_eps <- b_BB_AmPut(tol = tol, S = S, sigma = sigma + eps, t = t, r = r)
      
  # Boundary derivate via the incremental ratio
  dbnd <- (bnd_eps - bnd) / eps
    
  nalph <- length(alpha)
  nsigma <- length(sigma)
  
  # asymptotic lower confidence function
  bnd.low <- qnorm(rep(alpha / 2, each = (N + 1) * nsigma), 
                   mean = bnd, sd = abs(dbnd) * sigma / sqrt(2 * nsamp))
  dim(bnd.low) <- c(nsigma, N + 1, nalph) 
  
  # asymptotic upper confidence function
  bnd.up <- qnorm(rep(1 - alpha / 2, each = (N + 1) * nsigma), 
                  mean = bnd, sd = abs(dbnd) * sigma / sqrt(2 * nsamp)) 
  dim(bnd.up) <- c(nsigma, N + 1, nalph)
    
  return(list(bBB = bnd, bBB.low = drop(bnd.low), bBB.up = drop(bnd.up)))
    
}
###### Test ######
# set.seed(88)
# n <- 1
# alpha = 0.05
# N <- 200
# T <- 1
# a <- 0
# S <- 0
# t <- seq(0, T, by = T / N)

# samp <- rBB(n = n, N = N, a = 0, b = S, sigma = 1)
# N1 <- floor(N/3)
# Sigma <- SigMl(samp = samp, T = T, N1 = N1)
# boundary <- b_BB_AmPut(S = S, sigma = drop(Sigma), T = T, N = N)
# boundary <- b_BB_AmPut_Conf(bnd = boundary, S = S, 
#             sigma = drop(Sigma), T = T, N = N, nsamp = N1, alpha = alpha)
# boundary.true <- S - 0.8399 * sigma * sqrt(T - t) # TRUE boundary for r = 0
# plot(t, samp, type = "n", xlab = "", ylab = "", 
#      ylim = c(min(samp, boundary$bBB.low[N1]), 0.75), xlim = c(0,1), 
#      bty = "n", yaxt = "n", xaxt = "n")
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
#                   TeX("$\\tilde{b}_{\\widehat{\\sigma}_{n}}$"), 
#                   TeX("$\\tilde{c}_{1, \\widehat{\\sigma}_{n}}$"), 
#                   TeX("$\\tilde{c}_{2, \\widehat{\\sigma}_{n}}$")), 
#        lty = c(1, 1, 1, 1), lwd = c(1, 1 ,1, 1), 
#        col = c("red", "blue", "orange", "green"), bty = "n")
