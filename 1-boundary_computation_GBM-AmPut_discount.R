###### bBB - boundary computation ######
## This function computes n optimal boundaries of GBMs, 
## with drifts and volatilities given in the vectors mu and sigma respectively.
## * The gain function is G(x) = (S - x)^+, where S is the strike price.  
## * [0, T] can be discretized: 
##     (1) in N + 1 points corresponding to the times Delta * i, 
##         with Delta = T / N and i = 0, 1, ...., N.
##     (2) or simply acording to t if it is provided.
## * Kernel: is a function that takes two real vectors "t" and "bt" and 
##   two matrix "t_u" and "bt_u", and return a matrix which (i, j)th element
##   is K(t[i], bt[i], t_u[i, j], bt_u[i, j]), 
##   where K is the kernel of the Volterra integral equation characterizing the boundary
b_GBM_AmPut <- function (tol = 1e-3, mu = 0, sigma = 1, 
                         T = 1, N = 2e2, t = NULL, S = 10, 
                         Kernel = K_GBM_Amput) {
  
  # How many boundaries to compute
  n <- length(sigma)
  
  # Check if the discretization is given
  if (is.null(t)) {
    
    Delta <- rep(T / N, N + 1)          # Length step matrix
    t <- seq(0, T, by = T / N)          # Discretization
    
  } else {
    
    N <- length(t) - 1                  # Number of subintervals
    T <- t[N + 1]                       # Expiration date
    Delta <- t[2:(N + 1)] - t[1:N]      # Length step matrix
    
  }
  
  # Preallocating boundary
  bnd <- matrix(rep(S, (N + 1) * n),  ncol = N + 1)
  
  # Computing boundary's second last point 
  # bnd[, N] <- sigma / 2 * sqrt(pi * (T - t[N]) / 2)  + y
  
  # Boundary computation
  for ( i in N:1 ) {
    
    bnd_old <- bnd[, i + 1] 
    e <- 1                    # error in the while-loop
    t_u <- t[(i + 1):(N + 1)]       # time from the (i + 1)-th element
    bt_u <- bnd[, (i + 1):(N + 1)]  # boundary from the (i + 1)-th element onwards
    
    # A term needed that doesn't depend on b(t_i)
    # aux <- 2 * (T - t[i]) / ( (T - t[i]) + (t[N] - t[i]) )
    # Evaluate integral int_{t_{N - 1}}^{T} sqrt((u - t_i) / (T - u)) du
    # H <- ( T - t[i] ) * ( pi / 2 - atan(sqrt((t[N] - t[i]) / (T - t[N]))) ) +
    #  ( T - t[N] ) * sqrt((t[N] - t[i]) / (T - t[N]))
    
    # Fixed point algorithm
    while (e > tol) {  
      
      # Evaluate the kernel
      K <- Kernel(t[i], bnd_old, t_u, bt_u, y, T, mu, sigma)
      
      # Updating the boundary
      I <- integrate(function(z) pnorm((log((S - z) / bnd_old) - (mu - sigma^2 / 2) * (T - t[i]))
                                       / (sigma * sqrt(T - t[i])), 
                                       mean = 0, sd = 1, lower.tail = T), 0, S, rel.tol = 1e-8)
      bnd[, i] <- S - exp(-mu * (T - t[i])) * I$value - 
                  mu *  S * rowSums( t(Delta[i:N] * t(K)) ) 
      
      # Relative error
      e <- max( abs((bnd[, i] - bnd_old) / bnd[, i]) ) 
      # Caution: the while will keep running until all errors associated to 
      # each of the sigma's are smaller than the tolerange
      
      # Update
      bnd_old <- bnd[, i]
      
    }
    
  }
  
  return(bnd)
  
}

# Kernel of the Volterra integral equation using the Mean as optimality criterion
K_GBM_Amput <- function(t, x, t_u, bt_u, y, T, mu, sigma){
  
  # Number of sigmas and length of the time grid
  n <- length(sigma)
  N <- length(t_u)
  
  # Compute z
  mu <- matrix(rep(mu, N), ncol = N)
  sigma <- matrix(rep(sigma, N), ncol = N)
  x <- matrix(rep(x, N), ncol = N)
  t_u <-matrix(rep(t_u, each = n), nrow = n)
  
  # Evaluate Kernel
  z <- (log(bt_u/ x) - (mu - sigma^2/2) * (t_u - t)) / (sigma * sqrt(t_u - t))
  K <- exp(-mu * (t_u - t)) * pnorm(z, mean = 0, sd = 1, lower.tail = T)
  
  return(K)
}

###### Test ######
# Settings
# y <- 0        # expiration price
# T <- 1        # expiration date
# S <- 10       # starting point
# mu <- 0.7     # drift
# sigma <- 0.5  # volatility
# N <- 2e2      # number of points (N + 1) used to discretize the interval [0, T]

# Efficiency
# require(microbenchmark)
# microbenchmark(bBB(y = y, sigma = sigma, T = T, N = N), times = 100)

# Efect of the partition uon the boundary computation
# t <- seq(0, T, by = T / N)
# t <- log10(seq(10^0, 10^1, l = N))
# t <- log(seq(exp(0), exp(1), l = N))
# boundary <- b_GBM_Amput(mu = mu, sigma = sigma, t = t, S = 10)

# Some BB simulations
# source(file = "0-BB_simulations.R")
# GBMs <- rGBM(n = 5, x = S, mu = mu, sigma = sigma, t = t)

# Plot zone
# if(FigGen == TRUE) pdf("boundary_GBM_mu07_sigma05.pdf")

# matplot(t, t(GBMs), type = "l", lty = 1, xlab = "Time", ylab = "Boundary", ylim = c(min(boundary, GBMs), 18), bty = "n")
# lines(t, boundary, lty = 1, lwd = 1, col = "red")
# lines(c(0, 1), c(10, 10), lty = 2)
# title(main = "GBM | G(x) = (10 - x)^+", line = 2.8, cex.main = 0.8)
# title(main = paste("Sigma: ", as.character(sigma), sep = ""),
#       line = 1.5, cex.main = 0.8)
# title(main = paste("Mu: ", as.character(mu), sep = ""),
#       line = 0.5, cex.main = 0.8)

# if (FigGen == TRUE) dev.off()

# max(abs(boundary.true - boundary))
