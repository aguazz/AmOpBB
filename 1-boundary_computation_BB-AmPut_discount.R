###### bBB - boundary computation ######
## Implementation of Algorithm 1 from the paper 
## "Optimal exercise of American options under stock pinning".
## Computes the n optimal boundaries of the discounted OSPs defined by
## Brownian bridges ending at X_T = S with volatilities given by the vector sigma
## as the underlying processes, and the gain fucntion g(x) = (S - x)^+.
## * r: discount rate (r >= 0)
## * T, N, t: [0, T] can be discretized: 
##     (1) using N + 1 points corresponding to the times Delta * i, 
##         with Delta = T / N and i = 0, 1, ...., N.
##     (2) or simply acording to t if it is provided.
## * Kernel: is a function that takes two real vectors "t" and "bt" and 
##   two matrix "t_u" and "bt_u", and return a matrix which (i, j)-th element
##   is K(t[i], bt[i], t_u[i, j], bt_u[i, j]), where K is the kernel of 
##   Voletrra integral equation characterizing the OSB
## * tol: error threshold of the fix point algorithm
b_BB_AmPut <- function (tol = 1e-3, S = 0, sigma = 1, T = 1, N = 2e2, t = NULL, 
                    r = 0.5, Kernel = BB_I_K) {
  
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
  bnd[, N] <- (S * (1 - r * (T - t[N]) / 2) / 2 - 
               sigma * sqrt((T - t[N]) / (2 * pi)) * (1 + r * (T - t[N]) / 3)) / 
              ((1 / 2) - (r * (T - t[N]) / 4)) 
  
  # Boundary computation
  for ( i in (N - 1):1 ) {
    
    bnd_old <- bnd[, i + 1] 
    e <- 1                    # error in the while-loop
    t_u <- t[(i + 1):N]       # time from the (i + 1)-th element
    bt_u <- bnd[, (i + 1):N]  # boundary from the (i + 1)-th element
    
    # Fixed point algorithm
    while (e > tol) {  
      
      # Evaluate the kernel
      K <- Kernel(t[i], bnd_old, t_u, bt_u, S, T, sigma, r)
      
      # Updating the boundary
      T_N <- T - t[N]
      T_i <- T - t[i]
      c <- 1 - 0.5 * exp(-r * (T_i - T_N)) * T_N * (1 + r * T_N / 2) / T_i
      bnd[, i] <- S - rowSums( t(Delta[i:(N - 1)] * t(K)) ) - 
                      0.5 * exp(-r * (T_i - T_N)) * (S * T_N * (1 + r * T_N / 2) / T_i + 
                                                     sigma * sqrt(2 * T_N / pi) * (1 + r * T_N / 3)) 
      bnd[, i] <- bnd[, i] / c
      
      # Relative error
      e <- max( abs((bnd[, i] - bnd_old) / bnd[, i]) ) 
      # Caution: the while will keep running until all errors associated to 
      # each of the sigma's are smaller than the tolerance

      # Update
      bnd_old <- bnd[, i]
      
    }
    
  }
  
  return(bnd)
  
}

# Kernel of the Volterra integral equation using the Mean as optimality criterion
BB_I_K <- function(t, x, t_u, bt_u, S, T, sigma, r){
  
  # Number of sigmas and length of the time grid
  n <- length(sigma)
  N <- length(t_u)
  
  # Compute z
  sigma <- matrix(rep(sigma, N), ncol = N)
  x <- matrix(rep(x, N), ncol = N)
  t_u <-matrix(rep(t_u, each = n), nrow = n)
  
  mu <- (x * (T - t_u) + S * (t_u - t)) / (T - t)
  nu <- sigma * sqrt((t_u - t) * (T - t_u) / (T - t))
  z <- (bt_u - mu) / nu
  
  # Evaluate Kernel
  K <- ((S - x) / (T - t) + r * (S - mu)) * pnorm(z, mean = 0, sd = 1, lower.tail = T) + 
    (r + (T - t_u) ^ -1) * nu * dnorm(z, mean = 0, sd = 1)
  K <- exp(-r * (t_u - t)) * K
  
  return(K)
}

###### Test ######
# Settings
# T <- 1        # expiration date
# S <- 10       # strike price
# r <- 0        # discount rate
# sigma <- .01  # volatility
# N <- 200      # number of points (N + 1) used to discretize the interval [0, T]

# Efect of the partition on the boundary computation
# t <- seq(0, T, by = T / N)
# t <- log(seq(exp(0), exp(1), l = N))
# boundary <- drop(b_BB_AmPut(S = S, sigma = sigma, t = t, r = r))
# boundary.true <- S - 0.8399 * sigma * sqrt(T - t) # TRUE boundary for r = 0
# plot(t, boundary, type = "l")
# lines(t, boundary.true, lty = 1, col = "red")