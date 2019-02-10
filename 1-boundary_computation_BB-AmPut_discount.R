###### bBB - boundary computation ######
## This function computes n optimal boundaries of BBs 
## ending at X_T = S with volatilities given by the vector sigma.
## * discount: TRUE
## * The gain function is g(x) = x 
## * [0, T] can be discretized: 
##     (1) in N + 1 points corresponding to the times Delta * i, 
##         with Delta = T / N and i = 0, 1, ...., N.
##     (2) or simply acording to t if it is provided.
## * Kernel: is a function that takes two real vectors "t" and "bt" and 
##   two matrix "t_u" and "bt_u", and return a matrix which (i, j)th element
##   is K(t[i], bt[i], t_u[i, j], bt_u[i, j]), 
##   where K is the kernel of the Volterra integral equation characterizing the boundary
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
      # each of the sigma's are smaller than the tolerange

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
# S <- 10        # expiration price
# r <- 0        # discount rate
# sigma <- .01  # volatility
# N <- 200      # number of points (N + 1) used to discretize the interval [0, T]

# Efficiency
# require(microbenchmark)
# microbenchmark(bBB(y = y, sigma = sigma, T = T, N = N), times = 100)

# Efect of the partition on the boundary computation
# t <- seq(0, T, by = T / N)
# t <- log10(seq(10^0, 10^1, l = N))
# t <- log(seq(exp(0), exp(1), l = N))
# boundary <- drop(b_BB_AmPut(S = S, sigma = sigma, t = t, r = r))
# boundary.true <- S - 0.8399 * sigma * sqrt(T - t)
# plot(t, boundary, type = "l")
# lines(t, boundary.true, lty = 1, col = "red")
 
# Some BB simulations
# source(file = "0-BB_simulations.R")
# BBs <- rBB(n = 3, a = y, b = y, sigma = sigma, t = t)

# Covert "." into "*"
# yy <- strsplit(as.character(y), split = "")[[1]]
# if(any(yy == ".")) yy[which(yy == ".")] <- "*"
# yy <- paste(yy, collapse = "")

# rr <- strsplit(as.character(r), split = "")[[1]]
# if(any(rr == ".")) rr[which(rr == ".")] <- "*"
# rr <- paste(rr, collapse = "")

# ss <- strsplit(as.character(sigma), split = "")[[1]]
# if(any(ss == ".")) ss[which(ss == ".")] <- "*"
# ss <- paste(ss, collapse = "")

# Plot zone
# if(FigGen == TRUE) pdf(paste("boundary_BB_y", yy, "_r", rr, "_sigma", ss, ".pdf", sep = ""))

# matplot(t, t(BBs), type = "l", lty = 1.5, xlab = "Time", 
#         ylab = "Boundary", 
#         ylim = c(min(BBs, boundary), max(rbind(BBs, boundary))),
#         bty = "n", col = c("blue", "red", "green"))
# matlines(t, exp(-r * t) * t(BBs), lty = 2, lwd = 0.75, col = c("blue", "red", "green"))
# lines(t, boundary, lty = 1, lwd = 2.5, col = "red")
# lines(c(0, 1), c(y, y), lty = 2)
# title(main = "Bb | G(x) = x", line = 1.4, cex.main = 0.8)
# mtext(side = 3, line = -2.1, 
#       text = paste("                                y =", y), adj = 0, outer = T, cex = 0.8) 
# mtext(side = 3, line = -3, 
#       text = paste("                                r =", r), adj = 0, outer = T, cex = 0.8) 
# mtext(side = 3, line = -3.9, 
#       text = paste("                         sigma =", sigma), adj = 0, outer = T, cex = 0.8, col = "red") 

# if (FigGen == TRUE) dev.off()

#### Discounted process

# if(FigGen == TRUE) pdf(paste("boundary_BB_y", yy, "_r", rr, "_sigma", ss, "_discounted.pdf", sep = ""))

# matplot(t, exp(-r * t) * t(BBs), type = "l", lty = 1, xlab = "Time", 
#         ylab = "Boundary", ylim = c(min(BBs, boundary), max(rbind(BBs, boundary))),
#         bty = "n")
# lines(t, boundary, lty = 1, lwd = 2, col = "red")
# lines(c(0, 1), c(y, y), lty = 2)
# title(main = "Bb | G(x) = x (discounted process)", line = 1.4, cex.main = 0.8)
# mtext(side = 3, line = -2.1, 
#       text = paste("                                y =", y), adj = 0, outer = T, cex = 0.8) 
# mtext(side = 3, line = -3, 
#       text = paste("                                r =", r), adj = 0, outer = T, cex = 0.8) 
# mtext(side = 3, line = -3.9, 
#       text = paste("                         sigma =", sigma), adj = 0, outer = T, cex = 0.8, col = "red") 

# if (FigGen == TRUE) dev.off()


