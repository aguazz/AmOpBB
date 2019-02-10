source(file = "0-simulations.R")

# - Loglikelihood
# Note: X is a list such that:
# X$path is the path f the GBM
# X$time is the partitoon where the path was sampled
# X$drift is the prehanded fixed drift of th GBM
GBM.lik <- function(vol, X){
  
  mu <- X$drift
  n <- length(X$path) - 1
  Delta_x <- diff(log(X$path))
  Delta_t <- diff(X$time)
  
  # d_logl <- sum( (Delta_x - (mu - vol^2 / 2) * Delta_x) ^ 2 / Delta_t ) - 
  #        (log(X$path[n + 1]) - log(X$path[1]) - 
  #          (mu - vol^2 / 2) * (X$time[n + 1] - X$time[n + 1]) + n) / vol  
  
  logl <- - n * log(vol) - 
    sum( (Delta_x - (mu - vol^2 / 2) * Delta_t)^2 / Delta_t ) / (2 * vol^2)  
  
  return(-logl)
  
}
##### Test ######
## rGBM
T <- 1
N <- 2e2
t_line <- seq(0, T, by = T / N)
GBM <- rGBM(n = 1, x = 10, T = T, N = N, mu = 1, sigma = 0.05)
matplot(t_line, t(GBM), type = "l", lty = 1)

output <- optim(par = 1, lower = 0.0001, fn = GBM.lik, 
                X = list(drift = 1, path = drop(GBM), time = t_line), 
                method = "L-BFGS-B")
print(paste("estimated volatility:", output$par))
