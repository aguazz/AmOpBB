## FigGen = TRUE is for generating the figures in the current path
# FigGen = FALSE makes the images show up in the plot zone
FigGen <- TRUE

## Compile C++ functions

library(Rcpp)
library(RcppArmadillo)
Sys.setenv("PKG_CXXFLAGS" = "-std=c++11") # Needs C++11

Rcpp::sourceCpp('1-boundary_computation.cpp')

## Load R functions

source("0-BB_simulations.R")
source("2-inference.R")

###### Libraries ######
library(latex2exp)
library(mvtnorm)
require(viridis)

###### Settings ######
a <- 0        # strike price
y <- 0        # expiration price
T <- 1        # expiration date
sigma <- 1    # volatility
N <- 2e2      # number of points (N + 1) used to discretize the interval [0, T]
h <- T / N    # lenght step used to discretize the interval [0, T]
t <- seq(0, T, by = h)  # discretization of the interval [0, T]
boundary.true <- sigma * 0.8399 * sqrt(rev(t))

###### Figures ######
## Figure 0: Brownian bridges
# n <- 3
# set.seed(89)
# BBs <- rBB(n = n, a = a, b = y, T = T, N = 400)
# t.line400 <- seq(0, T, by = T / 400)

# if(FigGen == TRUE) {
#   pdf("BB.pdf", width = 7.5, height = 5)
#   par( mfrow=c(1,1), mai=c(0.25, 0.25, 0.25,0.25), pin=c(6.8, 4.1), fin=c(7.5, 5), tcl=-0.3)
# }

# matplot(t.line400, t(BBs), type = "l", lty = 1, lwd = 1.5, bty="n", yaxt='n', ylab = "", xlab = "", xaxt='n', ylim = c(-1,1))
# axis(1, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), padj = -0.75)
# axis(2, at = c(-1, -0.5, 0, 0.5, 1), labels = TeX(c("-1","-0.5", "0", "0.5", "1")), padj = 0.3)

# if (FigGen == TRUE) dev.off()

## Figure 1:
## Asset's evolution
# set.seed(55)
# n <- 1
# samp <- rBB(n = n, N = N, b = 0.25)

# if(FigGen == TRUE) {
#   pdf("asset_evol.pdf", width = 6, height = 3)
#   par( mfrow=c(1,1), mai=c(0.25, 0.25, 0.25,0.25), pin=c(5.3, 2.3), fin=c(6, 3), tcl=-0.3)
# }

# matplot(t, t(samp), type = "l", lty = 1, bty="n", yaxt='n', xaxt='n', ylab = "", xlab = "")
# axis(1, at = c(0, 1), labels = TeX(c("$t$","$T$")), padj = -0.5, line = 0.3)
# axis(1, at = c(0.5), labels = TeX("$\\rightarrow s$"), tick = F, padj = -0.5, line = 0.3)
# axis(2, at = c(-0.3,0,0.5), labels = TeX(c("","$x$","")), padj = 1)
# title(main = "Asset evolution", line = 0.5, cex.main = 0.8)

# if(FigGen == TRUE) dev.off()

## Figure 2:
## Stopping/Continuation set
set.seed(71)
samp <- rBB(n = 1, N = N, b = y)
ost <- min(which(samp > boundary.true))

eps1 <- 0.05
t_left1 <- seq(min(samp) - eps1, boundary.true[1], l = 20)
t_right1 <- seq(boundary.true[N + 1], min(samp) - eps1, l = 6)
tt <- seq(0, T, l = 60)
fttC <- c( min(samp) - eps1, rep(c(min(samp) - eps1 - 0.01, min(samp) - eps1 + 0.01), 29), min(samp) - eps1)
fttD <- c( max(samp) + eps1, rep(c(max(samp) + eps1 - 0.01, max(samp) + eps1 + 0.01), 29), max(samp) + eps1)

color1 <- viridis(256, alpha = 0.5, option = "inferno")[125]
color2 <- viridis(256, alpha = 0.5, option = "viridis")[200]

if(FigGen == TRUE) {
  pdf("OSB.pdf", width = 7.5, height = 4.2)
  par( mfrow=c(1,1), mai=c(0.25, 0.25, 0.25,0.25), pin=c(6.7,3.6), fin=c(7.5, 4.2), tcl=-0.3)
}

plot(t, t(samp), type = "l", lty = 1, bty="n" , yaxt='n', ylab = "", xlab = "", xaxt='n', ylim = c(min(samp) - eps1, max(samp) + eps1), col = "blue", lwd = 1)
polygon(c(rev(tt), t), c(rev(fttD), boundary.true), col = color1, border = NA)
polygon(c(rev(tt), t), c(rev(fttC), boundary.true), col = color2, border = NA)
lines(t, boundary.true, col = "red", lwd = 1)
lines(t, t(samp), type = "l", lty = 1, col = "black", lwd = 1)
lines(c(t[ost], t[ost]), c(min(samp), samp[ost]), lty = 3, lwd = 1)
axis(1, at = c(0, t[ost], 1), labels = TeX(c("t", "$\\tau^{*}$","T")), padj = -0.8, line = 0.1)
points(0.798, 0.775, pch = "_")
points(0.8, 0.7, pch = "D")
points(0.1, 0.5, pch = "C")
axis(2, at = c(min(samp),0,max(samp)), labels = c("","x",""), padj = 1, line = -0.4)
# title(main = "Optimal Stopping Boundary (OSB)", line = 0.5, cex.main = 0.8)

if (FigGen == TRUE) dev.off()

## Figure 3:
## bBB performance (Log-time)
t20 <- log(seq(exp(0), exp(1), l = 21))
t50 <- log(seq(exp(0), exp(1), l = 51))
t100 <- log(seq(exp(0), exp(1), l = 101))
t200 <- log(seq(exp(0), exp(1), l = 201))

boundary20 <- bBBCpp(tol = 1e-3,y = y, sigma = sigma, t = t20)
boundary50 <- bBBCpp(tol = 1e-3,y = y, sigma = sigma, t = t50)
boundary100 <- bBBCpp(tol = 1e-3,y = y, sigma = sigma, t = t100)
boundary200 <- bBBCpp(tol = 1e-3,y = y, sigma = sigma, t = t200)

# N = 20
if(FigGen == TRUE) {
  pdf("algOSB_20.pdf", width = 5, height = 2.8)
  par( mfrow=c(1, 1), mai=c(0.25, 0.25, 0.25, 0.25), pin=c(4.2, 2.2), fin=c(5, 2.8), tcl = -0.3, cex.axis = 1)
}

plot(t, boundary.true, type = "l", lwd = 3, col = "red", bty="n", yaxt='n', xaxt='n', ylab = "", xlab = "")
lines(t20, boundary20, lty = 1, lwd = 1.5, col = "blue")
axis(1, at = c(0, 0.5, 1), labels = TeX(c("0", "0.5", "1")), padj = -0.8)
axis(2, at = c(0, 0.4, 0.8), labels = c("0", "0.4", "0.8"), padj = 1)
legend("topright", legend = c(TeX("$b$"), TeX("$\\tilde{b}$")), lty = c(1, 1), lwd = c(3, 1.5), bty = "n", col = c("red", "blue"), cex = 1.5)

if (FigGen == TRUE) dev.off()

# N = 50
if(FigGen == TRUE) {
  pdf("algOSB_50.pdf", width = 5, height = 2.8)
  par( mfrow=c(1, 1), mai=c(0.25, 0.25, 0.25, 0.25), pin=c(4.2, 2.2), fin=c(5, 2.8), tcl = -0.3)
}

plot(t, boundary.true, type = "l", lwd = 3, col = "red", bty="n", yaxt='n', xaxt='n', ylab = "", xlab = "")
lines(t50, boundary50, lty = 1, lwd = 1.5, col = "blue")
axis(1, at = c(0, 0.5, 1), labels = TeX(c("0", "0.5", "1")), padj = -0.8)
axis(2, at = c(0, 0.4, 0.8), labels = c("0", "0.4", "0.8"), padj = 1)
legend("topright", legend = c(TeX("$b$"), TeX("$\\tilde{b}$")), lty = c(1, 1), lwd = c(3, 1.5), bty = "n", col = c("red", "blue"), cex = 1.5)

if (FigGen == TRUE) dev.off()

# N = 100
if(FigGen == TRUE) {
  pdf("algOSB_100.pdf", width = 5, height = 2.8)
  par( mfrow=c(1, 1), mai=c(0.25, 0.25, 0.25, 0.25), pin=c(4.2, 2.2), fin=c(5, 2.8), tcl = -0.3)
}

plot(t, boundary.true, type = "l", lwd = 3, col = "red", bty="n", yaxt='n', xaxt='n', ylab = "", xlab = "")
lines(t100, boundary100, lty = 1, lwd = 1.5, col = "blue")
axis(1, at = c(0, 0.5, 1), labels = TeX(c("0", "0.5", "1")), padj = -0.8)
axis(2, at = c(0, 0.4, 0.8), labels = c("0", "0.4", "0.8"), padj = 1)
legend("topright", legend = c(TeX("$b$"), TeX("$\\tilde{b}$")), lty = c(1, 1), lwd = c(3, 1.5), bty = "n", col = c("red", "blue"), cex = 1.5)

if (FigGen == TRUE) dev.off()

# N = 200
if(FigGen == TRUE) {
  pdf("algOSB_200.pdf", width = 5, height = 2.8)
  par( mfrow=c(1, 1), mai=c(0.25, 0.25, 0.25, 0.25), pin=c(4.2, 2.2), fin=c(5, 2.8), tcl = -0.3)
}

plot(t, boundary.true, type = "l", lwd = 3, col = "red", bty="n", yaxt='n', xaxt='n', ylab = "", xlab = "")
lines(t200, boundary200, lty = 1, lwd = 1.5, col = "blue")
axis(1, at = c(0, 0.5, 1), labels = TeX(c("0", "0.5", "1")), padj = -0.8)
axis(2, at = c(0, 0.4, 0.8), labels = c("0", "0.4", "0.8"), padj = 1)
legend("topright", legend = c(TeX("$b$"), TeX("$\\tilde{b}$")), lty = c(1, 1), lwd = c(3, 1.5), bty = "n", col = c("red", "blue"), cex = 1.5)

if (FigGen == TRUE) dev.off()


## Figure 4:
## MLE for the volatility based on the sample size using n BB paths
set.seed(43) # ... forcing beauty
n <- 5

samp <- rBB(n = n, N = N, a = a, b = y, sigma = 1)

Sigma <- sapply(2:N, function(x) SigMl(samp = samp, N1 = x, T = T))

if(FigGen == TRUE) {
  pdf("vol_est.pdf", width = 7.5, height = 4.6)
  par( mfrow=c(1, 1), mai=c(0.25, 0.25, 0.25, 0.25), pin=c(6.7,3.8), fin=c(7.5, 4.6), tcl = -0.3, cex.axis = 1.2, cex.lab = 0.8)
}

matplot(1:(N - 1), t(Sigma), type = "l", lty = 1, lwd = 1.5, col = "blue", xlab = TeX("$k$"), ylab = "", bty="n", yaxt = "n", xaxt = "n")
lines(1:(N - 1),rep(sigma, N-1), lty = 2, lwd = 3, col = "red")
axis(1, at = c(0, 100, 200), labels = c("0", "100", "200"), padj = 0)
axis(2, at = c(0.25, 1, 3), labels = c("", "1", "3"), padj = 0)
legend("topright", legend = c(TeX("$\\widehat{\\sigma}_{k}$"), TeX("$\\sigma = 1$")), lty = c(1,2), lwd = c(1.5,2), col = c("blue", "red"), bty = "n", cex = 1.6)
# title(xlab = TeX("Sample size k"), line = 1.5)

if (FigGen == TRUE) dev.off()

## Figure 5:
## Inferring the boundary (including the confidence functions) using
## the first N1 observations of the path 
set.seed(88) # ... forcing beauty
n <- 1
alpha = 0.05
samp <- rBB(n = n, N = N, a = a, b = y, sigma = 1)
N1 <- floor(N / 3)
Sigma <- SigMl(samp = samp, T = T, N1 = N1)

tlog <- log(seq(exp(0), exp(1), l = N + 1)) 
  
boundary <- bBBCpp(tol = 1e-3, y = y, sigma = drop(Sigma), t = tlog)
boundary <- bBBConf(tol = 1e-3, bnd = boundary, y = y, sigma = drop(Sigma),
                    t = tlog, nsamp = N1, alpha = alpha, eps = 1e-2)

boundary$bBB <- splinefun(tlog, boundary$bBB)(t)
boundary$bBB.low <- splinefun(tlog, boundary$bBB.low)(t)
boundary$bBB.up <- splinefun(tlog, boundary$bBB.up)(t)

if(FigGen == TRUE) {
  pdf("boundary_inf.pdf", width = 5, height = 2.5)
  par( mfrow=c(1, 1), mai=c(0.25, 0.25, 0.25, 0.25), pin=c(4.2, 1.9), fin=c(5, 2.5), tcl = -0.3, cex.axis = 0.7)
}

plot(t, samp, type = "n", xlab = "", ylab = "", ylim = c(-0.5,0.75), xlim = c(0,1), bty = "n", yaxt = "n", xaxt = "n")
lines(t[1:N1], samp[1:N1], lty = 1, lwd = 1)
lines(t[N1:(N + 1)], boundary.true[N1:(N + 1)], lty = 1, lwd = 1, col = "red")
lines(t[N1:(N + 1)], samp[N1:(N + 1)], lty = 1, lwd = 0.5)
lines(t[N1:(N + 1)], boundary$bBB[N1:(N + 1)], lty = 1, lwd = 1, col = "blue")
lines(t[N1:(N + 1)], boundary$bBB.up[N1:(N + 1)], lty = 1, lwd = 1, col = "orange")
lines(t[N1:(N + 1)], boundary$bBB.low[N1:(N + 1)], lty = 1, lwd = 1, col = "green")
lines(c(t[N1],t[N1]), c(-0.4, 0.85), lty = 1)
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"), padj = -1,  line = -0.2)
axis(2, at = c(-0.4, 0, 0.4, 0.8), labels = c("-0.4", "0", "0.4", "0.8"), padj = 1)
legend(x = 0.005, y = 0.9, legend = c(TeX("$b_{\\sigma}$"), TeX("$\\tilde{b}_{\\widehat{\\sigma}_{k}}$"), TeX("$\\tilde{c}_{1,\\widehat{\\sigma}_{k}}$"), TeX("$\\tilde{c}_{2,\\widehat{\\sigma}_{k}}$")), lty = c(1, 1, 1, 1), lwd = c(1, 1 , 1, 1), col = c("red", "blue", "orange", "green"), bty = "n")

if (FigGen == TRUE) dev.off()

## Figure 6:
## BB percentile paths
prc02 <- qnorm( 0.2, mean = (a + (y - a) * t) / T, sd = sigma * sqrt(t * (T - t) / T ) ) 
prc04 <- qnorm( 0.4, mean = (a + (y - a) * t) / T, sd = sigma * sqrt(t * (T - t) / T ) ) 
prc06 <- qnorm( 0.6, mean = (a + (y - a) * t) / T, sd = sigma * sqrt(t * (T - t) / T ) ) 
prc08 <- qnorm( 0.8, mean = (a + (y - a) * t) / T, sd = sigma * sqrt(t * (T - t) / T ) ) 

n <- 4
set.seed(26)
BBs <- rBB(n = n, a = a, b = y, T = T, N = N)

if (FigGen == TRUE){
  pdf("percentiles.pdf", width = 4.5, height = 3)
  par( mfrow=c(1, 1), mai=c(0.25, 0.25, 0.25, 0.25), pin=c(4, 2.9), fin=c(4.5, 3), tcl = -0.3, cex.axis = 0.7, cex.lab = 0.8)
}

matplot(t, t(BBs), type = "l", col = "grey", lty = 1, xlab = "Time", ylab = "Percentile", ylim = c(y - 1, y + 1), yaxt = "n", xaxt = "n", bty = "n")
lines(t[1:90], prc08[1:90], lty = 1, lwd = 2)
lines(t[110:201], prc08[110:201], lty = 1, lwd = 2)
lines(t[1:90], prc06[1:90], lty = 1, lwd = 2)
lines(t[110:201], prc06[110:201], lty = 1, lwd = 2)
lines(t[1:90], prc04[1:90], lty = 1, lwd = 2)
lines(t[110:201], prc04[110:201], lty = 1, lwd = 2)
lines(t[1:90], prc02[1:90], lty = 1, lwd = 2)
lines(t[110:201], prc02[110:201], lty = 1, lwd = 2)
points(rbind(c(t[97], prc08[100]), c(t[97], prc06[100]), c(t[97], prc04[100]), c(t[97], prc02[100])), 
       pch = c("0","0","0","0"))
points(rbind(c(t[100] - 0.002, prc08[100] - .035), c(t[100] - 0.002, prc06[100] - .035), c(t[100] - 0.002, prc04[100] - .035), c(t[100] - 0.002, prc02[100] - .035)), 
       pch = c(".",".",".","."))
points(rbind(c(t[103], prc08[100]), c(t[103], prc06[100]), c(t[103], prc04[100]), c(t[103], prc02[100])), 
       pch = c("8","6","4","2"))
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"), padj = -0.7,  line = -1.5)
axis(2, at = c(y -0.75, y, y + 0.75), labels = c("-0.75", "0", "0.75"), padj = 0.6,  line = -0.5)

if (FigGen == TRUE) dev.off()

## Figure 7.1:
## Mean (Sigma = 1)
load(file = "payoff_alpha05_sigma1.RData")
x.axis <- t[-1]
percentiles <- c(0.2, 0.4, 0.6, 0.8)
y.lim <- c(.5, .5, .5, .5)

for (k in 1:4) {
  
  if(FigGen == TRUE) {
    pdf(paste("mean0", 2 * k,".pdf", sep = ""), width = 2.9, height = 1.7)
    par( mfrow=c(1, 1), mai=c(0.25, 0.25, 0.25, 0.25), pin=c(2.5, 1.1), fin=c(2.9, 1.7), tcl = -0.3, cex.axis = 0.6, cex.lab = 0.6)
  }
  
  plot(x.axis, colMeans(payoff_alpha05_sigma1$true[, , k]), type = "l", ylab = "", xlab = "", col = "red", lwd = 1.5, ylim = c(0, y.lim[k]), xlim = c(0,1), bty = "n", yaxt = "n", xaxt = "n")
  lines(x.axis, colMeans(payoff_alpha05_sigma1$est[, , k]), lty = 1, lwd = 1.5, col = "blue")
  lines(x.axis, colMeans(payoff_alpha05_sigma1$up[, , k]), lty = 1, lwd = 1.5,  col = "orange")
  lines(x.axis, colMeans(payoff_alpha05_sigma1$low[, , k]), lty = 1, lwd = 1.5, col = "green")
  axis(1, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), padj = -2.5,  line = 0.5)
  axis(2, at = c(0, 0.25, 0.5), labels = c("0", "0.25", "0.5"), padj = 1.6, line = 0)
  title(main = TeX(paste("$q =$", as.character(percentiles[k]))), cex.main = 0.7, line = 0.25)
  # legend("topright", legend = c(TeX("$b$"), TeX("$\\widetilde{b}_{\\widehat_{\\sigma}_{k,N}}$"), 
  #                              TeX("$\\widetilde{c}_{1,\\widehat{\\sigma}_{k,N},\\alpha}$"),
  #                              TeX("$\\widetilde{c}_{2,\\widehat{\\sigma}_{k,N},\\alpha}$")),
  #       lty = c(1,1,2,3), lwd = c(1, 1, 1.25, 1.25) , col = c("red", "blue", "blue", "blue"), bty = "n", cex = 0.7)
  
  if (FigGen == TRUE) dev.off()
  
}

## Figure 7.2:
## Variance (Sigma = 1)
# colVars (Auxiliary function)
colVars <- function(X){
  colSums( ( X - rep(1, nrow(X)) %*% t(colMeans(X)) ) ^ 2 ) / (nrow(X) - 1)
}

y.lim <- c(.15, .15, .15, .15)

for (k in 1:4) {
  
  if(FigGen == TRUE) {
    pdf(paste("variance0", 2 * k,".pdf", sep = ""), width = 2.9, height = 1.7)
    par( mfrow=c(1, 1), mai=c(0.25, 0.25, 0.25, 0.25), pin=c(2.5, 1.1), fin=c(2.9, 1.7), tcl = -0.3, cex.axis = 0.6, cex.lab = 0.6)
  }
  
  plot(x.axis, colVars(payoff_alpha05_sigma1$true[, , k]), type = "l", ylab = "", xlab = "", col = "red", lwd = 1.5, ylim = c(0, y.lim[k]), xlim = c(0,1), bty = "n", yaxt = "n", xaxt = "n")
  lines(x.axis, colVars(payoff_alpha05_sigma1$est[, , k]), lty = 1, lwd = 1.5, col = "blue")
  lines(x.axis, colVars(payoff_alpha05_sigma1$up[, , k]), lty = 1, lwd = 1.5,  col = "orange")
  lines(x.axis, colVars(payoff_alpha05_sigma1$low[, , k]), lty = 1, lwd = 1.5, col = "green")
  axis(1, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), padj = -2.5,  line = 0.5)
  axis(2, at = c(0, .05, 0.1, .15), labels = c("0", "0.05", "0.1", "0.15"), padj = 1.6, line = 0)
  title(main = TeX(paste("$q  =$", as.character(percentiles[k]))), cex.main = 0.7, line = 0.25)
  # legend("topright", legend = c(TeX("$b$"), TeX("$\\widetilde{b}_{\\widehat_{\\sigma}_{k,N}}$"), 
  #                              TeX("$\\widetilde{c}_{1,\\widehat{\\sigma}_{k,N},\\alpha}$"),
  #                              TeX("$\\widetilde{c}_{2,\\widehat{\\sigma}_{k,N},\\alpha}$")),
  #       lty = c(1,1,2,3), lwd = c(1, 1, 1.25, 1.25) , col = c("red", "blue", "blue", "blue"), bty = "n", cex = 0.7)
  
  if (FigGen == TRUE) dev.off()
  
}

## Figure 8
## Empirical validation of the confidence functions
## by showing the sample proportion of inclusion of the true boundary.
# settings
n <- 1000 # number of BB paths
alpha <- c(0.1, 0.05) # confidence level
nalpha <- length(alpha)
# Number of points for the sigma estimation.
N1 <- floor(N/3) 
N2 <- 2*floor(N/3)

# CI for the proportion of contentions of the boundary
up.bound <- alpha + qnorm(1 - alpha/2, mean = 0, sd = sqrt(alpha * (1 - alpha) / n)) 
low.bound <- alpha - qnorm(1 - alpha/2, mean = 0, sd = sqrt(alpha * (1 - alpha) / n))

set.seed(1234) #... forcing beauty
# BB paths
samp <- rBB(n = n, a = a, b = y, sigma = sigma, T = 1, N = N)  

Sigma1 <- drop(SigMl(samp = samp, T = T, N1 = N1))
Sigma2 <- drop(SigMl(samp = samp, T = T, N1 = N2))

# Defining logarithmic partition
tlog <- log(seq(exp(0), exp(1), l = N + 1))

# true boundary at the logarithmic partition
boundary.true <- 0.8399*sqrt(T - tlog)

# Generate the boundary and confidence functions...
boundary1 <- bBBCpp(y = y, sigma = Sigma1, t = tlog)
boundary1 <- bBBConf(bnd = boundary1, y = y, sigma = Sigma1, t = tlog, nsamp = N1, alpha = alpha)

boundary2 <- bBBCpp(y = y, sigma = Sigma2, t = tlog)
boundary2 <- bBBConf(bnd = boundary2, y = y, sigma = Sigma2, t = tlog, nsamp = N2, alpha = alpha)

# Estimated proportion of successes (boundary within the confidence interval (pointwise))
Prop1 <- sapply(1:nalpha, function(i) {
  rowMeans( t(boundary1$bBB.low[, , i]) <= boundary.true & t(boundary1$bBB.up[, , i]) >= boundary.true)
})

Prop2 <- sapply(1:nalpha, function(i) {
  rowMeans( t(boundary2$bBB.low[, , i]) <= boundary.true & t(boundary2$bBB.up[, , i]) >= boundary.true)
})

# Save the data
# save(Prop1, file = "Prop1.RData")
# save(Prop2, file = "Prop2.RData")
# ... or load it 
load(file = "Prop1.RData")
load(file = "Prop2.RData")

# Images 8.1 (N1)
for (i in 1:2) {

  if(FigGen == TRUE) {
    pdf(paste("confidence_1third_alpha", alpha[i],".pdf", sep = ""), width = 5, height = 2.8)
    par( mfrow=c(1, 1), mai=c(0.25, 0.25, 0.25, 0.25), pin=c(2.5, 1.1), fin=c(2.9, 1.7), tcl = -0.3, cex.axis = 0.6, cex.lab = 0.6)
  }

  plot(tlog, rep(alpha[i], N + 1), type = "l", lty = 2, lwd = 3, ylim = c(0, 2*alpha[i]), ylab = "Proportion", xlab = "Time", bty = "n", yaxt = "n", xaxt = "n")
  lines(tlog[1:(N + 1)], 1 - Prop1[1:(N + 1), i], lty = 1, lwd = 1.5)
  lines(c(0, T), rep(min(1, up.bound[i]), 2), lty = 3, lwd = 2)
  lines(c(0, T), rep(max(0 , low.bound[i]), 2), lty = 3, lwd = 2)
  axis(1, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), padj = -0.8)
  axis(2, at = c(0, alpha[i], 2*alpha[i]), labels = c("0", as.character(alpha[i]), as.character(2*alpha[i])), padj = 1)  
  if (FigGen == TRUE) dev.off()
  
}

# Images 8.2 (N2)
for (i in 1:2) {
  
  if(FigGen == TRUE) {
    pdf(paste("confidence_2third_alpha", alpha[i],".pdf", sep = ""), width = 5, height = 2.8)
    par( mfrow=c(1, 1), mai=c(0.25, 0.25, 0.25, 0.25), pin=c(4.2, 2.2), fin=c(5, 2.8), tcl = -0.3)
  }
  
  plot(t, rep(alpha[i], N + 1), type = "l", lty = 2, lwd = 3, ylim = c(0, 2*alpha[i]), ylab = "Proportion", xlab = "Time", bty = "n", yaxt = "n", xaxt = "n")
  lines(tlog[1:(N + 1)], 1 - Prop2[1:(N + 1), i], lty = 1, lwd = 1.5)
  lines(c(0, T), rep(min(1, up.bound[i]), 2), lty = 3, lwd = 2)
  lines(c(0, T), rep(max(0 , low.bound[i]), 2), lty = 3, lwd = 2)
  axis(1, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), padj = -0.8)
  axis(2, at = c(0, alpha[i], 2*alpha[i]), labels = c("0", as.character(alpha[i]), as.character(2*alpha[i])), padj = 1)
  
  if (FigGen == TRUE) dev.off()
  
}
