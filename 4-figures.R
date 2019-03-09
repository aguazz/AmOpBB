# FigGen = TRUE is for generating the figures in the current path
# FigGen = FALSE makes the images show up in the plot zone
FigGen <- FALSE

## Load R functions
source("0-simulations.R")
source("1-boundary_computation_BB-AmPut_discount.R")
source("2-inference_BB.R")

###### Libraries ######
library(latex2exp)
library(mvtnorm)
require(viridis)

###### Settings ######
a <- 10       # initial price
S <- 10       # Strike price
r <- 0        # discount rate
T <- 1        # expiration date
sigma <- 1    # volatility
N <- 2e2      # number of points (N + 1) used to discretize the interval [0, T]
h <- T / N    # lenght step used to discretize the interval [0, T]
t <- seq(0, T, by = h)  # discretization of the interval [0, T]
boundary.true <- S - sigma * 0.8399 * sqrt(T - t)

###### Figures ######
## Figure 1:
## Stopping/Continuation set
T <- 1        # expiration date
sigma <- 1    # volatility
N <- 5e2      # number of points (N + 1) used to discretize the interval [0, T]
t <- seq(0, T, by = 1/N)  # discretization of the interval [0, T]
S <- 0        # Strike price

boundary <- sigma * 0.8399 * sqrt(T - t) # true boundary

set.seed(2)
samp <- rBB(n = 1, N = N, b = S)
ost <- min(which(samp > boundary))


eps1 <- 0.05
t_left1 <- seq(min(samp) - eps1, boundary[1], l = 20)
t_right1 <- seq(boundary[N + 1], min(samp) - eps1, l = 6)
tt <- seq(0, T, l = 60)
fttC <- c( min(samp, boundary) - eps1, 
           rep(c(min(samp, boundary) - eps1 - 0.01, min(samp, boundary) - eps1 + 0.01), 29), 
           min(samp,boundary) - eps1)
fttD <- c( max(samp, boundary) + eps1, 
           rep(c(max(samp, boundary) + eps1 - 0.01, max(samp, boundary) + eps1 + 0.01), 29), 
           max(samp, boundary) + eps1)

color1 <- viridis(256, alpha = 0.5, option = "inferno")[125]
color2 <- viridis(256, alpha = 0.5, option = "viridis")[200]

if(FigGen == TRUE) {
  pdf("OSB.pdf", width = 7.5, height = 4.2)
  par(mfrow=c(1,1), mai=c(0.25, 0.25, 0.25,0.25), 
       pin=c(6.7,3.6), fin=c(7.5, 4.2), tcl=-0.3)
}

plot(t, t(samp), type = "l", lty = 1, bty="n" , 
     yaxt='n', ylab = "", xlab = "", xaxt='n', 
     ylim = c(min(samp, boundary) - eps1, max(samp, boundary) + eps1), 
     col = "blue", lwd = 1)
polygon(c(rev(tt), t), c(rev(fttD), boundary), col = color1, border = NA)
polygon(c(rev(tt), t), c(rev(fttC), boundary), col = color2, border = NA)
lines(t, boundary, col = "red", lwd = 1)
lines(t, t(samp), type = "l", lty = 1, col = "black", lwd = 1)
lines(c(t[ost], t[ost]), c(min(samp), samp[ost]), lty = 3, lwd = 1)
axis(1, at = c(0, t[ost], 1), labels = TeX(c("$t$", "$\\tau^{*}$","$T$")), padj = -0.8, line = 0.1)
points(0.89, 0.9, pch = "D")
points(0.6, 0.15, pch = "C")
axis(2, at = c(min(samp),0,max(samp)), labels = TeX(c("","$x$","")), padj = 1, line = -0.4)

if (FigGen == TRUE) dev.off()

## Figure 2:
## bBB performance (Log-time)
T <- 1        # expiration date
sigma <- 1    # volatility
S <- 10       # strike price

t20 <- log(seq(exp(0), exp(1), l = 21))
t50 <- log(seq(exp(0), exp(1), l = 51))
t100 <- log(seq(exp(0), exp(1), l = 101))
t200 <- log(seq(exp(0), exp(1), l = 201))

boundary20 <- b_BB_AmPut(S = 10, sigma = 1, t = t20, r = 0)
boundary50 <- b_BB_AmPut(S = 10, sigma = 1, t = t50, r = 0)
boundary100 <- b_BB_AmPut(S = 10, sigma = 1, t = t100, r = 0)
boundary200 <- b_BB_AmPut(S = 10, sigma = 1, t = t200, r = 0)

tt <- seq(0, T, l = 1000)
boundary.true <- S - 0.8399 * sigma * sqrt(T - tt)

# N = 20
if(FigGen == TRUE) {
  pdf("algOSB_20.pdf", width = 5, height = 2.8)
  par( mfrow=c(1, 1), mai=c(0.25, 0.25, 0.25, 0.25), pin=c(4.2, 2.2), fin=c(5, 2.8), tcl = -0.3, cex.axis = 1)
}

plot(tt, boundary.true, type = "l", lwd = 3, 
     col = "red", bty = "n", ylim = c(9, 10.2), yaxt = 'n', xaxt='n', 
     ylab = "", xlab = "")
lines(t20, boundary20, lty = 1, lwd = 1.5, col = "blue")
lines(c(0, 1), c(S, S), lty = 3)
axis(1, at = c(0, 0.5, 1), labels = TeX(c("0", "0.5", "1")), padj = -0.8)
axis(2, at = c(9, 9.5, 10, 10.2), labels = c("9", "9.5", "10", ""), padj = 1)
rug(t20, ticksize = 0.05)
legend(x = 0.82, y = 9.55, legend = c(TeX("$b$"), TeX("$\\tilde{b}$")), 
       lty = c(1, 1), lwd = c(3, 1.5), col = c("red", "blue"), cex = 1)

if (FigGen == TRUE) dev.off()

# N = 50
if(FigGen == TRUE) {
  pdf("algOSB_50.pdf", width = 5, height = 2.8)
  par( mfrow=c(1, 1), mai=c(0.25, 0.25, 0.25, 0.25), pin=c(4.2, 2.2), fin=c(5, 2.8), tcl = -0.3)
}

plot(tt, boundary.true, type = "l", lwd = 3, 
     col = "red", bty = "n", ylim = c(9, 10.2), yaxt = 'n', xaxt='n', 
     ylab = "", xlab = "")
lines(t50, boundary50, lty = 1, lwd = 1.5, col = "blue")
lines(c(0, 1), c(S, S), lty = 3)
axis(1, at = c(0, 0.5, 1), labels = TeX(c("0", "0.5", "1")), padj = -0.8)
axis(2, at = c(9, 9.5, 10, 10.2), labels = c("9", "9.5", "10", ""), padj = 1)
rug(t50, ticksize = 0.05)

if (FigGen == TRUE) dev.off()

# N = 100
if(FigGen == TRUE) {
  pdf("algOSB_100.pdf", width = 5, height = 2.8)
  par( mfrow=c(1, 1), mai=c(0.25, 0.25, 0.25, 0.25), pin=c(4.2, 2.2), fin=c(5, 2.8), tcl = -0.3)
}

plot(tt, boundary.true, type = "l", lwd = 3, 
     col = "red", bty = "n", ylim = c(9, 10.2), yaxt = 'n', xaxt='n', 
     ylab = "", xlab = "")
lines(t100, boundary100, lty = 1, lwd = 1.5, col = "blue")
lines(c(0, 1), c(S, S), lty = 3)
axis(1, at = c(0, 0.5, 1), labels = TeX(c("0", "0.5", "1")), padj = -0.8)
axis(2, at = c(9, 9.5, 10, 10.2), labels = c("9", "9.5", "10", ""), padj = 1)
rug(t100, ticksize = 0.05)

if (FigGen == TRUE) dev.off()

# N = 200
if(FigGen == TRUE) {
  pdf("algOSB_200.pdf", width = 5, height = 2.8)
  par( mfrow=c(1, 1), mai=c(0.25, 0.25, 0.25, 0.25), pin=c(4.2, 2.2), fin=c(5, 2.8), tcl = -0.3)
}

plot(tt, boundary.true, type = "l", lwd = 3, 
     col = "red", bty = "n", ylim = c(9, 10.2), yaxt = 'n', xaxt='n', 
     ylab = "", xlab = "")
lines(t200, boundary200, lty = 1, lwd = 1.5, col = "blue")
lines(c(0, 1), c(S, S), lty = 3)
axis(1, at = c(0, 0.5, 1), labels = TeX(c("0", "0.5", "1")), padj = -0.8)
axis(2, at = c(9, 9.5, 10, 10.2), labels = c("9", "9.5", "10", ""), padj = 1)
rug(t200, ticksize = 0.05)

if (FigGen == TRUE) dev.off()

## Figure 3:
## Inferring the boundary (including the confidence functions) using
## the first N1 observations of the path 
set.seed(49) # ... forcing beauty

T <- 1        # expiration date
sigma <- 1    # volatility
N <- 2e2      # number of points (N + 1) used to discretize the interval [0, T]
S <- 10       # Strike price
t <- seq(0, T, by = 1/N)  # discretization of the interval [0, T]
boundary.true <- S - 0.8399 * sigma * sqrt(T - t)

n <- 1 # number of paths
alpha = 0.05 # confidence level
samp <- rBB(n = n, N = N, a = S, b = S, sigma = 1)
N1 <- floor(N / 3)
Sigma <- SigMl(samp = samp, T = T, N1 = N1)

tlog <- log(seq(exp(0), exp(1), l = N + 1)) 
  
boundary <- b_BB_AmPut(S = S, sigma = drop(Sigma), t = tlog, r = 0)
boundary <- b_BB_AmPut_Conf(bnd = boundary, S = S, r = 0, sigma = drop(Sigma), 
                            t = tlog, nsamp = N1, alpha = alpha, eps = 0.001)

boundary$bBB <- splinefun(tlog, boundary$bBB)(t)
boundary$bBB.low <- splinefun(tlog, boundary$bBB.low)(t)
boundary$bBB.up <- splinefun(tlog, boundary$bBB.up)(t)

if(FigGen == TRUE) {
  pdf("boundary_inf.pdf", width = 5, height = 2.5)
  par( mfrow=c(1, 1), mai=c(0.25, 0.25, 0.25, 0.25), pin=c(4.2, 1.9), fin=c(5, 2.5), tcl = -0.3, cex.axis = 0.7)
}

plot(t, samp, type = "n", xlab = "", ylab = "", ylim = c(8.5, 10.5), 
     xlim = c(0, 1), bty = "n", yaxt = "n", xaxt = "n")
lines(t[1:N1], samp[1:N1], lty = 1, lwd = 1.5)
lines(t[N1:(N + 1)], boundary.true[N1:(N + 1)], lty = 1, lwd = 1, col = "red")
lines(t[N1:(N + 1)], samp[N1:(N + 1)], lty = 1, lwd = 0.5)
lines(t[N1:(N + 1)], boundary$bBB[N1:(N + 1)], lty = 1, lwd = 1, col = "blue")
lines(t[N1:(N + 1)], boundary$bBB.up[N1:(N + 1)], lty = 1, lwd = 1, col = "orange")
lines(t[N1:(N + 1)], boundary$bBB.low[N1:(N + 1)], lty = 1, lwd = 1, col = "green")
lines(c(t[N1], t[N1]), c(8.8, 10.5), lty = 1)
lines(c(0, 1), c(S, S), lty = 2)
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
     labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"), padj = -1,  line = -0.2)
axis(2, at = c(8.5, 9, 9.5, 10, 10.5, 11), 
     labels = c("8.5", "9", "9.5", "10", "10.5", "11"), padj = 1)
legend(x = 0.81, y = 9.46, legend = c(TeX("$b_{\\sigma}$"), 
                                 TeX("$\\tilde{b}_{\\widehat{\\sigma}_{n}}$"), 
                                 TeX("$\\tilde{c}_{1,\\widehat{\\sigma}_{n}}$"), 
                                 TeX("$\\tilde{c}_{2,\\widehat{\\sigma}_{n}}$")), 
       lty = c(1, 1, 1, 1), lwd = c(1, 1 , 1, 1), 
       col = c("red", "blue", "orange", "green"), cex = 0.75)

if (FigGen == TRUE) dev.off()

## Figure 4
## Empirical validation of the confidence functions
## by showing the sample proportion of inclusion of the true boundary.
# settings
a <- 10   # initial value
S <- 10   # strike price
sigma <- 1 # volatility
T <- 1     # expiration date
N <- 5e2   # number of subintervals of the interval [0, T]


n <- 1000 # number of BB paths
alpha <- 0.05 # confidence level
nalpha <- length(alpha)
# Number of points for the sigma estimation.
N1 <- floor(N/3) 
N2 <- 2 * floor(N/3)

# CI (with confidence level of 0.05) for the proportion of contentions of the boundary
# within the confidence curves
up.bound <- alpha + qnorm(1 - 0.05/2, mean = 0, 
                          sd = sqrt(alpha * (1 - alpha) / n)) 
low.bound <- alpha - qnorm(1 - 0.05/2, mean = 0, 
                           sd = sqrt(alpha * (1 - alpha) / n))

set.seed(234)
# BB paths
samp <- rBB(n = n, a = a, b = S, sigma = sigma, T = 1, N = N)  

Sigma1 <- drop(SigMl(samp = samp, T = T, N1 = N1))
Sigma2 <- drop(SigMl(samp = samp, T = T, N1 = N2))

# Defining logarithmic partition
tlog <- log(seq(exp(0), exp(1), l = N + 1))

# true boundary at the logarithmic partition
boundary.true <- S - 0.8399 * sigma * sqrt(T - tlog)

# Generate the boundary and confidence functions...
# boundary1 <- b_BB_AmPut(S = S, sigma = Sigma1, t = tlog, r = 0)
# boundary1 <- b_BB_AmPut_Conf(bnd = boundary1, S = S, sigma = Sigma1, 
#                              t = tlog, nsamp = N1, alpha = alpha, r = 0)

# boundary2 <- b_BB_AmPut(S = S, sigma = Sigma2, t = tlog, r = 0)
# boundary2 <- b_BB_AmPut_Conf(bnd = boundary2, S = S, sigma = Sigma2, 
#                              t = tlog, nsamp = N2, alpha = alpha, r = 0)

# Estimated proportion of successes (boundary within the confidence interval (pointwise))
# Prop1 <- rowMeans( t(boundary1$bBB.low) <= boundary.true & t(boundary1$bBB.up) >= boundary.true )
# Prop2 <- rowMeans( t(boundary2$bBB.low) <= boundary.true & t(boundary2$bBB.up) >= boundary.true )

# Save the data
# save(Prop1, file = "Prop1.RData")
# save(Prop2, file = "Prop2.RData")
# ... or load it 
# Load the data: proportions of inclusions
load(file = "Prop1.RData")
load(file = "Prop2.RData")

alpha_after_dot <- strsplit(as.character(alpha), "")[[1]]
alpha_after_dot <- alpha_after_dot[(which(alpha_after_dot == ".") + 1):length(alpha_after_dot)]
alpha_after_dot <- paste(alpha_after_dot, collapse = "")

# Images 4.1 (N1)
if(FigGen == TRUE) {
  pdf(paste("confidence_1third_alpha", alpha_after_dot,".pdf", sep = ""), width = 5, height = 2.8)
  par( mfrow=c(1, 1), mai=c(0.25, 0.25, 0.25, 0.25), pin=c(4.2, 2.2), fin=c(5, 2.8), tcl = -0.3, cex.axis = 1)
}

plot(tlog, rep(alpha, N + 1), type = "l", lty = 2, lwd = 3, ylim = c(0, 2 * alpha), ylab = "", xlab = "", bty = "n", yaxt = "n", xaxt = "n")
lines(tlog[1:(N + 1)], 1 - Prop1[1:(N + 1)], lty = 1, lwd = 1.5)
lines(c(0, T), rep(min(1, up.bound), 2), lty = 3, lwd = 2)
lines(c(0, T), rep(max(0 , low.bound), 2), lty = 3, lwd = 2)
axis(1, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), padj = -0.8)
axis(2, at = c(0, alpha, 2*alpha), labels = c("0", as.character(alpha), as.character(2*alpha)), padj = 1)  
title(xlab = TeX("Time"), line = 2, cex.lab = 1.2)
title(ylab = TeX("Proportion of non containtion"), line = 2, cex.lab = 1.2)

if (FigGen == TRUE) dev.off()

# Images 4.2 (N2)
if(FigGen == TRUE) {
  pdf(paste("confidence_2third_alpha", alpha_after_dot,".pdf", sep = ""),  width = 5, height = 2.8)
  par( mfrow=c(1, 1), mai=c(0.25, 0.25, 0.25, 0.25), pin=c(4.2, 2.2), fin=c(5, 2.8), tcl = -0.3, cex.axis = 1)
}

plot(tlog, rep(alpha, N + 1), type = "l", lty = 2, lwd = 3, ylim = c(0, 2 * alpha), ylab = "", xlab = "", bty = "n", yaxt = "n", xaxt = "n")
lines(tlog[1:(N + 1)], 1 - Prop2[1:(N + 1)], lty = 1, lwd = 1.5)
lines(c(0, T), rep(min(1, up.bound), 2), lty = 3, lwd = 2)
lines(c(0, T), rep(max(0 , low.bound), 2), lty = 3, lwd = 2)
axis(1, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), padj = -0.8)
axis(2, at = c(0, alpha, 2*alpha), labels = c("0", as.character(alpha), as.character(2*alpha)), padj = 1)  
title(xlab = TeX("Time"), line = 2.9, cex.lab = 1.2)
title(ylab = TeX("Proportion of non containtion"), line = 2.6, cex.lab = 1.2)

if (FigGen == TRUE) dev.off()

## Figure 5:
## BB percentile paths
a <- 10   # initial value
S <- 10   # final price
T <- 1    # expiration date
N <- 4e2 # number of subintervals
t <- seq(0, T, by = T / N)  # discrite partiton of the interval [0, T]
sigma <- 1

prc02 <- qnorm( 0.2, mean = (a + (S - a) * t) / T, sd = sigma * sqrt(t * (T - t) / T ) ) 
prc04 <- qnorm( 0.4, mean = (a + (S - a) * t) / T, sd = sigma * sqrt(t * (T - t) / T ) ) 
prc06 <- qnorm( 0.6, mean = (a + (S - a) * t) / T, sd = sigma * sqrt(t * (T - t) / T ) ) 
prc08 <- qnorm( 0.8, mean = (a + (S - a) * t) / T, sd = sigma * sqrt(t * (T - t) / T ) ) 

n <- 3

set.seed(134)
idx1 <- floor(length(t) / 5)
guided_BBs1 <- rBBB(n = n, a = a, b = prc02[idx1], c = S, sigma = 1, t = t, N1 = idx1)

set.seed(444)
idx2 <- floor(4 * length(t) / 5)
guided_BBs2 <- rBBB(n = n, a = a, b = prc08[idx2], c = S, sigma = 1, t = t, N1 = idx2)


if (FigGen == TRUE){
  pdf("percentiles.pdf", width = 4.5, height = 3.5)
  par( mfrow=c(1, 1), mai=c(0.25, 0.25, 0.25, 0.25), pin=c(4, 3), fin=c(4.5, 3.5), tcl = -0.3, cex.axis = 0.7, cex.lab = 0.8)
}

matplot(t, t(guided_BBs1), type = "l", col = viridis::viridis(200)[150], 
        lty = 1, xlab = "Time", ylab = "Percentile", lwd = 1, 
        ylim = c(S - 1, S + 1), yaxt = "n", xaxt = "n", bty = "n")
matlines(t, t(guided_BBs2), lty = 1, lwd = 1, col = viridis::magma(200)[150])
points(t[idx1], prc02[idx1], pch = 20)
points(t[idx2], prc08[idx2], pch = 20)
lines(c(0, T), c(S, S), lty = 3)
lines(t[1:185], prc08[1:185], lty = 1, lwd = 1)
lines(t[215:401], prc08[215:401], lty = 1, lwd = 1)
lines(t[1:185], prc06[1:185], lty = 1, lwd = 1)
lines(t[215:401], prc06[215:401], lty = 1, lwd = 1)
lines(t[1:185], prc04[1:185], lty = 1, lwd = 1)
lines(t[215:401], prc04[215:401], lty = 1, lwd = 1)
lines(t[1:185], prc02[1:185], lty = 1, lwd = 1)
lines(t[215:401], prc02[215:401], lty = 1, lwd = 1)
points(rbind(c(t[195], prc08[200]), c(t[195], prc06[200]), 
             c(t[195], prc04[200]), c(t[195], prc02[200])), 
       pch = c("0","0","0","0"))
points(rbind(c(t[200] - 0.002, prc08[200] - .035), 
             c(t[200] - 0.002, prc06[200] - .035), 
             c(t[200] - 0.002, prc04[200] - .035), 
             c(t[200] - 0.002, prc02[200] - .035)), 
       pch = c(".",".",".","."))
points(rbind(c(t[205], prc08[200]), c(t[205], prc06[200]), 
             c(t[205], prc04[200]), c(t[205], prc02[200])), 
       pch = c("8","6","4","2"))
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"), padj = -0.7,  line = -1)
axis(2, at = c(S -0.75, S, S + 0.75), labels = c("9.25", "10", "10.75"), padj = 0.6,  line = -0.5)

if (FigGen == TRUE) dev.off()

## Figure 6.1:
## Mean (Sigma = 1)
load(file = "payoff_N200_ratio1.RData")
# load(file = "payoff_N200_ratio25.RData")
r <- 1   # ratio used for the simulated data
T <- 1   # expiration date
N <- 2e2 # number of subintervals
t <- seq(0, T, by = T / N)  # discrite partiton of the interval [0, T]
x.axis <- t[-1]
percentiles <- c(0.2, 0.4, 0.6, 0.8)
y.lim <- c(.5, .5, .5, .5)

for (k in 1:4) {
  
  if(FigGen == TRUE) {
    pdf(paste("mean0", 2 * k,"LOW.pdf", sep = ""), width = 2.9, height = 1.7)
    par( mfrow=c(1, 1), mai=c(0.25, 0.25, 0.25, 0.25), pin=c(2.5, 1.1), 
         fin=c(2.9, 1.7), tcl = -0.3, cex.axis = 0.6, cex.lab = 0.6)
  }
  
  plot(x.axis, colMeans(payoff_N200_ratio1$true[, , k]), type = "l", ylab = "", 
       xlab = "", col = "red", lwd = 1.5, ylim = c(0, y.lim[k]), 
       xlim = c(0,1), bty = "n", yaxt = "n", xaxt = "n")
  lines(x.axis, colMeans(payoff_N200_ratio1$est[, , k]), 
        lty = 1, lwd = 1.5, col = "blue")
  lines(x.axis, colMeans(payoff_N200_ratio1$up[, , k]), 
        lty = 1, lwd = 1.5,  col = "orange")
  lines(x.axis, colMeans(payoff_N200_ratio1$low[, , k]), 
        lty = 1, lwd = 1.5, col = "green")
  axis(1, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), padj = -2.5,  line = 0.5)
  axis(2, at = c(0, 0.25, 0.5), labels = c("0", "0.25", "0.5"), padj = 1.6, line = 0)
  title(main = TeX(paste("$q =$", as.character(percentiles[k]))), cex.main = 0.7, line = 0.25)
  if(k == 1 && r == 1){
  legend(x = 0.19, y = 0.38, legend = c(TeX("$b_{\\sigma}$"), 
                                   TeX("$\\tilde{b}_{\\widehat{\\sigma}_{n}}$"), 
                                   TeX("$\\tilde{c}_{1,\\widehat{\\sigma}_{n}}$"), 
                                   TeX("$\\tilde{c}_{2,\\widehat{\\sigma}_{n}}$")), 
         lty = c(1, 1, 1, 1), lwd = c(1, 1 , 1, 1), 
         col = c("red", "blue", "orange", "green"), cex = 0.75)
  }
  
  if (FigGen == TRUE) dev.off()
  
}

## Figure 6.2:
## Variance (Sigma = 1)
# colVars (Auxiliary function)
colVars <- function(X){
  colSums( ( X - rep(1, nrow(X)) %*% t(colMeans(X)) ) ^ 2 ) / (nrow(X) - 1)
}

y.lim <- c(.15, .15, .15, .15)

for (k in 1:4) {
  
  if(FigGen == TRUE) {
    pdf(paste("variance0", 2 * k,"LOW.pdf", sep = ""), width = 2.9, height = 1.7)
    par( mfrow=c(1, 1), mai=c(0.25, 0.25, 0.25, 0.25), pin=c(2.5, 1.1), 
         fin=c(2.9, 1.7), tcl = -0.3, cex.axis = 0.6, cex.lab = 0.6)
  }
  
  plot(x.axis, colVars(payoff_N200_ratio1$true[, , k]), type = "l", ylab = "", xlab = "", col = "red", lwd = 1.5, ylim = c(0, y.lim[k]), xlim = c(0,1), bty = "n", yaxt = "n", xaxt = "n")
  lines(x.axis, colVars(payoff_N200_ratio1$est[, , k]), lty = 1, lwd = 1.5, col = "blue")
  lines(x.axis, colVars(payoff_N200_ratio1$up[, , k]), lty = 1, lwd = 1.5,  col = "orange")
  lines(x.axis, colVars(payoff_N200_ratio1$low[, , k]), lty = 1, lwd = 1.5, col = "green")
  axis(1, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), padj = -2.5,  line = 0.5)
  axis(2, at = c(0, .05, 0.1, .15), labels = c("0", "0.05", "0.1", "0.15"), padj = 1.6, line = 0)
  title(main = TeX(paste("$q  =$", as.character(percentiles[k]))), cex.main = 0.7, line = 0.25)
  if(k == 1 && r == 1){
    legend(x = 0.73, y = 0.15, legend = c(TeX("$b_{\\sigma}$"), 
                                       TeX("$\\tilde{b}_{\\widehat{\\sigma}_{n}}$"), 
                                       TeX("$\\tilde{c}_{1,\\widehat{\\sigma}_{n}}$"), 
                                       TeX("$\\tilde{c}_{2,\\widehat{\\sigma}_{n}}$")), 
           lty = c(1, 1, 1, 1), lwd = c(1, 1 , 1, 1), 
           col = c("red", "blue", "orange", "green"), cex = 0.65)
  }
  
  if (FigGen == TRUE) dev.off()
  
}
  