##### Aux function ######
subset_string <- function(str1, rng1, rng2) {
  
  if(rng2 == "end"){
    aux <- strsplit(str1, split = "")[[1]]
    return(paste(aux[rng1:length(aux)], collapse = ""))
  }
  return(paste(strsplit(str1, split = "")[[1]][rng1:rng2], collapse = ""))
  
}
###### ######
###### Sigma estimation ######
### GBM
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

### BB
source("2-inference_BB.R")
###### ######
###### Data processing ######
# Loading intraday data...
apple_intraday_prices <- read.csv(file = "apple_intraday_prices.csv", 
                              header = TRUE, stringsAsFactors = FALSE) # intraday prices

# Get the position of the rows with dates (in string format)
where_dates <- Map(apple_intraday_prices[, 1], f = function(x) {
  length(strsplit(x, split = "")[[1]]) == 25 
})
where_dates <- unname(unlist(where_dates)) # turn from list to vector
where_dates <- which(where_dates)          # get the position
# Get the dates (in string format)
intra_dates <- apple_intraday_prices[where_dates, 1]
intra_dates <- unname(unlist(Map(intra_dates, f = function(x) subset_string(x, 1, 9))))
# Get the dates (in date format)
Sys.setlocale("LC_TIME", "C") # recognize month's names in english
intra_dates <- as.Date(intra_dates, format = "%d%b%Y") # format example: 14JUN2017
intra_prices <- apple_intraday_prices[-where_dates, 2]
# Actualizing the position where a new day starts 
# in the vector intra_prices once removed the dates
where_dates <- where_dates - 0:(length(where_dates) - 1)

# months range of intraday prices
library("lubridate")
fyear <- subset_string(as.character(year(intra_dates[1])), 3, 4) # get the first year
fmonth <- toupper(month.abb[month(intra_dates[1])]) # get the first month
fdate <- paste(fmonth, fyear, sep = "") # bind first year with first month

# Loading data of option expiring between 2010 and 2017
# OI, daily mean prices, Strike prices, pinning degrees, etc
load(file = "apple_data_2010_2018.Rdata")
apple_data <- apple_data_2010_2018
# take only the options which expire between max(fmonth.fyear, JAN2010) + 1 and DEC2017
# Warning!! I'm assuming the intraday data goes up to 14/12/2017
if(as.numeric(fyear) >= 10){ 
  apple_data <- apple_data[(which(names(apple_data) == fdate) + 1):96]
}

# keeping only the most pinned options per month
most_pinned <- function (formated_list){
  
  most_pinned <- list()
  # delta <- 1 
  
  for(i in 1:length(formated_list)){
  
    # Shorthand for data  and its nrow, ncol
    df <- formated_list[[i]]
    nr <- nrow(df)
    nc <- ncol(df)
    # Shorthand for pinning degrees
    pinning <- attributes(df)$pinning
  
    # Taking what it is important...
    #... most pinned OI evolution and underlyinng stock's prices... 
    who_pinned <- which.min(pinning) # TODO: compute the pinning degrees with the 
                                     # intraday prices instead of the daily prices 
    nona <- min(which(!is.na(df[, 1 + who_pinned]))) # first non-NA indeex of the OI evolution of the most pinned
    most_pinned[[i]] <- df[nona:nr, c(1 + who_pinned, nc)] # OI evolution and prices of the most pinned
  
    #... time (based on daily prices) fitted into [0, 1]
    t <- attributes(df)$time[nona:nr]
    t <- (t - t[1]) / (t[length(t)] - t[1])
    most_pinned[[i]] <- cbind(t, most_pinned[[i]])
  
    # Attributes of the most pinned option (strike price, maturity and pinning degree)
    attr(most_pinned[[i]], "original_strike") <- attributes(df)$strike[who_pinned] # strike price before scaling
    attr(most_pinned[[i]], "maturity_date") <- df[nr, 1] # maturity date
    attr(most_pinned[[i]], "first_date") <- max(df[1, 1], intra_dates[1]) # oppening date
    attr(most_pinned[[i]], "pinning") <- attributes(df)$pinning[who_pinned] # pinning degree
    
    # option's timespan
    first_date <- attributes(most_pinned[[i]])$first_date
    last_date <- attributes(most_pinned[[i]])$maturity_date
    
    # Attribute (intraday prices and frequency)
    which_first <- which(intra_dates == first_date)
    if(length(which_first) == 0) which_first <- which(intra_dates == (first_date + 1))
    which_last <- which(intra_dates == last_date)
    if(length(which_last) == 0) which_last <- which(intra_dates == (last_date - 1))
    
    subset_intra_prices <- intra_prices[where_dates[which_first]:where_dates[which_last]]
    
    attr(most_pinned[[i]], "intra_prices") <- subset_intra_prices # intraday prices
    attr(most_pinned[[i]], "frec_intra_prices") <- "5min" # frequency
    
    # Scalling (divide by the strike price) prices (both daily and intraday) to make them comparable
    most_pinned[[i]][, 3] <- most_pinned[[i]][, 3] / attributes(df)$strike[who_pinned]
    attributes(most_pinned[[i]])$intra_prices <- attributes(most_pinned[[i]])$intra_prices / 
                                                 attributes(df)$strike[who_pinned]  
  
    # Storing the pinning_degree
    # pinning_degree <- c(pinning_degree, pinning[who_pinned])
  
    # Naming data
    colnames(most_pinned[[i]]) <- c("time", "OI", "daily_average_price")
    rownames(most_pinned[[i]]) <- NULL
    
    # Updating delta
    # delta <- min(delta, min(t[2:length(t)] - t[1:(length(t) - 1)]))
  
  }
  # Naming most_pinned list
  names(most_pinned) <- names(formated_list)
  
  # putting all the pinning degrees in a single vector 
  attr(most_pinned, "pinning_degree") <- unname(unlist(Map(most_pinned, f = function(x) attributes(x)$pinning)))

  # Creaating common time
  # for(i in 1:length(most_pinned)){
  
  #  t <- most_pinned[[i]][, 1]
  #  t <- t * delta / min(t[2:length(t)] - t[1:(length(t) - 1)])
  #  t <- t + (1 - t[length(t)])
  
  #  attr(most_pinned[[i]], "ctime") <- t
  
  #}
  
  return(most_pinned)

}
# Getting the most pinned options
most_p <- most_pinned(apple_data)

# Getting the LIBOR
# LIBOR <- read.csv(file = "LIBOR/USD1MTD156N.csv", header = TRUE, stringsAsFactors = FALSE)
# LIBOR[, 1] <- as.Date(LIBOR[, 1], format = "%Y-%m-%d")
# while (any(LIBOR[, 2] == ".")) {
#  LIBOR[which(LIBOR[, 2] == "."), 2] <- LIBOR[which(LIBOR[, 2] == ".") - 1, 2]    
# }
# LIBOR[, 2] <- as.numeric(LIBOR[, 2]) / 100

# Getting the T-Bill interest rate
TBill <- read.csv(file = "daily_T-Bill.csv", header = TRUE, stringsAsFactors = FALSE)
TBill[, 1] <- as.Date(TBill[, 1], format = "%m/%d/%Y")
# 13 weeks treasury bill rate
W13_TBill <- TBill[, c(1, 6)]
W13_TBill[W13_TBill[, 2] == "N/A ", 2] <- W13_TBill[which(W13_TBill[, 2] == "N/A ") - 1, 2] 
W13_TBill[, 2] <- as.numeric(W13_TBill[, 2]) / 100

# 52 (1Year) weeks treasury bill rate
W52_TBill <- TBill[, c(1, 6)]
W52_TBill[W52_TBill[, 2] == "N/A ", 2] <- W52_TBill[which(W52_TBill[, 2] == "N/A ") - 1, 2] 
W52_TBill[, 2] <- as.numeric(W52_TBill[, 2]) / 100
# Setting up the discount rate (it can be either LIBOR or TBill)
DiscountRate13W <- W13_TBill
# DiscountRate13W[, 1] <- as.Date(as.character(DiscountRate13W[, 1]), "%y-%m-%d") # If TBill were selected
DiscountRate52W <- W52_TBill
# DiscountRate52W[, 1] <- as.Date(as.character(DiscountRate52W[, 1]), "%y-%m-%d") # If TBill were selected

###### ######
###### Plot most pinned ######
for (i in 1:length(most_p)) {
  
  prices <- attributes(most_p[[i]])$intra_prices
  plot(1:length(prices)/length(prices), type = "l", prices, bty = 'n', 
       ylab = 'Intraday price (5min)', xlab = 'Time', xlim = c(0, 1), 
       ylim = c(min(prices, 1), max(prices, 1)))
  
  title(main = names(most_p)[i],
        line = 1.5, cex.main = 0.8)
  lines(c(0, 1), c(1, 1), lty = 2, col = "red") 
  mtext(side = 3, line = -2.8, 
        text = paste("                  Pinning degree:", as.character(round(attributes(most_p[[i]])$pinning, 4))), adj = 0, outer = T, cex = 0.8) 
  
}
###### ######

###### Computing profits (real data) ######
# Loading R functions
source("1-boundary_computation_BB-AmPut_discount.R")
source("1-boundary_computation_GBM-AmPut_discount.R")

# Warming times
Wtime <- seq(0.1, 0.9, by = 0.1)

# Preallocating profits
profitBB <- matrix(nrow = length(Wtime), ncol = length(most_p))
profitBB.low <- matrix(nrow = length(Wtime), ncol = length(most_p))
profitGBM <- matrix(nrow = length(Wtime), ncol = length(most_p))


# performing the comparison
for (k in 1:length(Wtime)) {
  
  # Setting up the warming time
  wtime <- Wtime[k]
  # Computing the profits for each strategy 
  # Note!! Set FigGen1 <- TRUE for saving the images
  FigGen1 <- FALSE
  for(i in 1:length(most_p)){
    print(i)
    # strike price
    strike <- 1 # All strike prices are standarized to 1
  
    # Shorthand data
    df <- most_p[[i]]
    
    # time and intraday prices
    prices <- attributes(df)$intra_prices
    prices <- prices[prices != 0]
    t <- seq(0, 1, l = length(prices)) # WARNING! I'm taking the time partition 
                                       # equally spaciated, without taking into account
                                       # the market inactivity (e.g., weekends)
    
    # Shorthands
    lp <- length(prices) # length intraday prices
    cp <- floor(wtime * lp) # cutpoint
    
    # find out the date associated to cp
    intraday_position <- where_dates[which(intra_dates == attributes(df)$first_date)] + cp
    intraday_date_position <- max(which(where_dates <= intraday_position))
    cp_date <- intra_dates[intraday_date_position]
    
    # drift GBM = 1Y Treasury Bill discount rate
    while (length(which(DiscountRate52W[, 1] == cp_date)) == 0) {
    
        cp_date <- cp_date - 1
    
        }
    
    driftGBM <- DiscountRate52W[which(DiscountRate52W[, 1] == cp_date), 2]
    driftGBM <- max(driftGBM, 0.0001)
    
    # Sigma estimation
    sigmaBB <- SigMl(samp = t(as.matrix(prices[1:cp])),  
                     t = t[1:cp], N1 <- cp)               # BB volatility
    
    sigmaGBM <- optim(par = 0.5, lower = 0.0001, fn = GBM.lik, 
                      X = list(drift = driftGBM, path = prices[1:cp], time = t[1:cp]), 
                      method = "L-BFGS-B") # GBM volatility
    sigmaGBM <- sigmaGBM$par
    
    print(sigmaBB)
    print(sigmaGBM)
    print(driftGBM)
    
    # Logarithmic time for computing the boundary
    tlog <- log(seq(exp(0), exp(1), l = 201))
    
    # Boundary estimation
    boundaryBB <- b_BB_AmPut(sigma = sigmaBB, t = tlog, r = driftGBM, S = strike) # BB
    boundaryGBM <- b_GBM_AmPut(mu = driftGBM, sigma = sigmaGBM, t = tlog, S = strike) # GBM
    
    # Confidence boundaries
    boundaryBB <- b_BB_AmPut_Conf(boundaryBB, S = strike, sigma = sigmaBB,
                                  r = driftGBM, t = tlog, nsamp = cp)
    boundaryBB.low <- boundaryBB$bBB.low
    boundaryBB <- boundaryBB$bBB
    
    # boundary at needed points
    boundaryBB <- splinefun(tlog, boundaryBB)(t) # BB
    boundaryBB.low <- splinefun(tlog, boundaryBB.low)(t) # GBM
    boundaryGBM <- splinefun(tlog, boundaryGBM)(t) # GBM
    
    # proportion 5minutes / 1Year
    # Assuming 250 working days, and 8 hours per day.
    m_y <- 1/(250 * 8 * 12)
    
    # Profit
    idx <- min(which(prices[cp:lp] <= boundaryBB[cp:lp]))
    profitBB[k, i] <- exp(-driftGBM * (idx - 1) * m_y) * 
                      (strike - prices[cp:lp][idx])
    # if(is.na(profitBB[length(profitBB)])) profitBB[length(profitBB)] <- 0 # BB
    
    # Profit
    idx <- min(which(prices[cp:lp] <= boundaryBB.low[cp:lp])) 
    profitBB.low[k, i] <- exp(-driftGBM * (idx - 1) * m_y) * 
                          (strike - prices[cp:lp][idx])
    # if(is.na(profitBB[length(profitBB)])) profitBB[length(profitBB)] <- 0 # BB
    
    idx <- min(which(prices[cp:lp] <= boundaryGBM[cp:lp]))
    profitGBM[k, i] <- exp(-driftGBM * (idx - 1) * m_y) * 
                       (strike - prices[cp:lp][idx])
    # if(is.na(profitGBM[length(profitGBM)])) profitGBM[length(profitGBM)] <- 0 # GBM
    
    # plot
    ylim <- c(min(c(boundaryBB[1], boundaryBB.low[1], boundaryGBM[1], prices)),
              max(c(strike, prices)))
    
    if(FigGen1 == TRUE) pdf(paste("real_data_BB_vs_GBM_", names(most_p)[i], "_", 
                                  "wtime_0*", 
                                  subset_string(as.character(round(wtime, 2)), 3, "end"), 
                                  ".pdf", sep = ""))
    
    plot(0, bty = 'n', pch = '', ylab = 'Price', xlab = 'Time', xlim = c(0, 1), ylim = ylim)
    lines(t[1:cp], prices[1:cp], lty = 1, lwd = 1.5)
    lines(t[cp:lp], prices[cp:lp], lty = 1, lwd = 1)
    lines(t[cp:lp], boundaryBB[cp:lp], lty = 1, lwd = 1, col = "red")
    lines(t[cp:lp], boundaryBB.low[cp:lp], lty = 2, lwd = 1, col = "red")
    lines(t[cp:lp], boundaryGBM[cp:lp], lty = 1, lwd = 1, col = "blue")
    lines(rep(t[cp], 2), ylim, lty = 2)
    lines(c(0, 1), rep(strike, 2), lty = 2)
    title(main = names(most_p)[i], line = 3.5, cex.main = 0.8) 
    title(main = paste("Sigma estimation: ", as.character(round(sigmaBB, 3)), sep = ""),
          line = 2.5, cex.main = 0.8, col.main = "red")
    title(main = paste("Sigma estimation: ", as.character(round(sigmaGBM, 3)), sep = ""),
          line = 1.5, cex.main = 0.8, col.main = "blue")
    title(main = paste("Drift estimation: ", as.character(round(driftGBM, 3)), sep = ""),
          line = 0.5, cex.main = 0.8, col.main = "blue")
    mtext(side = 3, line = -1.8, 
          text = paste("                  Pinning degree:", round(attributes(df)$pinning, 4)),
          adj = 0, outer = T, cex = 0.8) 
    mtext(side = 3, line = -2.8, 
          text = paste("                  Warming time:", round(wtime, 2)), 
          adj = 0, outer = T, cex = 0.8) 
    mtext(side = 3, line = -3.8, 
          text = paste("                  Sample size :", length(prices)), 
          adj = 0, outer = T, cex = 0.8) 
    
    if (FigGen1 == TRUE) dev.off()
    
  }

}

profitBB[is.na(profitBB)] <- 0
profitBB.low[is.na(profitBB.low)] <- 0
profitGBM[is.na(profitGBM)] <- 0
# ##### ######
###### Comparing (real data) ######
# rowVars function (compue the variance per row)
rowVars <- function(X) rowMeans((X - rowMeans(X)) ^ 2)

# Getting the pinning degrees
pinning_degree <- attributes(most_p)$pinning_degree
  
qvector <- seq(0.025, 1, by = 0.025)

# preallocating mean and var profits
mean_profitBB <- matrix(nrow = length(Wtime), ncol = length(qvector))
mean_profitBB.low <- matrix(nrow = length(Wtime), ncol = length(qvector))
mean_profitGBM <- matrix(nrow = length(Wtime), ncol = length(qvector))

var_profitBB <- matrix(nrow = length(Wtime), ncol = length(qvector))
var_profitBB.low <- matrix(nrow = length(Wtime), ncol = length(qvector))
var_profitGBM <- matrix(nrow = length(Wtime), ncol = length(qvector))

# Making the comparison
for(q in 1:length(qvector)){
# for(q in 0.45){
  print(qvector[q])
  # Declaring pinning threeshold and pinned options
  pinning_threeshold <- quantile(pinning_degree, qvector[q]) # pinning threeshold
  pinned <- which(pinning_degree <= pinning_threeshold) # options whose ponning degree is lower than the pinning threeshold

  # compute mean profit
  mean_profitBB[, q] <- rowMeans(profitBB[, pinned])
  mean_profitBB.low[, q] <- rowMeans(profitBB.low[, pinned])
  mean_profitGBM[, q] <- rowMeans(profitGBM[, pinned])

  # compute var profit
  var_profitBB[, q] <- rowVars(profitBB[, pinned])
  var_profitBB.low[, q] <- rowVars(profitBB.low[, pinned])
  var_profitGBM[, q] <- rowVars(profitGBM[, pinned])

}

mean_profitBB_wt_ag <- colMeans(mean_profitBB)
mean_profitBB.low_wt_ag <- colMeans(mean_profitBB.low)
mean_profitGBM_wt_ag <- colMeans(mean_profitGBM)
  
# Note!! Set FigGen2 <- TRUE for saving the images
# FigGen2 <- TRUE
# Mean profit plot (aggregated wtime)
if(FigGen2 == TRUE) pdf("real_data_mean_profit_BB_vs_GBM_wtime_aggregated.pdf")

ylim <- c(min(0, mean_profitBB_wt_ag, mean_profitGBM_wt_ag), 
          max(mean_profitBB_wt_ag, mean_profitGBM_wt_ag))

plot(0, bty = 'n', pch = '', ylab = 'Profit', xlab = "", xlim = c(0, 1), ylim = ylim)
lines(qvector, mean_profitBB_wt_ag, lty = 1, col = "red")
lines(qvector, mean_profitGBM_wt_ag, lty = 1, col = "blue")
title(main = "Mean profit", line = 1.8, cex.main = 0.8)
legend("bottomright", lty = c(1, 1), lwd = c(1.5, 1.5), 
        col = c("red", "blue"), legend = c("BB", "GBM"), bty = "n")
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
     labels = paste(rep("(", 4), round(quantile(pinning_degree, c(0, 0.2, 0.4, 0.6, 0.8, 1)), 4), rep(")", 4), sep = ""), 
     tick = FALSE, line = 1, pos = NA, lty = 3, lwd = 0.01)
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
     labels = paste(rep("[", 4), floor((length(most_p) * c(0, 0.2, 0.4, 0.6, 0.8, 1))), rep("]", 4), sep = ""), 
     tick = FALSE, line = 2, pos = NA, lty = 3, lwd = 0.01)
title(xlab = "quantile (pinning degree) [number of samples]", line = 4)
  
if (FigGen2 == TRUE) dev.off()
  
# Relative mean profit plot
  
if(FigGen2 == TRUE) pdf("real_data_relative_mean_profit_BB_vs_GBM_wtime_aggregated.pdf")
  
plot(qvector, (mean_profitBB_wt_ag - mean_profitGBM_wt_ag) / mean_profitGBM_wt_ag,
     bty = "n", type = 'l', pch = '', ylab = 'Relative mean profit', xlab = "", xlim = c(0, 1),
     ylim = c(min(0, (mean_profitBB_wt_ag - mean_profitGBM_wt_ag) / mean_profitGBM_wt_ag),
              max((mean_profitBB_wt_ag - mean_profitGBM_wt_ag) / mean_profitGBM_wt_ag)))
lines(c(0, 1), c(0, 0), lty = 2)
# title(main = "Relative mean profit", line = 0.5, cex.main = 0.8)
# title(main = "(MeanBB - MeanGBM) / MeanGBM", line = 0.5, cex.main = 0.8)

# axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
#     labels = paste(rep("(", 4), round(quantile(pinning_degree, c(0, 0.2, 0.4, 0.6, 0.8, 1)), 4), rep(")", 4), sep = ""), 
#     tick = FALSE, line = 1, pos = NA, lty = 3, lwd = 0.01)
# axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
#     labels = paste(rep("[", 4), floor((length(most_p) * c(0, 0.2, 0.4, 0.6, 0.8, 1))), rep("]", 4), sep = ""), 
#     tick = FALSE, line = 2, pos = NA, lty = 3, lwd = 0.01)
title(xlab = "quantile pinning deviance", line = 2.5)

if (FigGen2 == TRUE) dev.off()
  
# Var profit plot

var_profitBB_wt_ag <- colMeans(var_profitBB)
var_profitGBM_wt_ag <- colMeans(var_profitGBM)
  
if(FigGen2 == TRUE) pdf("real_data_var_profit_BB_vs_GBM_wtime_aggregated.pdf")
  
ylim <- c(min(var_profitBB_wt_ag, var_profitGBM_wt_ag, na.rm = TRUE), 
          max(var_profitBB_wt_ag, var_profitGBM_wt_ag, na.rm = TRUE))

plot(0, bty = 'n', pch = '', ylab = 'Profit', xlab = "", xlim = c(0, 1), ylim = ylim)
lines(qvector, var_profitBB_wt_ag, lty = 1, col = "red")
lines(qvector, var_profitGBM_wt_ag, lty = 1, col = "blue")
title(main = "Var profit", line = 1.8, cex.main = 0.8)
legend("bottomright", lty = c(1, 1), lwd = c(1.5, 1.5), 
       col = c("red", "blue"), legend = c("BB", "GBM"), bty = "n")
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
     labels = paste(rep("(", 4), round(quantile(pinning_degree, c(0, 0.2, 0.4, 0.6, 0.8, 1)), 4), rep(")", 4), sep = ""), 
     tick = FALSE, line = 1, pos = NA, lty = 3, lwd = 0.01)
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
     labels = paste(rep("[", 4), floor((length(most_p) * c(0, 0.2, 0.4, 0.6, 0.8, 1))), rep("]", 4), sep = ""), 
     tick = FALSE, line = 2, pos = NA, lty = 3, lwd = 0.01)
title(xlab = "quantile (pinning degree) [number of samples]", line = 4)
  
if (FigGen2 == TRUE) dev.off()
  
# Relative Var profit plot
  
if(FigGen2 == TRUE) pdf("real_data_relative_var_profit_BB_vs_GBM_wtime_aggregated.pdf")
  
plot(qvector, (var_profitGBM_wt_ag - var_profitBB_wt_ag) / var_profitGBM_wt_ag,
     bty = "n", type = 'l', pch = '', ylab = 'Profit', xlab = "", xlim = c(0, 1))
lines(c(0, 1), c(0, 0), lty = 2)
title(main = "Relative var profit", line = 1.8, cex.main = 0.8)
title(main = "(VarGBM - VarBB) / VarGBM", line = 0.5, cex.main = 0.8)
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
     labels = paste(rep("(", 4), round(quantile(pinning_degree, c(0, 0.2, 0.4, 0.6, 0.8, 1)), 4), rep(")", 4), sep = ""), 
     tick = FALSE, line = 1, pos = NA, lty = 3, lwd = 0.01)
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
     labels = paste(rep("[", 4), floor((length(most_p) * c(0, 0.2, 0.4, 0.6, 0.8, 1))), rep("]", 4), sep = ""), 
     tick = FALSE, line = 2, pos = NA, lty = 3, lwd = 0.01)
title(xlab = "quantile (pinning degree) [number of samples]", line = 4)

if (FigGen2 == TRUE) dev.off()
###### ######

###### Computing profits (simulated data) ######
library(Rcpp)
library(RcppArmadillo)
Sys.setenv("PKG_CXXFLAGS" = "-std=c++11") # Needs C++11
# Cpp script
Rcpp::sourceCpp('1-boundary_computation_BB-Identity.cpp')
# Loading R functions
source("1-boundary_computation_BB-AmPut_discount.R")
source("1-boundary_computation_GBM-AmPut_discount.R")
source("0-simulations.R")
# Simulating data (each one with 501 observations)
# set.seed(123)
# top <- 0.7
time_partition <- seq(0, 1, by = 1/500)
simulated_pinning_degree <- seq(0, top, by = 0.001)
# simulated_data <- sapply(simulated_pinning_degree, 
#                          function(x) {
#                             print(x)
#                             pinning_point <- 1 + sample(c(-1, 1), 1) * x
#                             rBB(n = 1, a = rnorm(1, mean = pinning_point, sd = 0.1),
#                                 b = pinning_point, t = time_partition,
#                                 sigma = runif(1, min = 0.01, max = 0.1 + 0.3 * (top - x)))
#                          })
# simulated_data <- t(simulated_data)
# save(simulated_data, file = "simulated_data.RData")
load(file = "simulated_data.RData")
# Setting up the warming time
wtime <- 1/4
# creating proffits
profitBB_simulated <- c()
profitGBM_simulated <- c()
# Computing the profits for each strategy 
# Note!! Set FigGen1 <- TRUE for saving the images
FigGen1 <- FALSE
for(i in 1:nrow(simulated_data)){
  print(i)
  # strike price
  strike <- 1 # All strike prices are standarized to 1
  
  # time and intraday prices
  prices <- simulated_data[i, ]
  t <- seq(0, 1, l = length(prices)) # WARNING! I'm taking the time partition 
  # equally spaciated, without taking into account
  # the market inactivity (e.g., weekends)
  
  # Shorthands
  lp <- length(prices) # length intraday prices
  cp <- floor(wtime * lp) # cutpoint
  
  # Sigma estimation
  sigmaBB <- SigMl(samp = t(as.matrix(prices)),  
                   t = t, N1 <- cp)               # BB volatility
  
  sigmaGBM <- unname(sigma_GBM(t[1:cp], prices[1:cp])[2]) # GBM volatility
  driftGBM <- max(sigma_GBM(t[1:cp], prices[1:cp])[1], 0.01) # GBM drift
  
  print(sigmaBB)
  print(sigmaGBM)
  print(driftGBM)
  
  # Logarithmic time for computing the boundary
  tlog <- log(seq(exp(0), exp(1), l = 201))
  
  # Boundary estimation
  boundaryBB <- b_BB_AmPut(sigma = sigmaBB, t = tlog, r = driftGBM, y = strike) # BB
  boundaryGBM <- b_GBM_AmPut(mu = driftGBM, sigma = sigmaGBM, t = tlog, S = strike) # GBM
  
  # boundary at needed points
  boundaryBB <- splinefun(tlog, boundaryBB)(t) # BB
  boundaryGBM <- splinefun(tlog, boundaryGBM)(t) # GBM
  
  # Profit
  profitBB_simulated <- c(profitBB_simulated, 
                          strike - prices[cp:lp][min(which(prices[cp:lp] <= boundaryBB[cp:lp]))])
  if(is.na(profitBB_simulated[length(profitBB_simulated)])){
    profitBB_simulated[length(profitBB_simulated)] <- 0 # BB
  } 
  
  profitGBM_simulated <- c(profitGBM_simulated, 
                           strike - prices[cp:lp][min(which(prices[cp:lp] <= boundaryGBM[cp:lp]))])
  if(is.na(profitGBM_simulated[length(profitGBM_simulated)])){
    profitGBM_simulated[length(profitGBM_simulated)] <- 0 # GBM
  }
  
  # plot
  ylim <- c(min(c(boundaryBB[1], boundaryGBM[1], prices)),
            max(c(strike, prices)))
  
  if(FigGen1 == TRUE) pdf(paste("simulated_data_BB_vs_GBM_", names(most_p)[i], "_", 
                                "wtime_0*", 
                                subset_string(as.character(round(wtime, 2)), 3, "end"), 
                                ".pdf", sep = ""))
  
  plot(0, bty = 'n', pch = '', ylab = 'Price', xlab = 'Time', xlim = c(0, 1), ylim = ylim)
  lines(t[1:cp], prices[1:cp], lty = 1, lwd = 1.5)
  lines(t[cp:lp], prices[cp:lp], lty = 1, lwd = 1)
  lines(t[cp:lp], boundaryBB[cp:lp], lty = 1, lwd = 1, col = "red")
  lines(t[cp:lp], boundaryGBM[cp:lp], lty = 1, lwd = 1, col = "blue")
  lines(rep(t[cp], 2), ylim, lty = 2)
  lines(c(0, 1), rep(strike, 2), lty = 2)
  title(main = paste("Sigma estimation: ", as.character(round(sigmaBB, 3)), sep = ""),
        line = 2.5, cex.main = 0.8, col.main = "red")
  title(main = paste("Sigma estimation: ", as.character(round(sigmaGBM, 3)), sep = ""),
        line = 1.5, cex.main = 0.8, col.main = "blue")
  title(main = paste("Drift estimation: ", as.character(round(driftGBM, 3)), sep = ""),
        line = 0.5, cex.main = 0.8, col.main = "blue")
  mtext(side = 3, line = -1.8, 
        text = paste("                  Pinning degree:", round(simulated_pinning_degree[i], 4)),
        adj = 0, outer = T, cex = 0.8) 
  mtext(side = 3, line = -2.8, 
        text = paste("                  Warming time:", round(wtime, 2)), 
        adj = 0, outer = T, cex = 0.8) 
  mtext(side = 3, line = -3.8, 
        text = paste("                  Sample size :", length(prices)), 
        adj = 0, outer = T, cex = 0.8) 
  
  if (FigGen1 == TRUE) dev.off()
  
}
###### ######
###### Comparing (simulated data) ######
# Allocating mean and var profits
mean_profitBB_simulated <- c()
mean_profitGBM_simulated <- c()

var_profitBB_simulated <- c()
var_profitGBM_simulated <- c()

qvector <- seq(0.01, 1, by = 0.01)
# Making the comparison
for(q in qvector){
  # for(q in 0.45){
  print(q)
  # Declaring pinning threeshold and pinned options
  pinning_threeshold <- quantile(simulated_pinning_degree, q) # pinning threeshold
  pinned <- which(simulated_pinning_degree <= pinning_threeshold) # options whose ponning degree is lower than the pinning threeshold
  
  # compute mean profit
  mean_profitBB_simulated <- c(mean_profitBB_simulated, mean(profitBB_simulated[pinned]))
  mean_profitGBM_simulated <- c(mean_profitGBM_simulated, mean(profitGBM_simulated[pinned]))
  
  # compute var profit
  var_profitBB_simulated <- c(var_profitBB_simulated, var(profitBB_simulated[pinned]))
  var_profitGBM_simulated <- c(var_profitGBM_simulated, var(profitGBM_simulated[pinned]))
  
}

# Note!! Set FigGen2 <- TRUE for saving the images
# FigGen2 <- TRUE
# Mean profit plot
if(FigGen2 == TRUE) pdf(paste("simulated_data_mean_profit_BB_vs_GBM_wtime_",
                              subset_string(as.character(round(wtime, 2)), 3, "end"),
                              ".pdf", sep = ""))

ylim <- c(min(mean_profitBB_simulated, mean_profitGBM_simulated), 
          max(mean_profitBB_simulated, mean_profitGBM_simulated))

plot(0, bty = 'n', pch = '', ylab = 'Profit', xlab = "", xlim = c(0, 1), ylim = ylim)
lines(qvector, mean_profitBB_simulated, lty = 1, col = "red")
lines(qvector, mean_profitGBM_simulated, lty = 1, col = "blue")
title(main = "Mean profit (Simulated data)", line = 1.8, cex.main = 0.8)
legend("bottomright", lty = c(1, 1), lwd = c(1.5, 1.5), 
       col = c("red", "blue"), legend = c("BB", "GBM"), bty = "n")
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
     labels = paste(rep("(", 4), round(quantile(simulated_pinning_degree, c(0, 0.2, 0.4, 0.6, 0.8, 1)), 4), rep(")", 4), sep = ""), 
     tick = FALSE, line = 1, pos = NA, lty = 3, lwd = 0.01)
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
     labels = paste(rep("[", 4), floor((nrow(simulated_data) * c(0, 0.2, 0.4, 0.6, 0.8, 1))), rep("]", 4), sep = ""), 
     tick = FALSE, line = 2, pos = NA, lty = 3, lwd = 0.01)
mtext(side = 3, line = -2.4, 
      text = paste("                  Warming time:", round(wtime, 2)), 
      adj = 0, outer = T, cex = 0.8) 
title(xlab = "quantile (pinning degree) [number of samples]", line = 4)

if (FigGen2 == TRUE) dev.off()

# Relative mean profit plot

if(FigGen2 == TRUE) pdf(paste("simulated_data_relative_mean_profit_BB_vs_GBM_wtime_",
                              subset_string(as.character(round(wtime, 2)), 3, "end"),
                              ".pdf", sep = ""))

plot(qvector, (mean_profitBB_simulated - mean_profitGBM_simulated) / mean_profitGBM_simulated,
     bty = "n", type = 'l', pch = '', ylab = 'Profit', xlab = "", xlim = c(0, 1))
lines(c(0, 1), c(0, 0), lty = 2)
title(main = "Relative mean profit (Simulated data)", line = 1.8, cex.main = 0.8)
title(main = "(MeanBB - MeanGBM) / MeanGBM", line = 0.5, cex.main = 0.8)
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
     labels = paste(rep("(", 4), round(quantile(simulated_pinning_degree, c(0, 0.2, 0.4, 0.6, 0.8, 1)), 4), rep(")", 4), sep = ""), 
     tick = FALSE, line = 1, pos = NA, lty = 3, lwd = 0.01)
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
     labels = paste(rep("[", 4), floor((nrow(simulated_data) * c(0, 0.2, 0.4, 0.6, 0.8, 1))), rep("]", 4), sep = ""), 
     tick = FALSE, line = 2, pos = NA, lty = 3, lwd = 0.01)
mtext(side = 3, line = -2.4, 
      text = paste("                  Warming time:", round(wtime, 2)), 
      adj = 0, outer = T, cex = 0.8) 
title(xlab = "quantile (pinning degree) [number of samples]", line = 4)

if (FigGen2 == TRUE) dev.off()

# Var profit plot

if(FigGen2 == TRUE) pdf(paste("simulated_data_var_profit_BB_vs_GBM_wtime_",
                              subset_string(as.character(round(wtime, 2)), 3, "end"),
                              ".pdf", sep = ""))

ylim <- c(min(var_profitBB_simulated, var_profitGBM_simulated, na.rm = TRUE), 
          max(var_profitBB_simulated, var_profitGBM_simulated, na.rm = TRUE))

plot(0, bty = 'n', pch = '', ylab = 'Profit', xlab = "", xlim = c(0, 1), ylim = ylim)
lines(qvector, var_profitBB_simulated, lty = 1, col = "red")
lines(qvector, var_profitGBM_simulated, lty = 1, col = "blue")
title(main = "Var profit (Simulated data)", line = 1.8, cex.main = 0.8)
legend("bottomright", lty = c(1, 1), lwd = c(1.5, 1.5), 
       col = c("red", "blue"), legend = c("BB", "GBM"), bty = "n")
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
     labels = paste(rep("(", 4), round(quantile(simulated_pinning_degree, c(0, 0.2, 0.4, 0.6, 0.8, 1)), 4), rep(")", 4), sep = ""), 
     tick = FALSE, line = 1, pos = NA, lty = 3, lwd = 0.01)
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
     labels = paste(rep("[", 4), floor((nrow(simulated_data) * c(0, 0.2, 0.4, 0.6, 0.8, 1))), rep("]", 4), sep = ""), 
     tick = FALSE, line = 2, pos = NA, lty = 3, lwd = 0.01)
mtext(side = 3, line = -2.4, 
      text = paste("                  Warming time:", round(wtime, 2)), 
      adj = 0, outer = T, cex = 0.8) 
title(xlab = "quantile (pinning degree) [number of samples]", line = 4)

if (FigGen2 == TRUE) dev.off()

# Relative Var profit plot

if(FigGen2 == TRUE) pdf(paste("simulated_data_relative_var_profit_BB_vs_GBM_wtime_",
                              subset_string(as.character(round(wtime, 2)), 3, "end"),
                              ".pdf", sep = ""))

plot(qvector, (var_profitGBM_simulated - var_profitBB_simulated) / var_profitGBM_simulated,
     bty = "n", type = 'l', pch = '', ylab = 'Profit', xlab = "", xlim = c(0, 1))
lines(c(0, 1), c(0, 0), lty = 2)
title(main = "Relative var profit (Simulated data)", line = 1.8, cex.main = 0.8)
title(main = "(VarGBM - VarBB) / VarGBM", line = 0.5, cex.main = 0.8)
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
     labels = paste(rep("(", 4), round(quantile(simulated_pinning_degree, c(0, 0.2, 0.4, 0.6, 0.8, 1)), 4), rep(")", 4), sep = ""), 
     tick = FALSE, line = 1, pos = NA, lty = 3, lwd = 0.01)
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
     labels = paste(rep("[", 4), floor((nrow(simulated_data) * c(0, 0.2, 0.4, 0.6, 0.8, 1))), rep("]", 4), sep = ""), 
     tick = FALSE, line = 2, pos = NA, lty = 3, lwd = 0.01)
mtext(side = 3, line = -2.4, 
      text = paste("                  Warming time:", round(wtime, 2)), 
      adj = 0, outer = T, cex = 0.8) 
title(xlab = "quantile (pinning degree) [number of samples]", line = 4)

if (FigGen2 == TRUE) dev.off()
###### ######