# Estimating population trends for US West Coast humpback stocks, 2006-2024
# Aug 2026, Jeff Moore

# ANALYSIS SUMMARY
# The code below uses the following data/inputs:
#   USWC time series of N estimates for USWC humpbacks 1989-2024 (Calambokidis and Barlow 2020, updated in 2025)
#   Curtis et al. (2026): estimate of N for Central America stock/DPS, and CASM trend estimate for 2004-2006 to 2022-2024
# A Monte Carlo approach is used
# USWC trend = MC distribution of mean lambdas (geometric means) from distributions of N(t)
# CASM = back-calculate distributions for N(t) from distns for N(2024) and trend, based on Curtis et al (2026)
# MX:  N(t).MX = N(t).CASM - N(t).USWC, then find geometric means of lambda(t)
# This approach initially results in some impossible outcomes, where N.CASM(2006) > N.USWC(2006)
# Address this by censoring all MC samples with these impossible outcomes and then re-summarizing the distributions for any derived parms

##################################
# packages
library(readxl)

###################################
# functions

logNormParms <- function(mean, sd=NA, cv=NA){
  if(is.na(cv)==FALSE) varlog <- log(cv * cv + 1)
  if(is.na(sd)==FALSE) varlog <- log((sd/mean)^2 + 1)
  mulog <- log(mean) - varlog/2
  sdlog <- sqrt(varlog)
  out <- list(mulog=mulog, sdlog=sdlog)
  out
}

rlognorm <- function(nsims, mean, sd=NA, cv=NA){
  if(is.na(cv)==FALSE) varlog <- log(cv * cv + 1)
  if(is.na(sd)==FALSE) varlog <- log((sd/mean)^2 + 1)
  mulog <- log(mean) - varlog/2
  sdlog <- sqrt(varlog)
  data <- rlnorm(nsims,mulog,sdlog)
  data
}

###################################
# ANALYSIS

set.seed(23)

# CRC abundance estimates (Chao) 1989:2024
year <- c(1989:2024)
nchao <- c(506, 517,                                                       # 1989-1990
           636, 766, 800, 830, 834, 922, 974, 1068, 1263, 1166,            # 1991-2000
           1220, 1277, 1566, 1618, 1878, 1985, 1826, 2123, 2131, 2324,     # 2001-2010
           2495, 2340, 2394, 2270, 2718, 4091, 4444, 5261, 5462, 5226,     # 2011-2020
           5558, 5577, 6208, 7021)                                         # 2021-2024
se.n <- c(40.5798, 36.5816,  
          42.3818, 45.4783, 39.4265, 37.7444, 36.5239, 47.8338, 52.9319, 53.6551, 68.1018, 65.4120,  
          70.4701,  82.4916, 106.3093, 102.9192, 119.9618, 132.6319, 135.3533, 155.1913, 148.0269, 145.0735, 
          144.5661, 123.9231, 120.1721, 105.6575, 124.2316, 182.3272, 171.5744, 174.6290, 168.3726, 162.9086, 
          179.5160, 177.7710, 179.6380, 181.4879)

nsims = 5000  # for monte carlo sampling

# distn of estimates for USWC
n.uswc <- matrix(NA,length(year),nsims)
row.names(n.uswc) <- year
for(i in 1:length(year)){
  n.uswc[i,] <- rlognorm(nsims=nsims,mean=nchao[i],sd=se.n[i])
}
# select 2006 - 2024 only (years since SPLASH surveys)
n.uswc.06.24 <- n.uswc[18:36,]
n.uswc.2.5 <- apply(n.uswc.06.24,1,quantile,0.025)
n.uswc.97.5 <- apply(n.uswc.06.24,1,quantile,0.975)

# distn of estimates for CASM
yr <- 2006:2024
tsl <- length(yr)  # time series length
casm.24 <- rlognorm(nsims,1453,cv=.108)  # from Curtis
trend <- rnorm(nsims,-0.004,0.028)  # from Curtis
n.casm <- matrix(data=NA,nrow=tsl,ncol=nsims)
rownames(n.casm) <- yr
n.casm[tsl,] <- casm.24
for(i in (tsl-1):1){
  n.casm[i,] <- n.casm[i+1,]/(1+trend)
}
n.casm.mean <- rowMeans(n.casm)
n.casm.median <- apply(n.casm,1,median)
n.casm.2.5 <- apply(n.casm,1,quantile,0.025)
n.casm.97.5 <- apply(n.casm,1,quantile,0.975)

# distn of estimates for MX = USWC - CASM
n.mx.06.24 <- n.uswc.06.24 - n.casm
n.mx.06.24.mean <- rowMeans(n.mx.06.24)
n.mx.06.24.upper <- apply(n.mx.06.24,1,quantile,.975)
n.mx.06.24.lower <- apply(n.mx.06.24,1,quantile,.025)

# eliminate impossible iterations (e.g., where CASM abundance greater than USWC)
idxMatrix.casm.exceeds.uswc <- n.casm > n.uswc.06.24
idxVec.casm.exceeds.uswc <- colSums(idxMatrix.casm.exceeds.uswc)>0
n.uswc.new <- n.uswc.06.24[,idxVec.casm.exceeds.uswc==FALSE]
n.casm.new <- n.casm[,idxVec.casm.exceeds.uswc==FALSE]
n.mx.new <- n.mx.06.24[,idxVec.casm.exceeds.uswc==FALSE]

# new distributions for N, 2006-2024, for USWC, CASM and MX
n.uswc.quants <- cbind(apply(n.uswc.new,1,quantile,c(0.025,.5,.975)))
n.casm.quants <- cbind(apply(n.casm.new,1,quantile,c(0.025,.5,.975)))
n.mx.quants <- cbind(apply(n.mx.new,1,quantile,c(0.025,.5,.975)))

# re-calcuate trends (after removing impossible iterations) 

# Approach 1 (don't use) -- fitting linear model to the log of the median abundance in each year
#newtrend.uswc <- lm(log(n.uswc.quants[2,])~yr); summary(newtrend.uswc)  # 7.9%
#newtrend.casm <- lm(log(n.casm.quants[2,])~yr); summary(newtrend.casm)  # 0.9%
#newtrend.mx <- lm(log(n.mx.quants[2,])~yr); summary(newtrend.mx) # 12.6%  # trend for entire time series
#newtrend.mx.1 <- lm(log(n.mx.quants[2,1:10])~yr[1:10]); summary(newtrend.mx.1) # 6.2%  # trend from 2006 - 2015, before large recruitment pulse from 2015-2016
#newtrend.mx.2 <- lm(log(n.mx.quants[2,11:19])~yr[11:19]); summary(newtrend.mx.2) # 7%  # trend frm 2016 - 2024, after large recruitment pulse from 2015-2016

# Approach 2 (use) -- find geometric mean of each iteration and then summarize that distribution
   # Do for entire time series but also two time periods: 2006-2015 and 2016-2024
   # This is because of pulse increase between 2015-1016
delta <- 2024-2006  # number of annual changes in N
delta.1 <- 2015-2006
delta.2 <- 2024-2016

# USWC = 7.28% (0.41%), 6.46% - 8.06%.  These are mean (SE) and 95% range
newtrend.uswc.geoMean <- (n.uswc.new[delta+1,]/n.uswc.new[1,])^(1/delta) - 1
newtrend.uswc.summary <- c(mean(newtrend.uswc.geoMean), sd(newtrend.uswc.geoMean), quantile(newtrend.uswc.geoMean,c(0.025, 0.5, 0.975))); newtrend.uswc.summary
# CASM = 1.30% (1.95%), -1.70% - 5.74%
newtrend.casm.geoMean <- (n.casm.new[delta+1,]/n.casm.new[1,])^(1/delta) - 1
newtrend.casm.summary <- c(mean(newtrend.casm.geoMean), sd(newtrend.casm.geoMean), quantile(newtrend.casm.geoMean,c(0.025, 0.5, 0.975))); newtrend.casm.summary
# MX = 12.62% (4.59%), 7.45% - 24.48%
newtrend.mx.geoMean <- (n.mx.new[delta+1,]/n.mx.new[1,])^(1/delta) - 1
newtrend.mx.summary <- c(mean(newtrend.mx.geoMean), sd(newtrend.mx.geoMean), quantile(newtrend.mx.geoMean,c(0.025, 0.5, 0.975))); newtrend.mx.summary
# MX from 2006 - 2015 = 8.85% (7.94%), 1.39% - 29.28%
newtrend.mx.geoMean.1 <- (n.mx.new[delta.1+1,]/n.mx.new[1,])^(1/delta.1) - 1
newtrend.mx.summary.1 <- c(mean(newtrend.mx.geoMean.1), sd(newtrend.mx.geoMean.1), quantile(newtrend.mx.geoMean.1,c(0.025, 0.5, 0.975))); newtrend.mx.summary.1
# MX from 2016 - 2024 = 9.16% (1.39%), 6.60% - 11.93%
newtrend.mx.geoMean.2 <- (n.mx.new[11+delta.2,]/n.mx.new[11,])^(1/delta.2) - 1
newtrend.mx.summary.2 <- c(mean(newtrend.mx.geoMean.2), sd(newtrend.mx.geoMean.2), quantile(newtrend.mx.geoMean.2,c(0.025, 0.5, 0.975))); newtrend.mx.summary.2
# MX from 2006 - 2024 (minus the 2015/16 pulse) = 8.94% (4.34%), 4.23% - 20.36%
newtrend.mx.geoMean.3 <- ((newtrend.mx.geoMean.1+1)^delta.1 * (newtrend.mx.geoMean.2+1)^delta.2) ^ (1/(delta.1 + delta.2))
newtrend.mx.summary.3 <- c(mean(newtrend.mx.geoMean.3), sd(newtrend.mx.geoMean.3), quantile(newtrend.mx.geoMean.3,c(0.025, 0.5, 0.975))); newtrend.mx.summary.3


########################
# PLOT RESULTS

# plot revised trends/estimates, where impossible estimates have been removed
x11(width = 8, height = 6)
par(mar = c(5.1, 7.1, 4.1, 2.1))

plot(yr,n.casm.quants[2,],type="l",col="red",ylim=c(0,8000), xlim=c(2005,2025), xlab="Year", ylab="Population size",bty="n", cex.axis=2, cex.lab=2)
polygon(c(yr, rev(yr)), c(n.casm.quants[3,], rev(n.casm.quants[1,])), col = adjustcolor("pink", alpha.f = 0.3), border = NA)
lines(yr, n.casm.quants[2,],type="l", lwd=2, col="red")
#lines(yr,n.casm.quants[1,], lty=2, col="red")
#lines(yr,n.casm.quants[3,], lty=2, col="red")

polygon(c(yr, rev(yr)), c(n.mx.quants[3,], rev(n.mx.quants[1,])), col = adjustcolor("lightblue", alpha.f = 0.3), border = NA)
lines(yr,n.mx.quants[2,],type="l", lwd=2, col="blue",ylim=c(0,8000))
#lines(yr,n.mx.quants[1,], lty=2, col="blue")
#lines(yr,n.mx.quants[3,], lty=2, col="blue")

polygon(c(yr, rev(yr)), c(n.uswc.quants[3,], rev(n.uswc.quants[1,])), col = adjustcolor("gray", alpha.f = 0.3), border = NA)
lines(yr,n.uswc.quants[2,],type="l", lwd=2, ylim=c(0,8000))
#lines(yr,n.uswc.quants[1,], lty=2)
#lines(yr,n.uswc.quants[3,], lty=2)

legend("topleft", 
       legend = c("CASM", "MX", "USWC"),
       col = c("red", "blue", "black"),         # Line colors
       lty = 1,                                 # Solid line type
       lwd = 2,                                 # Line width
       fill = c(adjustcolor("pink", 0.3),       # Shaded box colors
                adjustcolor("lightblue", 0.3), 
                adjustcolor("gray", 0.3)),
       border = NA,                             # Remove borders around shaded boxes
       bty = "n")                               # Optional: Remove outer legend box border

#dev.off()


