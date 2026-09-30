# Trend estimation for CASM
# Author: K. Alexandra Curtis <alex.curtis@noaa.gov>
# - run ChaoMth on the two time periods, restricted to Central America sightings
# - sensitivity to subsampling second period to equivalent capture probability

# load libraries
library(CMRutils)

# data prep
source("01_prepdata_casm.r")

# create/edit capture history objects
nsam <- 2000
## 2004-06
y2d.cam.0406 <- array2binom(apply(y3d.cam.0406,MARGIN=c(1,3),sum))
## 2022-24
y3d.cam.2224 <- dim_drop0(dim_drop0(y3d.cam[,,dimnames(y3d.cam)$occ %in% 2022:2024], 1), 2)
y2d.cam.2224 <- array2binom(apply(y3d.cam.2224,MARGIN=c(1,3),sum))

# ChaoMth estimate for 2004-06
mrc.0406 <- run_CAPTURE(y2d.cam.0406, keep.out=FALSE)
mth.0406 <- mrc.0406$N.ChaoMth
mth.se.0406 <- mrc.0406$seN.ChaoMth
## distribution for n
lx <- log((mth.0406^2)/sqrt(mth.0406^2 + mth.se.0406^2))
ls <- sqrt(log(1 + mth.se.0406^2/mth.0406^2))
n.0406 <- rlnorm(nsam, meanlog=lx, sdlog=ls)

# ChaoMth estimate for 2022-24
mrc.2224 <- run_CAPTURE(y2d.cam.2224, keep.out=FALSE)
mth.2224 <- mrc.2224$N.ChaoMth
mth.se.2224 <- mrc.2224$seN.ChaoMth
## distribution for n
lx <- log((mth.2224^2)/sqrt(mth.2224^2 + mth.se.2224^2))
ls <- sqrt(log(1 + mth.se.2224^2/mth.2224^2))
n.2224 <- rlnorm(nsam, meanlog=lx, sdlog=ls)

# raw population growth rate
lambda <- (n.2224/n.0406)^(1/18)
100*mean(lambda-1)
100*sd(lambda)
# approximate results (~0% with ~3%sd) robust to other combinations of years in 2021:2024, including all four

# resample 2022-24 to equivalent capture probability and refit
# distributions of estimated mean p
pcum.0406 <- nrow(y2d.cam.0406) / n.0406
px.0406 <- 1 - (1-pcum.0406)^(1/3)
pcum.2224 <- nrow(y2d.cam.2224) / n.2224
px.2224 <- 1 - (1-pcum.2224)^(1/3)
## ratio of mean p between two sampling periods by which to subsample captures
p.sub <- px.0406/px.2224
## MC
df <- data.frame(n.sub=rep(NA, nsam), pcum.sub=rep(NA,nsam))
for (i in 1:nsam) {
  # subsample captures
  y2d.cam.2224.sam <- matrix(0, nrow=nrow(y2d.cam.2224), ncol=ncol(y2d.cam.2224))
  caps <- apply(y2d.cam.2224, MARGIN=2, function(x) which(x>0))
  # subsample captures where p.sub < 1
  for (oi in 1:ncol(y2d.cam.2224)) {
    if (p.sub[i]<1) {
      # # enforce a minimum sample size to ensure successful CAPTURE execution
      # isub <- sample(caps[[oi]], size =  max(30, round(p.sub[i] * length(caps[[oi]]))))
      isub <- sample(caps[[oi]], size = round(p.sub[i] * length(caps[[oi]])))
    } else { isub <- caps[[oi]] }
    y2d.cam.2224.sam[isub, oi] <- 1
  }
  # drop all-zero capture histories
  ri <- rowSums(y2d.cam.2224.sam)>0
  y2d.cam.2224.sam <- y2d.cam.2224.sam[ri,]
  # estimate
  ## catch errors resulting from run_CAPTURE
  mrc <- tryCatch({
    run_CAPTURE(y2d.cam.2224.sam, keep.out = TRUE)
    }, 
    error = function(e) {
      message("Warning: run_CAPTURE failed. Returning NULL. Error: ", e$message)
      setwd("C:/Users/alex.curtis/Data/Github/Mn-sandbox")
      return(NA)
    }
  )
  if(is.list(mrc)) {
    df$n.sub[i] <- mrc$N.ChaoMth
    ## estimated cumulative p
    df$pcum.sub[i] <- nrow(y2d.cam.2224.sam)/df$n.sub[i]
  } else {
    df$n.sub[i] <- NA
    df$pcum.sub[i] <- NA
  }
}

sum(!(is.na(df$n.sub))) # 1915
# compare df$psam with p.0406
pratio.new <- pcum.0406/df$pcum.sub
mean(pratio.new, na.rm=TRUE)
# if acceptable, recalc pop growth rate
lambda.sub <- (df$n.sub/n.0406)^(1/18)
(mean(lambda.sub-1, na.rm=TRUE))*100
sd(lambda.sub, na.rm=TRUE)
# 1915 successful executions, with new pratio of 1.09, similar mean (-0.41%) and much higher SE (3.7%)
save(df, lambda, lambda.sub, n.0406, n.2224, pcum.0406, pcum.2224, y2d.cam.0406, y2d.cam.2224, file="trend.2006.2024.rdata")
