# Post-process CASM SCR results, calculate Mexico
# Author: K. Alexandra Curtis <alex.curtis@noaa.gov>

# load libraries
library(nimble)
library(nimbleEcology)
library(coda)


###### POST-PROCESS CASM SCR ######

# post-process CASM abundance with boundary uncertainty
attach("nim.SCR1Ddist.partialeff.log.sigreg.casm.2124.QAB.19.rdata")
attach("nim.SCR1Ddist.partialeff.log.sigreg.casm.2124.QAB.21.rdata")
## check parameter indices
# colnames(as.matrix(samples.nim.pleff.sigreg.2124.QAB.19$samples2)[,c(40,2039,2042,4041)])
## model domain boundary at 19.2 N
z <- as.matrix(samples.nim.pleff.sigreg.2124.QAB.19$samples2)[,2042:4041]
s <- as.matrix(samples.nim.pleff.sigreg.2124.QAB.19$samples2)[,40:2039] 
nlive.peff.casm.19 <- matrix(NA,dim(z)[1],20)
nbounds <- seq(25,28.2, length.out=20) 
for (i in 1:dim(z)[1]) {
  slive <- s[i,][as.logical(z[i,])]
  for (j in 1:length(nbounds)) {
    nlive.peff.casm.19[i,j] <- sum(slive <= nbounds[j])   # southern end of state-space limit coincides with PL
  }
}
## model domain boundary at BAB (lower envelope simulation)
z <- as.matrix(samples.nim.pleff.sigreg.2124.QAB.21$samples2)[,2042:4041]
s <- as.matrix(samples.nim.pleff.sigreg.2124.QAB.21$samples2)[,40:2039]
nlive.peff.casm.21 <- matrix(NA,dim(z)[1],20)
nbounds <- seq(25,28.2, length.out=20)
for (i in 1:dim(z)[1]) {
  slive <- s[i,][as.logical(z[i,])]
  for (j in 1:length(nbounds)) {
    nlive.peff.casm.21[i,j] <- sum(slive <= nbounds[j])  
  }
}
## correction factor
#  CF = 1/(1+ %bias), with birthdeaths + nocalfs + breeding-area-only census rnorm(3.0,5.0), and sexhet rnorm(-10.9,6.3)
cf <- 1/(1+rnorm(5760,.030,.050)-rnorm(5760,.109,.063))
## summaries
### 19
mean(nlive.peff.casm.19*cf) # 1453, SD=156 (CV 10.8%)
quantile(nlive.peff.casm.19*cf, 0.2)  # 1320
### 21
mean(nlive.peff.casm.21*cf)   # 1340 (CV 10.7%)
quantile(nlive.peff.casm.21*cf, 0.2)  # 1219
## save
save(nlive.peff.casm.19, nlive.peff.casm.21, cf, file="Nprocessed.SCR1Dpeff.2124.rdata")


###### MAINLAND MEXICO ABUNDANCE ######

# Mainland Mexico by subtraction
# most recent estimate from C & B 2025 for 2021-2024: x=7021, se=181.5
attach("Nprocessed.SCR1Dpeff.2124.rdata")

## simulation-based
# lx <- log((7021^2)/sqrt(7021^2+181.5^2))
# ls <- sqrt(log(1+181.5^2/7021^2))
# n.wc <- rlnorm(length(nlive.casm.19), meanlog=lx, sdlog=ls)
# n.mmex <- n.wc - as.vector(cf * nlive.casm.19)
# mean(n.mmex)   # 5600
# sd(n.mmex)/mean(n.mmex)   # 4.3%

## analytical (used in Curtis et al 2022)
### base model
n.mmex <- 7021 - mean(nlive.peff.casm.19*cf)   # 5568
sd.mmex <- sqrt(181.5^2 + sd(nlive.peff.casm.19*cf)^2)   # 240; cv = 4.3%
lx <- log((n.mmex^2)/sqrt(n.mmex^2+sd.mmex^2))
ls <- sqrt(log(1+sd.mmex^2/n.mmex^2))
round(qlnorm(0.2, lx, ls))   # 5365
### upper envelope
n.mmex <- 7021 - mean(nlive.peff.casm.21*cf)   # 5681
sd.mmex <- sqrt(181.5^2 + sd(nlive.peff.casm.21*cf)^2)   # 231; cv = 4.1%
lx <- log((n.mmex^2)/sqrt(n.mmex^2+sd.mmex^2))
ls <- sqrt(log(1+sd.mmex^2/n.mmex^2))
round(qlnorm(0.2, lx, ls))   # 5486

detach(2)
