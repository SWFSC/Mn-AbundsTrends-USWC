# Run SCR Models for CASM stock
# Author: K. Alexandra Curtis <alex.curtis@noaa.gov>

###*** RUN TIME FOR THIS SINGLE SCRIPT ON A DESKTOP SYSTEM MIGHT TAKE UPWARDS OF 30 HOURS, DEPENDING ON HARDWARE.
###*** INDIVIDUAL MODELS WERE RUN SEPARATELY ON A GOOGLE CLOUD WORKSTATION USING THE FOLLOWING TWO SCRIPTS (CA 3 H RUNTIME EA) 
###     - 'Google  CLoud Workstation code/assess.SCR.partialeff.log.sigreg.2124.19.gc.r'
###     - 'Google  CLoud Workstation code/assess.SCR.partialeff.log.sigreg.2124.21.gc.r'

# fit spatial model with distance coordinate, partial effort covariate with random effect

# load libraries
library(nimble)
library(nimbleEcology)
library(coda)

# define model
source("SCR1D.partialeff.log.sigreg.nimble.r")

# run
## model boundary at 19.2 N (near Tenacatita,  Colima/Jalisco boundary)
dl <- list(Y = y3d.casm.2124, traplocs = traplocs.2124, 
           ylim = c(model.limits$cumdist[1], model.limits$cumdist[6]), 
           K = 4, effort = eff$dist.km)
out.nim.peff.2124.QAB.19 <- SCR1D.partialeff.fitfun(data = dl, M = 2000, ni = 50000, nb = 2000)
### summarize posteriors
sink("summary nim.SCR1Ddist.partialeff.log.sigreg.casm.2124.19.txt")
cat("gelman.diag(out.nim.peff.2124.QAB.19$fit1$samples, multivariate = F)\n")
print(gelman.diag(out.nim.peff.2124.QAB.19$fit1$samples, multivariate = F))
cat("effectiveSize(out.nim.peff.2124.QAB.19$fit1$samples)\n")
print(effectiveSize(out.nim.peff.2124.QAB.19$fit1$samples))
cat("summary(out.nim.peff.2124.QAB.19$fit1$samples)\n")
print(summary(out.nim.peff.2124.QAB.19$fit1$samples))
sink()
### plot posteriors
pdf(file = "plot post nim.SCR1Ddist.partialeff.log.sigreg.casm.2124.19.pdf", width = 8.5, height = 11)
plot(out.nim.peff.2124.QAB.19$fit1$samples)
dev.off()
### save results
# Extract samples as a standard matrix or mcmc.list
samples.nim.pleff.sigreg.2124.QAB.19 <- out.nim.peff.2124.QAB.19$fit1
save(samples.nim.pleff.sigreg.2124.QAB.19, file="nim.SCR1Ddist.partialeff.log.sigreg.casm.2124.QAB.19.rdata")

## model boundary at Bahia de Banderas
dl <- list(Y = y3d.casm.2124, traplocs = traplocs.2124, 
           ylim = c(model.limits$cumdist[1], model.limits$cumdist[8]), 
           K = 4, effort = eff$dist.km)
out.nim.peff.2124.QAB.21 <- SCR1D.partialeff.fitfun(data = dl, M = 2000, ni = 50000, nb = 2000)
### summarize posteriors
sink("summary nim.SCR1Ddist.partialeff.log.sigreg.casm.2124.21.txt")
cat("gelman.diag(out.nim.peff.2124.QAB.21$fit1$samples, multivariate = F)\n")
print(gelman.diag(out.nim.peff.2124.QAB.21$fit1$samples, multivariate = F))
cat("effectiveSize(out.nim.peff.2124.QAB.21$fit1$samples)\n")
print(effectiveSize(out.nim.peff.2124.QAB.21$fit1$samples))
cat("summary(out.nim.peff.2124.QAB.21$fit1$samples)\n")
print(summary(out.nim.peff.2124.QAB.21$fit1$samples))
sink()
### plot posteriors
pdf(file = "plot post nim.SCR1Ddist.partialeff.log.sigreg.casm.2124.21.pdf", width = 8.5, height = 11)
plot(out.nim.peff.2124.QAB.21$fit1$samples)
dev.off()
### save results
samples.nim.pleff.sigreg.2124.QAB.21 <- out.nim.peff.2124.QAB.21$fit1
save(samples.nim.pleff.sigreg.2124.QAB.21, file="nim.SCR1Ddist.partialeff.log.sigreg.casm.2124.QAB.21.rdata")
