# run SCR script Part 1 for Google Cloud Workstation
# Author: K. Alexandra Curtis <alex.curtis@noaa.gov>

# # (1) run terminal code required to use nimble
# system2("sudo", args = c("apt", "update"))
# system2("sudo", args = c("apt", "install", "libglpk40", "-y"))

##load libraries
library(nimble)
library(nimbleEcology)
library(coda)
library(dplyr)

# import data
load("portadata.gc.rdata")

# run model
## define model
source("SCR1D.partialeff.log.sigreg.nimble.r")
## model boundary at 19.2 N (near Tenacatita,  Colima/Jalisco boundary), partial effort covariate
dl <- list(Y = y3d.casm.2124, traplocs = traplocs.2124, 
           ylim = c(model.limits$cumdist[1], model.limits$cumdist[6]), 
           K=4, effort = eff$dist.km)
out.nim.peff.2124.QAB.19 <- SCR1D.partialeff.fitfun(data = dl, M = 2000, ni = 50000, nb = 2000)

sink("summary nim.SCR1Ddist.partialeff.log.sigreg.casm.2124.19.txt")
cat("gelman.diag(out.nim.peff.2124.QAB.19$fit1$samples, multivariate = F)\n")
print(gelman.diag(out.nim.peff.2124.QAB.19$fit1$samples, multivariate = F))
cat("effectiveSize(out.nim.peff.2124.QAB.19$fit1$samples)\n")
print(effectiveSize(out.nim.peff.2124.QAB.19$fit1$samples))
cat("summary(out.nim.peff.2124.QAB.19$fit1$samples)\n")
print(summary(out.nim.peff.2124.QAB.19$fit1$samples))
sink()

pdf(file = "plot post nim.SCR1Ddist.partialeff.log.sigreg.casm.2124.19.pdf", width = 8.5, height = 11)
plot(out.nim.peff.2124.QAB.19$fit1$samples)
dev.off()

# save results
# Extract samples as a standard matrix or mcmc.list
samples.nim.pleff.sigreg.2124.QAB.19 <- out.nim.peff.2124.QAB.19$fit1
save(samples.nim.pleff.sigreg.2124.QAB.19, file="nim.SCR1Ddist.partialeff.log.sigreg.casm.2124.QAB.19.rdata")
