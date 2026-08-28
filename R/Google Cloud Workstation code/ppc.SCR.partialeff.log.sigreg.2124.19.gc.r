# Posterior predictive check script for Google Cloud Workstation
# Author: K. Alexandra Curtis <alex.curtis@noaa.gov>

# # (1) run terminal code required to use nimble
# system2("sudo", args = c("apt", "update"))
# system2("sudo", args = c("apt", "install", "libglpk40", "-y"))

# load libraries
library(nimble)
library(nimbleEcology)
library(coda)

# define model
source("SCR1D.partialeff.log.sigreg.nimble.r")

# load data
attach("nim.SCR1Ddist.partialeff.log.sigreg.casm.2124.QAB.19.rdata")
attach("portadata.gc.rdata")

samples <- as.matrix(samples.nim.pleff.sigreg.2124.QAB.19$samples2)
niter <- dim(samples)[1]
## Ensure we have the nodes needed to simulate new datasets
M <- 2000
dl <- list(Y = y3d.casm.2124, traplocs = traplocs.2124, 
           ylim = c(model.limits$cumdist[1], model.limits$cumdist[6]), 
           K = 4, effort = eff$dist.km)
X <- dl$traplocs; J <- length(X)
nind <- dim(dl$Y)[1]
y <- array(0, c(M, J, dl$K))
y[1:nind,,] <- dl$Y
effcov <- array(dl$effort, dim = c(J, dl$K))
effcov.x <- mean(dl$effort, na.rm=TRUE)
constants <- list(X = X, K = dl$K, M = M, J = J, ylim = dl$ylim)
m.scr1d <- nimbleModel(SCR1D.partialeff.code, constants = constants, 
                       data = list(y = y, effcov = effcov, effcov.x = effcov.x))
cm.scr1d <- compileNimble(m.scr1d)
dataNodes <- m.scr1d$getNodeNames(dataOnly = TRUE, returnScalarComponents = TRUE)   # effcov, y(indexing increments for last=K, first=M, then middle=J)
parentNodes <- m.scr1d$getParents(dataNodes, stochOnly = TRUE)   # sigma, s, z, logit_p0
simNodes <- m.scr1d$getDependencies(parentNodes, self = FALSE)   # includes N, p0, p, y (latter two indexing in groups of K, then M, then J)
# var name order from mcmc in samples: alpha0, alpha1, p0, psi, s, sigma, sigmareg, z
vars <- parentNodes
nsbin <- 10
binlims <- dl$ylim[1] + (0:nsbin)*((dl$ylim[2]-dl$ylim[1])/nsbin)
ppc <- data.frame(Ts = rep(NA, niter), Tsnew = rep(NA, niter), 
                  Tij = rep(NA, niter), Tijnew = rep(NA, niter),
                  Ti = rep(NA, niter), Tinew = rep(NA, niter),
                  Tj = rep(NA, niter), Tjnew = rep(NA, niter))
set.seed(1)
system.time({
  for(ni in 1:niter) {
    # get values for parentNodes from samples
    ## logit_p0 indexing in parentNodes does not line up with indexing in samples, so have to rearrange latter first
    p.perm <- as.vector(aperm(array(samples[ni, 3:38], c(9,4)), perm=c(2,1)))
    values(cm.scr1d, vars) <- c(samples[ni, "sigma"], samples[ni, 40:2039], samples[ni, 2042:4041], boot::logit(p.perm))
    cm.scr1d$simulate(simNodes, includeData = TRUE)
    ysimflat <- values(cm.scr1d, simNodes[18038:36037])  # values(cm.scr1d, "y") returns unknown order of y nodes, so use known order
    ysim <- aperm(array(ysimflat, dim=c(dl$K, M, J)), perm = c(2,3,1))
    pflat <- values(cm.scr1d, simNodes[38:18037])
    p <- aperm(array(pflat, dim=c(dl$K, M, J)), perm = c(2,3,1))
    s <- values(cm.scr1d, "s")
    z <- values(cm.scr1d, "z")
    sran <- runif(M, dl$ylim[1], dl$ylim[2])
    # Goodness of fit components
    ## spatial randomness of activity centers
    sbin <- rep(NA, nsbin)
    snewbin <- rep(NA, nsbin)
    for (b in 1:nsbin) {
      sbin[b] <- sum((z > 0) & (s > binlims[b]) & (s <= binlims[b+1]))
      snewbin[b] <- sum((z > 0) & (sran > binlims[b]) & (sran <= binlims[b+1]))
    }
    esbin <- sum(z)/nsbin
    ppc$Ts[ni] <- sum(pow(sqrt(sbin[1:nsbin]) - sqrt(esbin), 2))
    ppc$Tsnew[ni] <- sum(pow(sqrt(snewbin[1:nsbin]) - sqrt(esbin), 2))
    ## observations
    esumyij <- matrix(NA, nrow=M, ncol=J)
    sumyij <- matrix(NA, nrow=M, ncol=J)
    sumynewij <- matrix(NA, nrow=M, ncol=J)
    for (j in 1:J) {
      for (i in 1:M) {
        esumyij[i,j] <- sum(z[i]*p[i,j,])
        sumyij[i,j] <- sum(y[i,j,])
        sumynewij[i,j] <- sum(ysim[i,j,])
      }
    }
    esumyj <- colSums(esumyij)
    sumyj <- colSums(sumyij)
    sumynewj <- colSums(sumynewij)
    esumyi <- rowSums(esumyij)
    sumyi <- rowSums(sumyij)
    sumynewi <- rowSums(sumynewij)
    ppc$Tij[ni] <- sum(pow(sqrt(sumyij) - sqrt(esumyij), 2))
    ppc$Tijnew[ni] <- sum(pow(sqrt(sumynewij) - sqrt(esumyij), 2))
    ppc$Tj[ni] <- sum(pow(sqrt(sumyj) - sqrt(esumyj), 2))
    ppc$Tjnew[ni] <- sum(pow(sqrt(sumynewj) - sqrt(esumyj), 2))
    ppc$Ti[ni] <- sum(pow(sqrt(sumyi) - sqrt(esumyi), 2))
    ppc$Tinew[ni] <- sum(pow(sqrt(sumynewi) - sqrt(esumyi), 2))
  }
})
sum(ppc$Ts < ppc$Tsnew)/niter
sum(ppc$Tj < ppc$Tjnew)/niter
sum(ppc$Tij < ppc$Tijnew)/niter
sum(ppc$Ti < ppc$Tinew)/niter

save(ppc, file = "ppc.SCR1Ddist.partialeff.log.sigreg.casm.2124.QAB.19.rdata")
