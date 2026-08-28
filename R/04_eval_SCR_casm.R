# Model evaluation
# Author: K. Alexandra Curtis <alex.curtis@noaa.gov>

###*** RUN TIME FOR THIS SCRIPT ON A DESKTOP SYSTEM MIGHT TAKE UPWARDS OF 16 HOURS, DEPENDING ON HARDWARE.
###*** WE RAN THE POSTERIOR PREDICTIVE CHECK PORTION ON A GOOGLE CLOUD WORKSTATION INSTEAD USING THIS SCRIPT:
###     - 'Google  Cloud Workstation code/ppc.SCR.partialeff.log.sigreg.2124.19.gc.r'

###### GOODNESS OF FIT (POSTERIOR PREDICTIVE CHECK) ######

# load libraries
library(nimble)
library(nimbleEcology)
library(coda)

# define model
source("SCR1D.partialeff.log.sigreg.nimble.r")

# load results and data
attach("nim.SCR1Ddist.partialeff.log.sigreg.casm.2124.QAB.19.rdata")
attach("portadata.gc.rdata")

# set up and run GOF simulations
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
    s <- samples[ni, 40:2039]
    z <- samples[ni, 2042:4041]
    values(cm.scr1d, vars) <- c(samples[ni, "sigma"], s, z, boot::logit(p.perm))
    cm.scr1d$simulate(simNodes, includeData = TRUE)
    ysimflat <- values(cm.scr1d, simNodes[18038:36037])  # values(cm.scr1d, "y") returns unknown (sorted?) order of y nodes, so use known order
    ysim <- aperm(array(ysimflat, dim=c(dl$K, M, J)), perm = c(2,3,1))
    pflat <- values(cm.scr1d, simNodes[38:18037])
    p <- aperm(array(pflat, dim=c(dl$K, M, J)), perm = c(2,3,1))
    
    # Spatial randomness of activity centers
    sran <- runif(M, dl$ylim[1], dl$ylim[2])
    sbin <- rep(NA, nsbin)
    snewbin <- rep(NA, nsbin)
    for (b in 1:nsbin) {
      sbin[b] <- sum((z > 0) & (s > binlims[b]) & (s <= binlims[b+1]))
      snewbin[b] <- sum((z > 0) & (sran > binlims[b]) & (sran <= binlims[b+1]))
    }
    esbin <- sum(z)/nsbin
    ppc$Ts[ni] <- sum(pow(sqrt(sbin[1:nsbin]) - sqrt(esbin), 2))
    ppc$Tsnew[ni] <- sum(pow(sqrt(snewbin[1:nsbin]) - sqrt(esbin), 2))
    
    # Goodness of fit components
    ## captures
    ### captures per individual and trap across occasions
    esumyij <- matrix(NA, nrow=M, ncol=J)
    sumyij <- matrix(NA, nrow=M, ncol=J)
    sumynewij <- matrix(NA, nrow=M, ncol=J)
    for (j in 1:J) {
      for (i in 1:M) {
        esumyij[i,j] <- sum(z[i]*p[i,j,])  # expected, given model parameter values
        sumyij[i,j] <- sum(y[i,j,])   # observed
        sumynewij[i,j] <- sum(ysim[i,j,])   # simulated from model parameter values
      }
    }
    ### captures per trap across occasions and individuals
    esumyj <- colSums(esumyij)
    sumyj <- colSums(sumyij)
    sumynewj <- colSums(sumynewij)
    ### captures per individual across occasions and traps
    esumyi <- rowSums(esumyij)
    sumyi <- rowSums(sumyij)
    sumynewi <- rowSums(sumynewij)
    ### summary statistics for comparison
    ppc$Tij[ni] <- sum(pow(sqrt(sumyij) - sqrt(esumyij), 2))
    ppc$Tijnew[ni] <- sum(pow(sqrt(sumynewij) - sqrt(esumyij), 2))
    ppc$Tj[ni] <- sum(pow(sqrt(sumyj) - sqrt(esumyj), 2))
    ppc$Tjnew[ni] <- sum(pow(sqrt(sumynewj) - sqrt(esumyj), 2))
    ppc$Ti[ni] <- sum(pow(sqrt(sumyi) - sqrt(esumyi), 2))
    ppc$Tinew[ni] <- sum(pow(sqrt(sumynewi) - sqrt(esumyi), 2))
  }
})

# evaluate statistics
sum(ppc$Ts < ppc$Tsnew)/niter # 0.06
sum(ppc$Tj < ppc$Tjnew)/niter # 0.38
sum(ppc$Tij < ppc$Tijnew)/niter # 0.48
sum(ppc$Ti < ppc$Tinew)/niter  # 0.31

# save
save(ppc, file = "ppc.SCR1Ddist.partialeff.log.sigreg.casm.2124.QAB.19.rdata")


###### VARIANCE EXPLAINED BY EFFORT ######

# load libraries
library(coda)

# load results
attach("nim.SCR1Ddist.partialeff.log.sigreg.casm.2124.QAB.19.rdata")

# extract posteriors
alpha0 <- as.matrix(samples.nim.pleff.sigreg.2124.QAB.19$samples)[,"alpha0"]
alpha1 <- as.matrix(samples.nim.pleff.sigreg.2124.QAB.19$samples)[,"alpha1"]
sigmareg <- as.matrix(samples.nim.pleff.sigreg.2124.QAB.19$samples)[,"sigmareg"]
effcov <- as.matrix(samples.nim.pleff.sigreg.2124.QAB.19$samples)[,4:39]

# calculate predicted p0 component 
lp0.pred <- alpha0 + alpha1 * log(effcov+1)

# calculate component variances
var.lp0.pred <- apply(lp0.pred, MARGIN = 1, FUN = var)
var.sr <- sigmareg^2

# percentage variance explained by effort
pvar.ex <- var.lp0.pred/(var.lp0.pred + var.sr)
mean(pvar.ex)  # 0.89
