# nimble 1D SCR model for closed populations, with partial effort covariate
# Author: K. Alexandra Curtis <alex.curtis@noaa.gov>

library(nimble)
library(nimbleEcology)

SCR1D.partialeff.code <- nimbleCode({
  for (j in 1:J) {
    for (k in 1:K) {
      lp0.pred[j,k] <- alpha0 + alpha1 * log(effcov[j,k]+1)
      logit(p0[j,k]) ~ dnorm(lp0.pred[j,k], sd = sigmareg)
      effcov[j,k] ~ dexp(1/effcov.x)   # effort prior mirrors data distribution
    }
  }
  alpha0 ~ dnorm(0, sd = 10)
  alpha1 ~ dnorm(0, sd = 5)
  sigmareg ~ dunif(0, 10)
  sigma ~ dunif(0,50)  # scale parameter of encounter function
  psi ~ dunif(0, 1)     # DA parameter: E(N) = M*psi
  
  for(i in 1:M) {
    z[i] ~ dbern(psi)  # Is individual real?
    s[i] ~ dunif(ylim[1], ylim[2]) # 1D spatial coordinate
    # # dist between activity center and trap
    # d[i,1:J] <- abs(s[i] - X[1:J])
    for(j in 1:J) {
      p[i,j,1:K] <- p0[j,1:K] * exp(-(s[i] - X[j])^2/(2*sigma^2))  # capture prob at trap j
      y[i,j,1:K] ~ dOcc_v(z[i], p[i,j,1:K], len = K)   # recode so vector is 1:M looped over J and K?
    }
  }
  
  N <- sum(z[1:M])   # realized abundance
  EN <- psi*M
})


SCR1D.partialeff.fitfun <- function(data, M, rand = NULL, ni = 2000, nb = 1000, cm.scr01d = NULL, cmcmc.scr01d = NULL, 
                                         model="SCR1D.partialeff.code") {
  
  if (!is.null(rand)) set.seed(rand)
  
  X <- data$traplocs
  K <- data$K
  J <- length(X)
  ylim <- data$ylim
  effcov <- array(data$effort, dim = c(J, K))
  effcov.x <- mean(data$effort, na.rm=TRUE)
  
  yobs <- data$Y
  nind <- dim(yobs)[1]
  y <- array(0, c(M,J,K))
  y[1:nind,,] <- yobs
  z <- c(rep(1, nind), rep(0, M - nind))
  sst <- runif(M, ylim[1], ylim[2])
  X2 <- array(X, dim = c(J, K))
  for (i in 1:nind) {
    if (sum(y[i, , ]) == 0) 
      next
    # allow for individuals only observed south of state-space
    Xobs <- mean(X2[y[i, , ] > 0])
    if (Xobs < ylim[1]) {
      sst[i] <- ylim[1]
    } else {
      sst[i] <- Xobs
    }
  }
  
  constants <- list(X = X, K = K, M = M, J = J, ylim = ylim)
  
  inits.scr01d <- function() {
    list(psi = runif(1, 0.5, 1),   # this seems to be ignored - maybe determined by z init?
         logit_p0 = array(runif(K*J, -5, 5), dim=c(J, K)),
         alpha0 = runif(1, -10, -3),
         alpha1 = runif(1, 0.2, 1.5),
         sigmareg = runif(1, 0.1, 1.5),
         sigma = runif(1, 3, 9),
         s = sst,
         z = z)
  }
  
  if (is.null(cm.scr01d)) {
    m.scr01d <- nimbleModel(get(model), constants = constants, data = list(y = y, effcov = effcov, effcov.x = effcov.x))   #, calculate = FALSE)
    conf.scr01d <- configureMCMC(m.scr01d, monitors = c("p0", "sigma", "N", "alpha0", "alpha1", "sigmareg" ,"effcov"), thin=1,
                                 monitors2=c("s","z","p0","sigma","psi","alpha0","alpha1","sigmareg"), thin2=25)
    conf.scr01d$removeSamplers(c("alpha0","alpha1")) 
    conf.scr01d$addSampler(target = c("alpha0","alpha1"), type = "AF_slice")
    mcmc.scr01d <- buildMCMC(conf.scr01d)
    cm.scr01d <- compileNimble(m.scr01d)
    cmcmc.scr01d <- compileNimble(mcmc.scr01d, project=m.scr01d)
  } else {
    cm.scr01d$y <- y
  }
  
  fit1 = runMCMC(cmcmc.scr01d, inits = inits.scr01d, niter = ni, nburnin = nb, nchains = 3, samplesAsCodaMCMC = TRUE)
  
  return(list(fit1 = fit1, cm.scr01d = cm.scr01d, cmcmc.scr01d = cmcmc.scr01d))
}
