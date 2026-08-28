# Simulate closed-population capture histories, with optional sex heterogeneity
# Author: K. Alexandra Curtis <alex.curtis@noaa.gov>

# hetfac: ratio of capture probabilities (e.g., F:M)
# hetprop: proportion of population constituted by numerator of hetfac (e.g., F)

sim.closed <- function (N = 100, K = 20, hetfac = 1, hetprop = 0.5, p = 0.2,
          array2d = FALSE, rnd = NULL) 
{
  if (!is.null(rnd)) 
    set.seed(rnd)
  if (length(p)==1)
    p <- rep(p, K)
  else if (length(p)!=K)
    stop("p must be length of 1 or K")
  P = matrix(c(rep(p, round(N*(1-hetprop))), rep(p*hetfac, round(N*hetprop))), nrow = N, ncol = K, byrow = TRUE)
  hetcov = c(rep(1, round(N*(1-hetprop))), rep(2, round(N*hetprop)))
  
  Y0 <- matrix(NA, nrow = N, ncol = K)
  for (i in 1:N) {
    Y0[i, 1:K] <- rbinom(K, 1, P[i,1:K])
  }
  ncaps <- rowSums(Y0)
  Y <- Y0[ncaps > 0, ]

  if(!array2d) {
    Y = rowSums(Y)
    Y0 = rowSums(Y0)
  }
    
  list(Y = Y, Y0 = Y0, N = N, p = p, hetfac = hetfac, hetprop = hetprop, K = K, 
       hetcov0 = hetcov, hetcov = hetcov[ncaps>0])
}

