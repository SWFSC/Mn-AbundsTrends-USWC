# Assess bias due to sex heterogeneity in sighting prob off CAm/SMex
# (updated 2/10/26 to use both genetic and field sex)
# Author: K. Alexandra Curtis <alex.curtis@noaa.gov>

# run 01_prepdata_casm.r through Line 93

# define function for simulating closed-population mark-recapture sampling and analysis with sex heterogeneity
source("sim_CRclosed.r")


############ QUANTIFY SEX HETEROGENEITY IN CAPTURE PROBABILITY ###############

# Import CCES 2024 sex data from MMGP
sex.mmgp <- readxl::read_excel("MMGP CCES 2024 sexes.xlsx") %>%
  filter(SPNAME=="M nov" & CONTACT=="Jeff Moore" & !is.na(`Genetic Sex`) & `Genetic Sex` %in% c("female", "male")) %>% 
  select(Sample.Number=FieldID, sex.mmgp = `Genetic Sex`, Year=Year, Month=Month, Day=Day) %>%
  mutate(sex.mmgp = if_else(sex.mmgp=="female", "F", "M"))

# Genetic sex (from CRC biop table)
## edit and subset biop
biop.ed <- biop %>% 
  # omit or edit individuals with >1 different haplotype result
  filter(Sample.Number!="CRC-20170809-ZIP-03") %>% 
  filter(!(id %in% c(16757, 17302))) %>% 
  # omit individuals with >1 genetic sex result and no haplo info to confirm identity
  filter(!(id %in% c(5043, 10026))) %>% 
  # subset relevant samples
  filter(!is.na(id)) %>% 
   filter(Certainty.of.Sample.ID %in% c("Positive","Probable") | is.na(Certainty.of.Sample.ID)) %>%
  left_join(sex.mmgp) %>%
  mutate(SWFSC.Sex = if_else(!is.na(sex.mmgp), sex.mmgp, SWFSC.Sex)) %>%
  filter(SWFSC.Sex %in% c("F","M") | OSU.Sex %in% c("F","M"))

# Behavioral ("field") sex (from sight table)
# Vast majority of behavioral sex data are female.
sight.sex <- sight %>% 
  mutate(Beh.Sex = if_else(Beh.Sex=="UK", NA_character_, Beh.Sex)) %>% 
  filter((!is.na(Beh.Sex)) | (!is.na(AGEREP) & grepl("M", AGEREP, fixed=TRUE))) %>% 
  mutate(Beh.Sex = substr(Beh.Sex, 1, 1), 
         Beh.Sex.AR = if_else(grepl("M", AGEREP, fixed=TRUE), "F", NA_character_)) %>%
  # edit CRC-confirmed unreliable entry for Beh.Sex
  mutate(Beh.Sex = if_else(id==15230, NA_character_, Beh.Sex))

# Filter and combine genetic and behavioral sex
# Get all indiv sexed either by gen or beh on West Coast or before study period (i.e., independent of study sightings)
# - all biop samples from WC and from pre-study-period CASM
# - beh sex from WC and from pre-study-period CASM 
#   (use only M+ for AGEREP to eliminate records borrowing info from other sightings that year)
# - eliminate individuals where gen and beh sex disagree
## Filter genetic sex data
sex.biop <- biop.ed %>% 
  # spatiotemporal filter
  filter(is.na(lat) | (lat>30 & lat<49) | Year<2021) %>%    # NA latitudes currently all in GEMM lab data, assume off OR
  # retain best sex info from biopsy table
  # (1) SWFSC.Sex for SWFSC.GENDETER 11, 2; (2) All other sources, dropping any further individuals with conflicts
  group_by(id) %>% 
  summarize(sex.g = case_when(any(!is.na(SWFSC.Sex) & SWFSC.GENDETER %in% c(2,11)) ~ SWFSC.Sex[!is.na(SWFSC.Sex) & SWFSC.GENDETER %in% c(2,11)][1],
                              sum(c("M","F") %in% c(SWFSC.Sex, OSU.Sex)) > 1 ~ NA_character_,
                              TRUE ~ setdiff(unique(na.omit(c(SWFSC.Sex, OSU.Sex))),"U")[1])) %>% 
  filter(!is.na(sex.g))
## Filter behavioral sex data
sex.beh <- sight.sex %>% 
  filter(!grepl("?", AGEREP, fixed=TRUE)) %>%
  # keep only M+ for AGEREP info
  mutate(Beh.Sex.AR = if_else(grepl("M+", AGEREP, fixed=TRUE), "F", NA_character_)) %>% 
  # spatiotemporal filter
  filter(is.na(lat) | (lat>30 & lat<49) | year<2021) %>%    # Checked that no Nov/Dec 2020 CASM sightings sex data included
  filter(rgn %in% c("Central America", "S Mexico", "US West Coast")) %>% 
  # return unique non-NA values per id for combined Beh.Sex and Beh.Sex.AR 
  group_by(id) %>% 
  reframe(sex.b = unique(na.omit(c(Beh.Sex, Beh.Sex.AR)))) %>% 
  # make sure no conflicting behavioral sex assignments remain
  group_by(id) %>% filter(n()==1) %>% ungroup()
## Combine genetic and behavioral sex info
sex.info <- full_join(sex.biop, sex.beh) %>% 
  # drop any conflicts between genetic and field sex
  filter(is.na(sex.g) | is.na(sex.b) | sex.g==sex.b) %>% 
  # combine in one field and add source of sex info
  rowwise() %>% 
  mutate(sex = na.omit(c(sex.g, sex.b))[1],
         source = case_when(!is.na(sex.g) & !is.na(sex.b) ~ "both",
                            !is.na(sex.g) ~ "biop",
                            !is.na(sex.b) ~ "field")) %>% 
  ungroup() %>% select(-(c(sex.g, sex.b)))

# estimate CASM heterogeneity from capture histories
## merge sex data with sightings
sight.casm.sex <- left_join(sight.casm, sex.info) %>% 
  filter(!is.na(sex))
### 2021-2024
sex.het <- sight.casm.sex %>% 
  filter(occ %in% 2021:2024) %>%
  group_by(id,sex) %>% summarize(nocc=n_distinct(occ)) %>% 
  group_by(sex) %>% summarize(n=n(), ns=mean(nocc), ns.se=sd(nocc)/sqrt(n), nalt=sum(nocc>1))
# sex.het
# sex       n    ns   ns.se   nalt
# F        96  1.22   0.0517   17
# M        90  1.52   0.0794   35 
#### approximate distribution for simulations (approximate ratio as lognormal)
sex.het.med <- (sex.het$ns[2]-1)/(sex.het$ns[1]-1) # 2.39
cv.analytic <- sqrt((sex.het$ns.se[1]/(sex.het$ns[1]-1))^2+(sex.het$ns.se[2]/(sex.het$ns[2]-1))^2)  # 0.281
sigma <- sqrt(log(cv.analytic^2+1))
mu <- log(sex.het.med)
mean <- exp(mu + 0.5*sigma^2)  # 2.480
quants <- qlnorm(p=c(0.16, 0.84), meanlog=mu, sdlog=sigma)  # 1.81, 3.14
hetfac.2124 <- c(quants[1], sex.het.med, quants[2])
### clean up
rm(sex.het.med, cv.analytic, sigma, mu, mean, quants)

############ ESTIMATE EXPECTED BIAS DUE TO SEX HETEROGENEITY THROUGH SIMULATION ###############

# estimate potential bias in population estimate ignoring sex heterogeneity
# (necessarily an iterative process to get customized p given bias)
library(Rcapture)
## simulate 
nsim <- 5000
N <- 1500   # not sensitive to N in relevant param space (neg bias for N=1000 is only 0.2-0.3% smaller than for N=2000)
df.biassims <- data.frame(
  period = rep("2021-2024", 3*nsim),
  ny = rep(4, 3*nsim),   # number of occasions
  hetfac.level = rep(c("low", "med", "hi"), each=nsim),   
  hetfac = rep(hetfac.2124*0.91, each=nsim),   # relative sex heterogeneity in capture probability
  ip = rep(1, 3*nsim),   # index capture probability options
  het.sim = rep(NA, 3*nsim),
  p.sim = rep(NA, 3*nsim),
  p.sim.fem = rep(NA, 3*nsim),
  p.sim.mal = rep(NA, 3*nsim),
  mt.bias = rep(NA, 3*nsim),
  lp.bias = rep(NA, 3*nsim)
)
# set capture probabilities for simulations
# psim <- data.frame(p.2123 = c(0.27,0.22,0.19,NA), p.2124=c(0.24,0.20,0.17,0.15))
# bias-corrected (one iteration):
psim <- data.frame(p.2124 = c(.21, .17, .15, .13)*1.25)
# stability of results: differences among runs of 5000 sims only 0.1-0.2% 
# sensitivity to p:
# halving p results in ca. 1.3% (for K=4) or 2% (for K=3) less neg bias
# doubling p also results in 0.15-.5% (K=3) or 0.4-0.9% (K=4) less neg bias

for (ri in 1:nrow(df.biassims)) {
  
  # set number of occasions
  K = df.biassims$ny[ri]
  
  # simulate sampling
  test <- sim.closed(N = N, K = K, 
                     p = psim[,df.biassims$ip[ri]][1:K], 
                     hetfac = 1/(df.biassims$hetfac[ri]), array2d=TRUE)
  
  # Mt and L-P population estimates
  mr <- closedp(test$Y)
  n.mt = mr$results[2,1]
  n.lp = sum(test$Y[,2])*sum(test$Y[,3])/table(rowSums(test$Y[,2:3]))[3]
  
  # simulated observed sex het
  ## get sexes by replicating sex assignment from sim.closed function
  sex <- c(rep("A", round(test$N*(1-test$hetprop))), rep("B", round(test$N*test$hetprop)))
  ## capture histories by sex (including individuals with zero caps)
  ch.a0 <- test$Y0[sex=="A",]   
  ch.b0 <- test$Y0[sex=="B",]   
  ## capture histories by sex excluding zeros
  ch.a <- ch.a0[rowSums(ch.a0) > 0, ] 
  ch.b <- ch.b0[rowSums(ch.b0) > 0, ]
  ## mean relative recapture frequency of sex A over B
  het.sim <- (mean(rowSums(ch.a))-1)/(mean(rowSums(ch.b))-1)
  ## realized capture probability by sex
  p.sim <- nrow(test$Y)/N
  p.sim.fem <- nrow(ch.b)/(N*(1-test$hetprop))
  p.sim.mal <- nrow(ch.a)/(N*test$hetprop)
  
  # store results
  df.biassims$het.sim[ri] <- het.sim
  df.biassims$p.sim[ri] <- p.sim
  df.biassims$p.sim.fem[ri] <- p.sim.fem
  df.biassims$p.sim.mal[ri] <- p.sim.mal
  df.biassims$mt.bias[ri] <- n.mt/N
  df.biassims$lp.bias[ri] <- n.lp/N
  
  rm(test)
}

# Adjusted hetfac for sims to account for the fact that probability of recapture given capture is not unbiased estimator of relative prob cap, 
# because already has a sex het component built into the conditional prob
df.biassims %>% group_by(period, hetfac, hetfac.level) %>% 
  summarize(hetfac.sim.cap=median(p.sim.mal/p.sim.fem), hetfac.sim.recap=median(het.sim)) %>% 
  ungroup() %>% select(-hetfac)
# Do the summary statistics for the simulations (hetfac.sim.recap) closely approximate the observed and its uncertainty interval?
# there is almost certainly an algebraic solution but I did this numerically and approximately through manual iteration
# period    hetfac.level hetfac.sim.cap hetfac.sim.recap
# 2021-2024 low                    1.46             1.77
# 2021-2024 med                    1.84             2.38
# 2021-2024 hi                     2.34             3.22

df.biassims.summary <- df.biassims %>% 
  group_by(period, ny, hetfac, hetfac.level) %>% 
  summarize(mt.bias.x = mean(mt.bias), mt.bias.med = median(mt.bias),
            lp.bias.x = mean(lp.bias), lp.bias.med = median(lp.bias)) %>% 
  ungroup() %>% select(-hetfac)

# df.biassims.summary
# period       ny hetfac.level mt.bias.x mt.bias.med lp.bias.x lp.bias.med
# 2021-2024     4 low              0.951       0.950     0.955       0.941
# 2021-2024     4 med              0.892       0.891     0.896       0.881
# 2021-2024     4 hi               0.827       0.824     0.829       0.813
round((df.biassims.summary$mt.bias.med*100-100)[2], 1)   # -10.9
round(mean(abs((df.biassims.summary$mt.bias.med*100-100)[2] - 
                 (df.biassims.summary$mt.bias.med*100-100)[c(1,3)])), 1)  # [1] 6.3 (approximately symmetrical)
# Bias estimated for +/-SE of estimated heterogeneity is approximately symmetrical about the bias estimated for the median.
# Parameterize uncertainty in bias from sex heterogeneity as 
# -10.9%(+/-6.3) -> bias cf estimated as 1/(1 - calf + openpop + rnorm(-0.109,0.063))

# LP similar in performance.
