# Spatial randomization test - results not included in report of 2021-2024 since 
# consistent with previously published, and only used to support SCR model approach 
# Author: K. Alexandra Curtis <alex.curtis@noaa.gov>

# load libraries
library(ggplot2)

# Compare all distances among within-individual inter-occasion-trap sightings to all distances among occasion-trap-level sightings over all individuals
## Randomize lat/lon/trapid within occasion after sampling single point per individual-occasion-trap, prior to selection of multi-year individuals
### remove missing locations, years
any(is.na(sight.casm$occ)) # FALSE
any(is.na(sight.casm$lat)) # FALSE
sightloc <- sight.casm %>% 
  # filter(!duplicated(paste(id, trapid, occ))) %>%   # what if keep only first sighting per occ-trapid of each individual? results are the same. 
  select(id, occ, lat, lon, trapid) %>%    # , traplat, traplon (add these and code below if missing lat/lon for some sightings)
  group_by(id) %>% mutate(ny=n_distinct(occ)) %>% 
  ungroup() 
  # mutate(lat = if_else(is.na(lat), traplat, lat),
  #        lon = if_else(is.na(lon), traplon, lon))

# define function for distance metric
dist.ind.mean <- function(x) {
  # filter multiyear individuals
  my.iyl <- x %>% filter(ny>1)                                    
  # calculate mean distance among within-individual sightings 
  d.ind <- my.iyl %>% 
    group_by(id) %>% tidyr::nest() %>% 
    mutate(dists = purrr::map(data, ~ geosphere::distm(data.frame(lon=.x$lon,lat=.x$lat))/1000),
           xdist = purrr::map(dists, function(x) mean(x[upper.tri(x)])))   # mean per individual
  return(mean(unlist(d.ind$xdist)))   
}

# permute and calculate metric
nrep <- 1000
truth <- as.numeric(rep(NA, nrep))
test <- as.numeric(rep(NA, nrep))
for (i in 1:nrep) {
  all.iyl <- sightloc %>% 
    group_by(id, occ) %>% sample_n(1) %>% ungroup() %>% arrange(occ, id) 
  truth[i] <- dist.ind.mean(all.iyl)
  perm <- all.iyl %>% 
    group_by(occ) %>% 
    mutate(r=sample(1:n()), id=id[r], ny=ny[r]) %>%
    ungroup() %>% select(-r)
  test[i] <- dist.ind.mean(perm)
}
rm(all.iyl, perm, i)

# plot results
df <- data.frame(d = c(truth, test), Type = factor(c(rep("true", nrep), rep("permuted", nrep))))
p1 <- 
  ggplot(df, aes(d)) +
  geom_histogram(aes(fill = Type), color = "black") +
  scale_fill_manual(values = c("true" = "gray50", "permuted" = "gray80")) +
  xlab("Distance (km)") + ylab("Frequency") +
  scale_x_continuous(limits = c(500, 1000), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0,1000), expand = c(0,0)) +
  theme_bw() + theme(panel.grid.major= element_line(colour="white"), panel.grid.minor = element_line(colour="white"),
                     plot.margin = margin(t = 5, r = 10, b = 5, l = 5, unit = "pt"),
                     text = element_text(size = 10),
                     legend.position = c(0.95, 0.95), legend.justification = c(1, 1)) +
  theme(legend.key.size = unit(0.4,"cm"))

ggsave(plot=p1, filename="Fig_PermutationTest_through2024_1000km.tiff", 
       width=100, height=80, units="mm", dpi = 1200, compression = "lzw")

sum(truth<test)/nrep  # 1.000
