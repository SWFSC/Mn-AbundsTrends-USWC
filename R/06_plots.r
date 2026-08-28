# summaries and plots for Tech Memo
# Author: K. Alexandra Curtis <alex.curtis@noaa.gov>

# create data sets
source("01_prepdata_casm.r")

# recaptures
## temporal
table(rowSums(array2binom(apply(y3d.casm.2124, sum, MARGIN=c(1,3))))) 
## spatial
table(rowSums(array2binom(apply(y3d.casm.2124, sum, MARGIN=c(1,2))))) 

# load general libraries
require(ggplot2)
require(gridExtra)


#################### Figure 2: Study area map ####################

# libraries
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# define plot region
lon.range <- c(-108, -77)
lat.range <- c(5.5, 23)

# Get plot labels/locations
## Get country data
world <- ne_countries(scale = "medium", returnclass = "sf")
## Get internal points for labels
target_countries <- c("Costa Rica", "Guatemala", "Honduras", "Nicaragua", "Mexico", "Panama", "El Salvador")
country_sf <- world[world$name %in% target_countries, ]
country_map <- data.frame(
  name = c("Mexico", "Guatemala", "Honduras", "El Salvador", "Nicaragua", "Costa Rica", "Panama"),
  short_name = c("MX", "GT", "HN", "SV", "NI", "CR", "PA") 
)
country_points <- merge(country_sf, country_map, by = "name")
country_coords <- st_coordinates(suppressWarnings(st_point_on_surface(country_points$geometry)))
country_points <- cbind(country_points, country_coords)
## Mexican states
states_sf_all <- ne_states(country = "mexico", returnclass = "sf")
target_states <- c("Chiapas", "Oaxaca", "Guerrero", "Michoacán", "Colima", "Jalisco", "Nayarit")
states_sf <- states_sf_all[states_sf_all$name %in% target_states, ]
state_map <- data.frame(
  name = target_states,
  short_state = c("CS", "OA", "GR", "MI", "CL", "JA", "NA")
)
state_points <- merge(states_sf, state_map, by = "name")
state_coords <- st_coordinates(suppressWarnings(st_point_on_surface(state_points$geometry)))
state_points <- cbind(state_points, state_coords)

# Plot
g1 <- ggplot(data = world) +
  geom_sf(fill = "white", color = "black", size = 0.4) +
  geom_sf(data = states_sf, fill = "gray95", color = "gray50", size = 0.2) +
  ## Country Labels with your specific adjustments
  geom_text(data = subset(country_points, !(short_name %in% c("MX", "GT", "SV", "CR", "PA"))), aes(X, Y, label=short_name), fontface = "italic", size = 3.3) +
  geom_text(data = subset(country_points, short_name == "MX"), aes(X, Y, label=short_name), fontface = "italic", size = 3.3, nudge_y = -4, nudge_x = 4) +
  geom_text(data = subset(country_points, short_name == "GT"), aes(X, Y, label=short_name), fontface = "italic", size = 3.3, nudge_y = -0.3) +
  geom_text(data = subset(country_points, short_name == "SV"), aes(X, Y, label=short_name), fontface = "italic", size = 3.3, nudge_x = -0.2, nudge_y = -0.1) +
  geom_text(data = subset(country_points, short_name == "CR"), aes(X, Y, label=short_name), fontface = "italic", size = 3.3, nudge_y = 0.3, nudge_x = -0.25) +
  geom_text(data = subset(country_points, short_name == "PA"), aes(X, Y, label=short_name), fontface = "italic", size = 3.3, nudge_x = 0.4) +
  # State Labels
  geom_text(data = subset(state_points, !(short_state %in% c("OA", "JA"))), aes(X, Y, label=short_state), color = "gray50", fontface = "italic", size = 2.5) +
  geom_text(data = subset(state_points, short_state == "JA"), aes(X, Y, label=short_state), color = "gray50", fontface = "italic", size = 2.5, nudge_x = -0.5, nudge_y = -0.4) +
  geom_text(data = subset(state_points, short_state == "OA"), aes(X, Y, label=short_state), color = "gray50", fontface = "italic", size = 2.5, nudge_y = -0.3) +
  coord_sf(xlim = lon.range, ylim = lat.range) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(fill = NA, color = "black"),
    axis.title = element_blank(),
    axis.text = element_text(size = 8),
    plot.margin = margin(0, 0, 0, 0, "cm")
  ) +
  # add traplocs and all sightings locs
  geom_point(data = subset(sight.casm, occ %in% 2021:2024), aes(x = lon, y = lat), color = "gray65", size = 1) +
  geom_point(data = traplocs, aes(x = lon, y = lat), size = 2, pch = 13) + 
  # add state space limits and population boundaries
  ## SPB/SSSL
  geom_segment(data = subset(model.limits, label=="S limit"), 
               aes(x = lon, y = lat, xend = lon - 2, yend = lat - 2), 
               color = "black", size = 0.8) +
  geom_text(data = subset(model.limits, label=="S limit"), 
            aes(x = lon - 1.2, y = lat - 1.2, label = "SPB/SSSL"), 
            angle = 45, vjust = -0.5, size = 3) +
  ## NPB/NSSL
  geom_segment(data = subset(model.limits, label=="N limit 1"), 
               aes(x = lon, y = lat, xend = lon - 2, yend = lat - 2), 
               color = "black", size = 0.8) +
  geom_segment(data = subset(model.limits, label=="N limit 5"), 
               aes(x = lon, y = lat, xend = lon - 2, yend = lat - 2), 
               color = "black", size = 0.8) +
  annotate("segment", 
           x = subset(model.limits, label=="N limit 1")$lon - 1.7, 
           y = subset(model.limits, label=="N limit 1")$lat - 1.5, 
           xend = subset(model.limits, label=="N limit 5")$lon - 1.1, 
           yend = subset(model.limits, label=="N limit 5")$lat - 1.3, 
           color = "black", size = 0.5, 
           arrow = arrow(ends = "both", length = unit(0.13, "cm"), type = "closed")) +
  annotate("text", 
           x = mean(c(subset(model.limits, label=="N limit 1")$lon, subset(model.limits, label=="N limit 5")$lon)) - 1.3, 
           y = mean(c(subset(model.limits, label=="N limit 1")$lat, subset(model.limits, label=="N limit 5")$lat)) - 1.3, 
           label = "NPB", angle = -40, vjust = -0.5, size = 3) +
  geom_text(data = subset(model.limits, label=="N limit 5"), 
            aes(x = lon - 1, y = lat - 1, label = "NSSL"), 
            angle = 45, vjust = -0.5, size = 3) +
  ## alternate NSSL
  geom_segment(data = subset(model.limits, label=="N limit 6"), 
               aes(x = lon, y = lat, xend = lon - 2, yend = lat - 2), 
               color = "gray50", size = 0.8) +
  geom_text(data = subset(model.limits, label=="N limit 6"), 
            aes(x = lon - 1, y = lat - 1, label = "alt NSSL"), 
            angle = 45, vjust = -0.5, size = 3, color = "gray50")
  
ggsave(plot=g1, filename="Fig2_Map.tiff", 
        width=155, height = 94, units="mm", dpi = 1200, compression = "lzw")


#################### Figure 3: Conceptual diagram ####################

opar <- par()

ac <- c(20,50,80,93,35,65)
b <- seq(0.5,2,length.out=6)
s <- 20
plim <- 90
pmem <- ac <= plim

tiff("Fig3_poplimconcept.tif", width = 2, height = 5, units = "in", compression="lzw", res=1200,
     pointsize = 10)
par(mai = c(0,0,0,0))
plot(0,0, xlim=c(0,10),ylim=c(0,110), bty="n", xlab = "", ylab = "", main = "", 
     type="n", xaxt = "n", yaxt = "n")
lines(x=c(0,10), y=c(5,5), lwd = 1.5)
lines(x=c(0,10), y=c(100,100), lwd = 1.5)
lines(x=c(0,10), y=c(5,5), lwd = 1.5, col=2, lty=2)
lines(x=c(0,10), y=c(plim, plim), lwd = 1.5, col=2)
segments(b, ac-s, b, ac+s, lwd = 1.5, col=1)
segments(b[pmem], ac[pmem]-s, b[pmem], ac[pmem]+s, lwd = 1.5, col=2)
points(b, ac, pch = 21, bg = 1)
points(b[pmem], ac[pmem], pch=21, bg=2)
text(10,102,"Northern SSL", pos=2)
text(10,plim + 2,"Northern PB", pos=2)
text(10,7,"Southern PB and SSL", pos=2)
segments(rep(10,5), c(0.75,5.75,26.75,48.75,65.75), rep(10,5), c(4.25,25.25,47.25,64.25,89.25), lwd = 3, col="gray")
text(10, c(0.75,15,35,55,75) + 2, paste0("Trap ", 1:5), pos=2)
dev.off()

par(opar)


#################### Figure 4: spatiotemporal distribution of captures ####################

library(ggpattern)
library(ggtext)

# calculate halfway points among trap centroids
traplocs <- traplocs %>% 
  mutate(diff_pre = cumdist - lag(cumdist), diff_nxt = lead(cumdist) - cumdist) %>%
  rowwise() %>%
  mutate(
    ymin = if_else(trapid==2, model.limits$cumdist[1], cumdist - (if_else(is.na(diff_pre), 0.5, diff_pre/2))),
    ymax = if_else(trapid==1, model.limits$cumdist[1], cumdist + (if_else(is.na(diff_nxt), 0.5, diff_nxt/2)))) %>%
  ungroup()

sight.casm.2124.ann <- 
  sight.casm %>% 
  filter(occ %in% 2021:2024) %>% 
  group_by(id, occ, trapid) %>% 
  summarize(lat = mean(lat, na.rm=TRUE), lon = mean(lon, na.rm=TRUE)) %>% ungroup() %>% 
  group_by(occ, trapid) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  tidyr::complete(occ, trapid, fill = list(n = 0)) %>% 
  left_join(traplocs %>% select(trapid, cumdist, ymin, ymax)) %>% 
  left_join(eff %>% select(occ, trapid, dist.km)) %>% 
  mutate(logeff = log(dist.km+1))

g1 <- 
  ggplot(data=sight.casm.2124.ann, aes(x=factor(occ), y=cumdist)) +
  
  # heatmap for effort
  geom_rect(data = subset(sight.casm.2124.ann, !is.na(logeff)),
            aes(xmin = as.numeric(factor(occ)) - 0.45, 
                xmax = as.numeric(factor(occ)) + 0.45,
                ymin = ymin, ymax = ymax, fill = logeff), 
            alpha = 0.6) +
  scale_fill_viridis_c(
    option = "mako", 
    direction = -1, 
    limits = c(0, max(sight.casm.2124.ann$logeff, na.rm=TRUE)),
    name = "Log Effort (km)", guide = guide_colorbar(alpha = 0.6)) +
  
  # hatching for NAs
  geom_rect_pattern(data = subset(sight.casm.2124.ann, is.na(logeff)),
                    aes(xmin = as.numeric(factor(occ)) - 0.45, 
                        xmax = as.numeric(factor(occ)) + 0.45,
                        ymin = ymin, ymax = ymax),
                    pattern          = 'stripe',
                    pattern_colour   = "grey80",
                    pattern_spacing  = 0.05, # Very fine
                    pattern_res      = 600,
                    fill             = "white") +
  
  # plot bubble sizes on these axes that are proportional to n
  geom_point(data = subset(sight.casm.2124.ann, n > 0),
             aes(x = factor(occ), y = cumdist, size = n),
             color = "black") +
  # Plus sign for n == 0
  geom_point(data = subset(sight.casm.2124.ann, n == 0),
             aes(x = factor(occ), y = cumdist),
             shape = 3, color = "black", size = 1.5) +
  # bubble legend
  scale_size_continuous(name = "Number of unique IDs",
    breaks = c(0, 10, 50, 100), limits = c(0, 110)) +
  guides(size = guide_legend(override.aes = list(shape = c(3, 19, 19, 19)))) +
  
  labs(x="Year",y="Distance along trap array (100 km)") +

  scale_y_continuous(expand = c(0,0),
                     # Create a secondary axis on the right side with the same mapping as the left axis
                     sec.axis = sec_axis(
                       transform = ~., 
                       breaks = traplocs$cumdist, 
                       labels = c("PA", "CR-S", "CR-N", "NI", "SV", "GT", "OAX-PA", "OAX-PE", "GUE"), 
                       name = ""
                     )) +
   theme_bw() + theme(panel.grid.major= element_line(colour="white"), panel.grid.minor = element_line(colour="white"),
                      text = element_text(size = 10)) +
   coord_cartesian(xlim = c(1,4), ylim = c(-1, 25), clip = "off") +
   theme(plot.margin = unit(c(0.01, 0.11, 0, 0), units="npc")) +
   geom_hline(yintercept = model.limits$cumdist[1], linetype = 2, color = "gray") 
  

ggsave(plot=g1, filename="Fig4_DistCaps_2124.tiff", 
       width=110, height=140, units="mm", dpi = 1200, compression = "lzw")


#################### # Figure 5: Posterior distributions ####################

library(coda)
library(tidyr)
library(dplyr)

# get posterior distributions for p0, effort, and s
attach("nim.SCR1Ddist.partialeff.log.sigreg.casm.2124.QAB.19.rdata")
## s
s <- as.matrix(samples.nim.pleff.sigreg.2124.QAB.19$samples2)[,40:2039] 
z <- as.matrix(samples.nim.pleff.sigreg.2124.QAB.19$samples2)[,2042:4041]
slive <- matrix(NA, dim(s)[1], dim(s)[2])
for (i in 1:dim(z)[1]) {
  live <- as.logical(z[i,])
  slive[i,live] <- s[i,live]
}
df.s <- data.frame(s = as.vector(slive), model = "base")
## p0 and effort
trap_labels <- c("PA", "CR-S", "CR-N", "NI", "SV", "GT", "OAX-PA", "OAX-PE", "GUE")
pd.p <- data.frame(as.matrix(samples.nim.pleff.sigreg.2124.QAB.19$samples)[,40:75]) %>% 
  pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
  extract(param, into = c("trap", "occ"), regex = "\\.(\\d+)\\.\\.(\\d+)\\.", convert = TRUE) %>%
  mutate(type = "p0", trap = factor(trap, levels = 1:9, labels = trap_labels))
pd.e <- data.frame(as.matrix(samples.nim.pleff.sigreg.2124.QAB.19$samples)[,4:39]) %>% 
  pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
  extract(param, into = c("trap", "occ"), regex = "\\.(\\d+)\\.\\.(\\d+)\\.", convert = TRUE) %>%
  mutate(type = "effort") %>% 
  mutate(value = log(value + 1), trap = factor(trap, levels = 1:9, labels = trap_labels))
detach(2)

# plot boxplots for e and p0, and a continuous distribution for s
e_plot <- ggplot(pd.e, aes(x = factor(occ), y = value)) +
  stat_summary(
    fun.data = function(x) {
      q <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
      data.frame(ymin = q[1], lower = q[2], middle = q[3], upper = q[4], ymax = q[5])
    },
    geom = "boxplot", size = 0.2, fill = "gray70") +
  facet_grid(~ trap) +
   scale_x_discrete(name = "Occasion", sec.axis = dup_axis(name = "Trap", breaks = NULL)) +
  scale_y_continuous(name = "Log Effort (km)") +
  theme_bw() + 
  theme(
    panel.grid.major = element_line(colour = "white"), 
    panel.grid.minor = element_line(colour = "white"),
    text = element_text(size = 10),
    strip.background = element_blank(),
    strip.placement = "outside"
  )

p_plot <- ggplot(pd.p, aes(x = factor(occ), y = value)) +
  stat_summary(
    fun.data = function(x) {
      q <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
      data.frame(ymin = q[1], lower = q[2], middle = q[3], upper = q[4], ymax = q[5])
    },
    geom = "boxplot", size = 0.2, fill = "gray70") +
  facet_grid(~ trap) +
  scale_x_discrete(name = "Occasion", sec.axis = dup_axis(name = "Trap", breaks = NULL)) +
  scale_y_continuous(name = expression(italic(p[paste("0,j,k")]))) +
  theme_bw() + 
  theme(
    panel.grid.major = element_line(colour = "white"), 
    panel.grid.minor = element_line(colour = "white"),
    text = element_text(size = 10),
    strip.background = element_blank(),
    strip.placement = "outside"
  )
  
s_plot <- ggplot(df.s) + geom_density(aes(s, color=model), color="black") +
  scale_x_continuous(name = "Cumulative distance along trap array (100 km)") +
  scale_y_continuous(name = expression(paste("Probability density of all ", italic(s[i])))) +
  theme_bw() + 
  theme(
    panel.grid.major = element_line(colour = "white"), 
    panel.grid.minor = element_line(colour = "white"),
    text = element_text(size = 10),
    strip.background = element_blank(),
    strip.placement = "outside"
  )

ggsave(plot=marrangeGrob(grobs = list(e_plot, p_plot, s_plot), nrow = 3, ncol=1, top = NULL), 
       filename="posteriors.tiff", 
       width=160, height=160, units="mm", compression = "lzw")

