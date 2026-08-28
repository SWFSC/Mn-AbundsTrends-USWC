# Minimum reproducible datasets

# trap information
traps.2124 %<>% 
  left_join(eff %>% select(trapid, Region) %>% distinct()) %>% 
  select(trapid, Region, lat, lon, cumdist)
write.csv(traps.2124, row.names = F, file="Trap information.csv")

# 2021-2024 Spatial capture histories
sight.casm.2124.scr <- 
  sight.casm %>% 
  filter(occ %in% 2021:2024) %>% 
  select(id, occ, trapid) %>% distinct()
write.csv(sight.casm.2124.scr, row.names = F, file="Spatial captures CASM 2021 to 2024.csv")

# 2004-2006 capture histories
sight.cam.0406.ann <- 
  sight.cam %>% 
  filter(occ %in% 2004:2006) %>% 
  select(id, occ) %>% distinct()
write.csv(sight.cam.0406.ann, row.names = F, file="Captures Central America 2004 to 2006.csv")

# sex information
sight.casm.sex %<>% 
  filter(occ %in% 2021:2024) %>% 
  select(id, sex, source) %>% distinct() %>% 
  mutate(source=if_else(source=="biop","genetic",source))
write.csv(sight.casm.sex, row.names = F, file="Individual sex information.csv")