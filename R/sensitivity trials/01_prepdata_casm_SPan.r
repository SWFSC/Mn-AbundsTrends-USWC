# Prepare data for sensitivity trial for southern boundary location (assess associated difference from Curtis et al 2022)
# Author: K. Alexandra Curtis <alex.curtis@noaa.gov>
# adjust model.limit call in run script to stay at lat=19.2 for N limit

############ CREATE CENTRAL AMERICA / SOUTHERN MEXICO DATA SETS ############

# load libraries
# pak::pkg_install("kacurtis/CMRutils")
library(CMRutils)
library(dplyr)

# import data
mnid <- import_crcmnid(datdir = "../../Humpbacks/data/external/CRC/", mnidfile = "MN ID 20260127.mdb")

# initial data filtering, editing, and selection
source("selectfields_crcdata.r")

# edit all sightings
sight %<>% 
  mutate(year = if_else(is.na(date), year, as.integer(year(date)))) %>% 
  # drop any IDs seen in southern hemisphere
  group_by(id) %>% 
  mutate(sohem = any(subarea=="Peru-Antarc")) %>%
  ungroup() %>% 
  filter(is.na(sohem) | !sohem) %>% select(-sohem) %>% 
  # eliminate any other sightings south of Panama
  filter(is.na(lat) | lat>7.2) %>% 
  # correct any negative longitudes after filtering out Japan and Russia
  filter(!(loc %in% c("6","91","97"))) %>% 
  mutate(lon = if_else(lon>0, -1*lon, lon))

# subset wintering area sightings
sight.win <- sight %>% 
  # filter winter 
  filter(!(mm %in% 6:10)) %>% 
  # create occ to shift Nov/Dec into next year (if no month, occ=year)
  mutate(occ = if_else(mm %in% 11:12, year+1, as.numeric(year))) %>% 
  # subset wintering area regions
  filter(lat < 32.2) %>% 
  # eliminate new calf IDs for that year from wintering areas
  group_by(id, occ) %>% 
  mutate(calf = any(grepl("C", AGEREP))) %>% 
  ungroup() %>% filter(!calf) %>% select(-calf)
  
# filter and edit sight for CAm/SMex 
sight.casm <- sight.win %>%
  # filter Central America and S Mexico sightings (note "S Mexico" includes lat<19)
  filter((lat>7.2 & lat<18 & lon>(-102.2) & lon<(-75))) %>% # |  
         # this line does not currently add any records
         # (is.na(lat) & rgn %in% c("Central America", "S Mexico"))) %>% 
  # drop Dead animals
  filter(!Dead) %>%
  # filter out is.na(occ) (ids for records without occ all have other CASM records)
  filter(!is.na(occ)) %>% 
  # add field indicating whether Record Source is HW
  mutate(hw=grepl("HW",Record.Source)) %>% 
  # filter out poor quality (Q=3,4 or PQ, idq=1 or 2 without mitigating Q from double scoring yrs, or HW source and idq=0)
  filter((!hw & (is.na(Quality) | !(Quality %in% c("3","4","5","PQ","3 PQ","3 pq")))) | 
           (hw & ((is.na(HW_media_idq) & (is.na(Quality) | Quality %in% c("A","B",1:2))) | 
                    ((HW_media_idq %in% 1:2) & Quality %in% c("A","B",1:2)) | 
                    HW_media_idq %in% 3:5))) %>% 
  # group by id and date for filtering and create variables with which to identify problems
  group_by(id, occ, date) %>%
  # # (no NA lats in current data) drop records with NA positions on days with known-position record
  # mutate(nidd=n()) %>% 
  # filter(nidd==1 | all(is.na(lat)) | !is.na(lat)) %>%
  # replace incorrect positions for Nicaragua with more accurate
  mutate(flag = round(lat,5)==11.06329,
         lat = if_else(flag, 12.73, lat),
         lon = if_else(flag, -87.48, lon)) %>% 
  select(-flag) %>% 
  # drop incorrect position for Nicaragua where better available
  filter(CRC.Record!=111836) %>% 
  # keep only "PRECISE" and "Precise" when "GENERAL" or "General" or "Estimated" in same day
  mutate(postype = any(Position.Type %in% c("GENERAL","General","Estimated")) & any(Position.Type %in% c("PRECISE","Precise"))) %>%
  filter((postype & Position.Type %in% c("PRECISE","Precise")) | !postype) %>% 
  ungroup() %>% 
  # # calculate records per day, calculate max lat difference per day
  # group_by(id, occ, date) %>%
  # mutate(nidd=n(), 
  #        maxddlat = if(sum(!is.na(lat))>1) max(dist(lat),na.rm=TRUE) else NA_real_,
  #        maxddlon = if(sum(!is.na(lon))>1) max(dist(lon),na.rm=TRUE) else NA_real_) %>% 
  # select(-c(nidd, maxddlat, maxddlon)) %>% 
  # assign trapids
  # sight.casm %>% ggplot(aes(occ, lat)) + geom_point()
  mutate(trapid = case_when(is.na(lat) ~ NA_character_,
                            lat <= 8.25 ~ "1",   # Panama
                            lat>8.25 & lat<=10 ~ "2",   # South Costa Rica
                            lat>10 & lat<=11.03 ~ "3",   # North Coast Rica; 11.03 serves as de facto boundary for both SPLASH 1 and SPLASH 2 eras
                            lat>11.03 & lat<=12 ~ "4",   # Nicaragua
                            # lat>12 & lat<=13.25 ~ "5",   # S El Salvador/N Nicaragua (only two IDs in 2022 -> collapse with trap 6 for this analysis)
                            lat>12 & lat<=13.7 ~ "5",   # El Salvador
                            lat>13.7 & lat<=14.5 ~ "6",   # Guatemala
                            # longitude split only with coastline distance parameterization of space
                            lat>14.5 & lat<=16.3 & lon>=(-96.7) ~ "7",   # Puerto Angel, Oaxaca
                            lat>14.5 & lat<=16.3 & lon<(-96.7) ~ "8",   # Puerto Escondido, Oaxaca
                            lat>16.3 ~ "9")) %>%   # Guerrero
  # collapse sightings to daily to eliminate remaining duplicates and avoid overweighting particular geographic coordinates
  # (note that loc is not necessarily reliable, e.g. Guerrero assigned both 18 and 19)
  group_by(id, occ, trapid, date) %>%
  summarize(lat = if_else(any(!is.na(lat)), mean(lat, na.rm=T), NA_real_),
            lon = if_else(any(!is.na(lon)), mean(lon, na.rm=T), NA_real_)) %>% 
  ungroup() %>% 
  # calculate trap coordinates for focal years (2021-2024)
  group_by(trapid) %>% 
  mutate(occ.2124 = occ %in% 2021:2024,
         traplat = mean(lat[occ.2124], na.rm=T),
         traplon = mean(lon[occ.2124], na.rm=T)) %>% 
  ungroup()

# subset CAm only
sight.cam <- sight.casm %>% filter(lat < 14.5)

# create spatial annualized capture histories
## trap information
traplocs <- sight.casm %>% select(trapid, traplat, traplon) %>% 
  distinct() %>% arrange(trapid) %>% 
  rename(lat=traplat, lon=traplon)
# traplocs <- sort(unique(sight.casm$traplat))
## add population and model boundaries
model.limits <- data.frame(label=c("alt S limit", "S limit", "N limit 1", "N limit 2", "N limit 3", "N limit 4", "N limit 5", "Jalisco corner", "N limit 6"),
                     lat = c(7.2, 8.0, 17.9, 18.3, 18.6, 19, 19.2, 20.4, 20.8),  
                     lon = c(-80.9, -82.9, -102.2, -103.5, -103.7, -104.4, -104.8, -105.7, -105.5))
### get straight-line distance in 100-km units among the trap and limit coordinates
library(geodist)
coords.model <- traplocs %>% select(lat, lon) %>% 
  bind_rows(model.limits %>% select(lat, lon)) %>% 
  arrange(lat)
coords.dist <- geodist::geodist(coords.model[,c("lat","lon")], measure="haversine")/100000
#### subset sequential distances among traps
coords.dist <- diag(coords.dist[2:dim(coords.dist)[1], 1:(dim(coords.dist)[1]-1)])
### get cumulative distances
coords.model$cumdist <- c(0, cumsum(coords.dist))
### merge back into traplocs and limits data frames
traplocs <- traplocs %>% left_join(coords.model)
model.limits <- model.limits %>% left_join(coords.model)
## spatial capture histories
y3d.casm <- array2binom(unclass(with(sight.casm, table(id, trapid, occ))))
y3d.cam <- array2binom(unclass(with(sight.cam, table(id, trapid, occ))))

# closed population data structures
# (trap information and spatial capture histories for each set of years)
## 2021-2024
### CASM
yy <- 2021:2024
traps.2124 <- sight.casm %>% filter(occ %in% yy) %>%
  select(trapid) %>% distinct() %>% arrange(trapid) %>% 
  left_join(traplocs)
traplocs.2124 <- array(traps.2124$cumdist, dim=dim(traps.2124)[1], dimnames=list(traps.2124$trapid))
y3d.casm.2124 <- dim_drop0(dim_drop0(y3d.casm[,,dimnames(y3d.casm)$occ %in% yy], 1), 2)
### CAm (no trap info needed)
y3d.cam.2124 <- dim_drop0(dim_drop0(y3d.cam[,,dimnames(y3d.cam)$occ %in% yy], 1), 2)
## 2004-2006 (SPLASH 1) - only interested in Central America (SMX only sampled opportunistically in passing)
## (no trap information needed)
yy <- 2004:2006
y3d.cam.0406 <- dim_drop0(dim_drop0(y3d.cam[,,dimnames(y3d.cam)$occ %in% yy], 1), 2)
rm(yy)

# import and arrange effort
eff <- readr::read_csv("Effort data for Central America.csv", col_types = "cciciddiddc")
eff <- eff %>% 
  arrange(occ, trapid) %>% 
  filter(trapid!="10" & occ %in% 2021:2024) %>% 
  select(-(days.ded:dist.ded)) %>% 
  rename(days=days.tot, hrs=hrs.tot, dist.km=dist.tot)
