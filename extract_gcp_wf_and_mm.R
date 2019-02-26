 # look at GCP wetland extent model (AGU 2018 version)

library(raster)
library(ncdf4)
library(lubridate)
detach("package:tidyverse", unload=TRUE)

gcp_ch4 <- brick("C:/Users/gavin/Box/MethaneFluxnetSynthesisAnalysis/Data/Upscaling_Analysis/gcp-ch4_wetlands_2000-2018_05deg.nc")

str(gcp_ch4)
crs(gcp_ch4)

sites <- read.csv("C:/Users/gavin/Box/MethaneFluxnetSynthesisAnalysis/Data/Upscaling_Analysis/190201_sites.csv")

# set filtered data as spatial object
sites.geom <- st_as_sf(sites, coords = c("site.long", "site.lat"))
st_crs(sites.geom) <- "+proj=latlong"
sites.geom <- st_transform(sites.geom, crs = " +proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
st_geometry(sites.geom)

site.wetf <- as_tibble(raster::extract(gcp_ch4, sites.geom))

site.wetf <- cbind(sites,site.wetf)
head(site.wetf)

library(tidyverse)
sites.wf <- site.wetf %>% 
  gather(key = Month, value = WetlandF, 8:length(site.wetf)) %>% 
  arrange(site.names, Month)  

sites.wf <- as_tibble(sites.wf)
names(sites.wf) <- c("Index", "ID", "Class", "Biome",
                     "PF","LAT", "LON","Month","WetlandF")

sites.wf <- sites.wf %>% 
  mutate(Month2 = Month) %>% 
  mutate(ID = as.factor(substr(ID,1,5)),
         Date = substr(Month2,2,11),
         Year = substr(Month2,2,5),
         Month = as.factor(as.integer(substr(Month2,7,8)))) %>% 
  select(-Month2) %>% 
  mutate(Date = as.Date(Date, format = "%Y.%m.%d"),
         DOY = yday(Date)) %>% 
  select(ID, Class, Biome, PF, LAT, LON, Year, Month, WetlandF)

sites.wf %>% 
  filter(ID == "ATNeu") %>% 
  ggplot(aes(Month, WetlandF)) +
  geom_point() + 
  facet_wrap(~Year, ncol = 6, scale = 'free')

write.csv(sites.wf,"C:/Users/gavin/Box/MethaneFluxnetSynthesisAnalysis/Data/Upscaling_Analysis/gcp_wetlandfraction_2019.csv" )


## take a look at where the sites are (are there any sites off the map?)
par(mfrow=c(1,1))
plot(gcp_ch4$X2000.08.16)
plot(st_geometry(sites.geom), add = TRUE, pch = 18)

x <- drawExtent()
gcp_ch4_small <- crop(gcp_ch4_small, x)

# yes there are 3/4 sites off the map! e.g. NZKOP
plot(gcp_ch4_small$X2000.08.16)
plot(st_geometry(sites.geom), add = TRUE, pch = 18)

# write out list of sites and amount of WetlandF data per site
sites.wf %>% 
  group_by(ID) %>% 
  summarize(count = sum(!is.na(WetlandF))) %>% 
  write.csv("C:/Users/gavin/Box/MethaneFluxnetSynthesisAnalysis/Data/Upscaling_Analysis/gcp_wetlandf_fluxnetch4_sites.csv")

## extract methane flux and downscale to wetland area

merra2 <- brick("C:/Users/gavin/Box/MethaneFluxnetSynthesisAnalysis/Data/Upscaling_Analysis/LPJ_mmch4e_2000-2017_MERRA2.nc")
cru <- brick("C:/Users/gavin/Box/MethaneFluxnetSynthesisAnalysis/Data/Upscaling_Analysis/LPJ_mmch4e_2000-2017_CRU.nc")

site.merra2 <- as_tibble(raster::extract(merra2, sites.geom))
site.cru <- as_tibble(raster::extract(cru, sites.geom))

site.merra2 <- cbind(sites,site.merra2)
site.cru <- cbind(sites,site.cru)

library(tidyverse)

sites.merra2 <- site.merra2 %>% 
  gather(key = Month, value = merra2_ch4, 8:223) %>% 
  arrange(site.names, Month)  
sites.cru <- site.cru %>% 
  gather(key = Month, value = cru_ch4, 8:223) %>% 
  arrange(site.names, Month)  

sites.merra2 <- as_tibble(sites.merra2)
sites.cru <- as_tibble(sites.cru)
names(sites.merra2) <- c("Index", "ID", "Class", "Biome",
                     "PF","LAT", "LON","Month","merra2_ch4")
names(sites.cru) <- c("Index", "ID", "Class", "Biome",
                         "PF","LAT", "LON","Month","cru_ch4")

sites.wf.mm <- cbind(sites.wf[1:length(sites.merra2$Index),],sites.merra2[,9],sites.cru[,9])

sites.wf.mm.final <- sites.wf.mm %>% 
  mutate(merra2_ch4_dscaled = merra2_ch4/WetlandF,
         cru_ch4_dscaled = cru_ch4/WetlandF,
         DayinMonth = case_when(Month %in% c(1,3,5,7,8,10,12) ~ 31, 
                                Month %in% c(4,6,9,11) ~ 30,
                                Month == 2 ~ 28)) %>% 
  mutate(merra2_ch4_nmolm2s1 = merra2_ch4_dscaled*(1/16.04)*(10^9)*(1/(DayinMonth*24*60*60)),
         cru_ch4_nmolm2s1 = cru_ch4_dscaled*(1/16.04)*(10^9)*(1/(DayinMonth*24*60*60))) %>% 
  select(ID, Year, Month, WetlandF, Merra2 = merra2_ch4_nmolm2s1, CRU = cru_ch4_nmolm2s1)

sites.wf.mm.final %>% 
  filter(ID == "BCBog") %>% 
  ggplot(aes(Month, Merra2)) +
  geom_point() + 
  geom_point(aes(Month, CRU), col = 'pink') +
    facet_wrap(~Year, ncol = 6)


