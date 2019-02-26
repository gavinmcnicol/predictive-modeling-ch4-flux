 # look at GCP wetland extent model (AGU 2018 version)

library(raster)
library(ncdf4)
library(lubridate)
detach("package:tidyverse", unload=TRUE)

gcp_ch4 <- brick("C:/Users/gavin/Box/MethaneFluxnetSynthesisAnalysis/Data/Upscaling_Analysis/gcp-ch4_wetlands_2000-2017_025deg.nc")

str(gcp_ch4)
crs(gcp_ch4)

sites <- read.csv("C:/Users/gavin/Box/MethaneFluxnetSynthesisAnalysis/Data/Upscaling_Analysis/190201_sites.csv")

# set filtered data as spatial object
sites.geom <- st_as_sf(sites, coords = c("site.long", "site.lat"))
st_crs(sites.geom) <- "+proj=latlong"
sites.geom <- st_transform(sites.geom, crs = " +proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
st_geometry(sites.geom)

site.wetf <- as_tibble(extract(gcp_ch4, sites.geom))

site.wetf <- cbind(sites,site.wetf)
head(site.wetf)

library(tidyverse)
sites.wf <- site.wetf %>% 
  gather(key = Month, value = WetlandF, 8:223) %>% 
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

write.csv(sites.wf,"C:/Users/gavin/Box/MethaneFluxnetSynthesisAnalysis/Data/Upscaling_Analysis/gcp_wetlandfraction.csv" )
