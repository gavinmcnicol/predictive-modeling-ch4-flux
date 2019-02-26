## Script for extracting WorldClim variables
# Gavin McNicol
# Feb 2019

## load Bio1 (MAT), Bio3 (Isothermality (Diurnal range/Annual Range)), 
##      Bio4 (Temperature Seasonality (SD*100)), Bio12 (MAP), Bio15 (Precipitation Seasonality (CV))

## load site list for worldclim (ID, Category, LAT ,LONG)
sites <- read_csv("C:/Users/gavin/Box/MODIS Data/AppEEARS Data/V2.0_sites.csv")
site.coords <- cbind(sites$Longitude,sites$Latitude)


## spatial libraries
library(raster)
library(sp)
library(sf)
library(rgdal)

# load rasters
Bio1 <- raster("C:/Users/gavin/OneDrive/Desktop/wc2.0_30s_bio/wc2.0_bio_30s_01.tif")
Bio3 <- raster("C:/Users/gavin/OneDrive/Desktop/wc2.0_30s_bio/wc2.0_bio_30s_03.tif")
Bio4 <- raster("C:/Users/gavin/OneDrive/Desktop/wc2.0_30s_bio/wc2.0_bio_30s_04.tif")
Bio12 <- raster("C:/Users/gavin/OneDrive/Desktop/wc2.0_30s_bio/wc2.0_bio_30s_12.tif")
Bio15 <- raster("C:/Users/gavin/OneDrive/Desktop/wc2.0_30s_bio/wc2.0_bio_30s_15.tif")

# create climate stack
clim <- stack(Bio1, Bio3, Bio4, Bio12, Bio15)

plot(Bio4)


# extract at points
clim.data <- extract(clim, site.coords)
clim.data

# bind final values
sites.clim.data <- as_tibble(cbind(sites, clim.data))
names(sites.clim.data) <- c("ID","Class","LAT","LONG","Bio1","Bio3","Bio4","Bio12",'Bio15')

# write csv
write.csv(sites.clim.data, "C:/Users/gavin/Box/Upscaling Resources/WorldClim Data/V2.0_sites_worldclim.csv")
