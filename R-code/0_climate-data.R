library(raster)
library(sp)


# read in data 
list.files()
dat_abun <- read.csv("Master List of Analysis Studies  - Abundance Stats.csv")
dat_dive <- read.csv("Master List of Analysis Studies  - Diversity Stats.csv")

# make dataframe with lat and long 
xy_abun <- dat_abun[,c("Longitude..E", "Latitude..N")]
xy_dive <- dat_dive[,c("Longitude..E", "Latitude..N")]

# make spatial points dataframelist 
spdf_abun <- SpatialPointsDataFrame(coords = xy_abun, 
                                    data = dat_abun,
                                    proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
spdf_dive <- SpatialPointsDataFrame(coords = xy_dive, 
                                    data = dat_dive,
                                    proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))



# load worldclim data 
r <- raster::getData("worldclim",var="bio",res=10)
r <- r[[c(1,12)]]
names(r) <- c("Temp","Prec")


# temp and precip values at my coordinates 
values_abun <- raster::extract(r,spdf_abun)
values_dive <- raster::extract(r,spdf_dive)


# add temp and precip values to dataframe 
df_abun <- cbind.data.frame(dat_abun,values_abun)
df_dive <- cbind.data.frame(dat_dive,values_dive)


#fix temp to get degrees C 
df_abun$Temp <- df_abun$Temp/10
df_dive$Temp <- df_dive$Temp/10


# export dataframe 
write.csv(df_abun, "Master List of Analysis Studies  - Abundance Stats_climdat.csv")
write.csv(df_dive, "Master List of Analysis Studies  - Diversity Stats_climdat.csv")
