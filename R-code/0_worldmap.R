library(tidyverse)
library(gridExtra)
library(rnaturalearth)
library(tmap)
library(sf) 
library(shinyjs)
library(sjmisc)
library(ggplot2)
library(ggmap)


# read data 
dat_abun <- read.csv("Master List of Analysis Studies  - Abundance Stats_climdat.csv")
dat_dive <- read.csv("Master List of Analysis Studies  - Diversity Stats_climdat.csv")


# get unique locations 
dat_abun$lat.long <- paste(dat_abun$Latitude..N, dat_abun$Longitude..E)
length(unique(dat_abun$lat.long))
# there are 23 unique locations for abun dat 

dat_dive$lat.long <- paste(dat_dive$Latitude..N, dat_dive$Longitude..E)
length(unique(dat_dive$lat.long))
# there are 36 unique locations for dive dat 


# make dataframe for locations and add info like number of observations 
locdat_abun <- data.frame(lat.long = unique(dat_abun$lat.long),
                          lat = rep(NA, length(unique(dat_abun$lat.long))),
                          long = rep(NA, length(unique(dat_abun$lat.long))),
                          nobs = rep(NA, length(unique(dat_abun$lat.long))),
                          nstudy = rep(NA, length(unique(dat_abun$lat.long))),
                          StudyID = rep(NA, length(unique(dat_abun$lat.long))))
for(i in 1:dim(locdat_abun)[1]){
  locdat_abun$lat[i] <- strsplit(locdat_abun$lat.long[i], " ")[[1]][1]
  locdat_abun$long[i] <- strsplit(locdat_abun$lat.long[i], " ")[[1]][2]
  locdat_abun$nobs[i] <- length(which(dat_abun$Longitude..E==locdat_abun$long[i] & dat_abun$Latitude..N==locdat_abun$lat[i]))
  study <- unique(dat_abun$Reference.ID[which(dat_abun$Longitude..E==locdat_abun$long[i] & dat_abun$Latitude..N==locdat_abun$lat[i])])
  locdat_abun$nstudy[i] <- length(study)
  if(length(study)==1){locdat_abun$StudyID[i] <- study[1]}
  if(length(study)==2){locdat_abun$StudyID[i] <- paste(study[1], study[2], sep=",")}
  if(length(study)==3){locdat_abun$StudyID[i] <- paste(study[1], study[2], study[3],sep=",")}
}
locdat_abun$lat <- as.numeric(locdat_abun$lat)
locdat_abun$long <- as.numeric(locdat_abun$long)

locdat_dive <- data.frame(lat.long = unique(dat_dive$lat.long),
                          lat = rep(NA, length(unique(dat_dive$lat.long))),
                          long = rep(NA, length(unique(dat_dive$lat.long))),
                          nobs = rep(NA, length(unique(dat_dive$lat.long))),
                          nstudy = rep(NA, length(unique(dat_dive$lat.long))),
                          StudyID = rep(NA, length(unique(dat_dive$lat.long))))
for(i in 1:dim(locdat_dive)[1]){
  locdat_dive$lat[i] <- strsplit(locdat_dive$lat.long[i], " ")[[1]][1]
  locdat_dive$long[i] <- strsplit(locdat_dive$lat.long[i], " ")[[1]][2]
  locdat_dive$nobs[i] <- length(which(dat_dive$Longitude..E==locdat_dive$long[i] & dat_dive$Latitude..N==locdat_dive$lat[i]))
  study <- unique(dat_dive$Reference.ID[which(dat_dive$Longitude..E==locdat_dive$long[i] & dat_dive$Latitude..N==locdat_dive$lat[i])])
  locdat_dive$nstudy[i] <- length(study)
  if(length(study)==1){locdat_dive$StudyID[i] <- study[1]}
  if(length(study)==2){locdat_dive$StudyID[i] <- paste(study[1], study[2], sep=",")}
  if(length(study)==3){locdat_dive$StudyID[i] <- paste(study[1], study[2], study[3],sep=",")}
}
locdat_dive$lat <- as.numeric(locdat_dive$lat)
locdat_dive$long <- as.numeric(locdat_dive$long)


# get world map 
world = ne_countries(returnclass = "sf") 
countries <- unique(world$name_long)
states <- unique(map_data("state")$region)


# create map for number of studies 
sf_use_s2(FALSE)
my_map_abun <- ggplot() +
  geom_sf(data = world, inherit.aes = F) +
  geom_sf(data = world, fill="white", alpha=0.2, colour = "gray75") +
  geom_point(data = locdat_abun, aes(x = long, y = lat, size=nstudy),  shape = 21, fill = "turquoise3", alpha=0.5) +
  theme(legend.position = c(0.8, 0.2)) +
  labs(x="Longitude", y="Latitude") +
  theme_bw() +
  coord_sf(expand = F) +
  scale_x_continuous(breaks = seq(-180, 180, by=60)) +
  scale_y_continuous(breaks = seq(-90, 90, by = 30)) +
  scale_size_continuous(name = "Studies", breaks = c(1,2,3), range = c(1,3))
my_map_abun

my_map_dive <- ggplot() +
  geom_sf(data = world, inherit.aes = F) +
  geom_sf(data = world, fill="white", alpha=0.2, colour = "gray75") +
  geom_point(data = locdat_dive, aes(x = long, y = lat, size=nstudy),  shape = 21, fill = "turquoise3", alpha=0.5) +
  theme(legend.position = c(0.8, 0.2)) +
  labs(x="Longitude", y="Latitude") +
  theme_bw() +
  coord_sf(expand = F) +
  scale_x_continuous(breaks = seq(-180, 180, by=60)) +
  scale_y_continuous(breaks = seq(-90, 90, by = 30)) +
  scale_size_continuous(name = "Studies", breaks = c(1,2,3), range = c(1,3))
my_map_dive


my_map <- ggpubr::ggarrange(my_map_dive, my_map_abun, 
                            label.x = -0.125, label.y=1.075,
                            labels=c("A) Pollinator diversity", "B) Pollinator abundance"),
                            ncol=1, common.legend=TRUE, legend = "bottom") +
  theme(plot.margin = margin(0.75,0,0,0, "cm")) 


png("WorldMap_nstudy.png", height=4000, width=3600,  units="px", res=600)
my_map
dev.off()


write.csv(dat_abun, "Master List of Analysis Studies  - Abundance Stats_climdat2.csv")
write.csv(dat_dive, "Master List of Analysis Studies  - Diversity Stats_climdat2.csv")

