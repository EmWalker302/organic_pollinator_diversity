#load ggplot2
library(ggplot2)
library(NatParksPalettes)




##### abundance effects 
dat_abun <- read.csv("Model-output4_Single-moderator_effects-by-category_abundance.csv")
dat_abun <- dat_abun[c(1:19),]
dat_abun$moderator <- c("All", rep("Taxonomic.Category", 5), rep("Sampling.Location", 4), rep("Cropping.System", 3), rep("Sampling.Method", 3), rep("Landscape.Complexity", 3))
dat_abun$level <- c("All", "Bees", "Bumblebees", "Butterflies", "Hoverflies", "Solitary Bees", "Across System", "Center", "Edge", "Margin", "Cereal", "Mixed", "Perennial", "Combination", "Pan traps", "Transects", "Complex", "Gradient", "Simple")


##### diversity effects
dat_dive <- read.csv("Model-output4_Single-moderator_effects-by-category_diversity.csv")
dat_dive <- dat_dive[c(1:23),]
dat_dive$moderator <- c("All", rep("Taxonomic.Category", 6), rep("Sampling.Location", 4), rep("Cropping.System", 4), rep("Sampling.Method", 5), rep("Landscape.Complexity", 3))
dat_dive$level <- c("All", "Bees", "Bumblebees", "Butterflies", "Hoverflies", "Moths", "Solitary Bees", "Across System", "Center", "Edge", "Margin", "Cereal", "Mixed", "Pasture", "Perennial", "Combination", "Other", "Pan Traps", "Transects", "Trap Nests", "Complex", "Gradient", "Simple")



mdat <- merge(dat_abun, dat_dive, all=TRUE, by="group", suffixes = c(".abun", ".dive"))
#mdat <- mdat[c(13, 19:24, 1:4, 5:8, 14:18, 9:12),]
mdat <- mdat[c(12, 9:11, 1:4, 18:19, 23, 20:22,5:8, 13, 15:17, 14),]
mdat$index <- as.numeric(c(dim(mdat)[1]:1))

# get colors 
cols <- natparks.pals("SmokyMtns", 14)


# Fig. 2 in manuscript
xminim <- min(na.omit(c(mdat$lci.abun, mdat$lci.dive)))
xmaxim <- max(na.omit(c(mdat$uci.abun, mdat$uci.dive))) + 2

mdat$signif.abun <- rep("", dim(mdat)[1])
mdat$signif.abun[which(mdat$lci.abun>0)] <- "*"

mdat$signif.dive <- rep("", dim(mdat)[1])
mdat$signif.dive[which(mdat$lci.dive>0)] <- "*" 


forestplot <- ggplot() +
  # add dividers
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  geom_hline(yintercept=c(22.5, 19.5, 15.5, 9.5, 5.5), color='black', linetype='solid', alpha=.5) +
  # add diversity responses
  geom_point(data=mdat, col=cols[2], aes(y=index+0.15, x=mean.dive)) +
  geom_errorbarh(data=mdat, height=.2, col=cols[2], aes(y=index+0.15, xmin=lci.dive, xmax=uci.dive)) +
  # add diversity responses
  geom_point(data=mdat, col=cols[6], aes(y=index-0.15, x=mean.abun)) +
  geom_errorbarh(data=mdat, height=.15, col=cols[6], aes(y=index-0.15, xmin=lci.abun, xmax=uci.abun)) +
  # specify axes
  scale_y_continuous(name = "", breaks=c(22:1), expand = c(0, 0),
                        labels=c( "Complex", "Gradient", "Simple",
                                  "Cereal", "Mixed", "Pasture", "Perennial",
                                  "Bees", "Bumble Bees",  "Solitary Bees", "Butterflies", "Hoverflies", "Moths", 
                                  "Across System", "Center", "Edge", "Margin",   
                                  "Combination", "Pan Traps", "Transects", "Trap Nests", "Other"),
                        limits=c(0.5,24)) + 
  annotate("text", x = -1.5, y = 23.25, label = "Diversity", hjust=1, color=cols[2]) +  
  annotate("text", x = -1.5, y = 22.75, label = "Abundance", hjust=1, color=cols[6]) +  
  coord_cartesian(clip = "off", xlim=c(-1.1, 5.1)) +
  scale_x_continuous(breaks = c(-1, 0, 2, 4)) +
  labs(title='', x="Hedges' g", y = '') +
  # set theme
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_line(size = 0.4)) +
  # add labels for moderators
  annotate(geom= "text", x=rep(xmaxim,6), y=c(23.25, 21, 17.5, 12.5, 7.5, 2.5), hjust = 1,vjust=0.5,fontface =3,
           label=c("Overall", "Landscape\nComplexity", "Crop\nType", "Pollinator\nGroup", "Sampling\nLocation", "Sampling\nMethod")) +
  # add sample sizes
  annotate(geom= "text", x=mdat$uci.dive+0.5, y=c(23:1+0.15), 
           hjust = 0, size=2.25, label=paste0(mdat$nobs.dive, " / ", mdat$nstudy.dive," ", mdat$signif.dive)) +
  annotate(geom= "text", x=mdat$uci.abun+0.5, y=c(23:1-0.15), 
           hjust = 0, size=2.25, label=paste0(mdat$nobs.abun, " / ", mdat$nstudy.abun," ", mdat$signif.abun)) 

png("Figures-2_Forestplot.png", width=5, height=7.5, units="in",res=600)
forestplot
dev.off()
