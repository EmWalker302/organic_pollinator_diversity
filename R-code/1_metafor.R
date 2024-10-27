# Use metafor function rma.mv to model effects of organic farming
# on abundance and diversity of pollinators while allowing for variation at three levels:  
# 1. within observations, 
# 2. among observations of a single study, and 
# 3. among studies

# Our hypotheses:
# 1) Organic farms will have an overall higher species richness and abundance of pollinators compared to conventional farms
# 2) Center samplings will exhibit the most drastic difference due to increased isolation from nearby wild floral resources
# 3) Bees, as obligate pollinators, will experience the largest shifts in diversity and abundance across farming systems 
# 4) Cereal crop fields will benefit most from organic farming due to increased importance of non-crop resources
# 5) Transect sampling will capture the largest difference due to the absence of a lure attracting pollinators from outside of fields 
# 6) Farms in simple landscapes offering less semi-natural habitat and forage will exhibit the largest biodiversity shifts 


library(metafor)

#read data 
dat_abun <- read.csv("Master List of Analysis Studies  - Abundance Stats_climdat2.csv")
dat_dive <- read.csv("Master List of Analysis Studies  - Diversity Stats_climdat2.csv")


##############################
### moderators 

## abundance:

# taxonomic group 
dat_abun$Taxonomic.Category <- as.factor(dat_abun$Taxonomic.Category)
levels(dat_abun$Taxonomic.Category)

# sampling location 
dat_abun$Edge.Center.Margin <- as.factor(dat_abun$Edge.Center.Margin)
levels(dat_abun$Edge.Center.Margin)

# crop type
dat_abun$Cropping.System <- as.factor(dat_abun$Cropping.System)
levels(dat_abun$Cropping.System)

# sampling method
dat_abun$Sampling.Method <- as.factor(dat_abun$Sampling.Method)
levels(dat_abun$Sampling.Method)

# landscape complexity
dat_abun$Landscape.Complexity <- as.factor(dat_abun$Landscape.Complexity)
levels(dat_abun$Landscape.Complexity)


## diversity:

# taxonomic group
dat_dive$Taxonomic.Category <- as.factor(dat_dive$Taxonomic.Category)
levels(dat_dive$Taxonomic.Category)

# sampling location
dat_dive$Edge.Center.Margin <- as.factor(dat_dive$Edge.Center.Margin)
levels(dat_dive$Edge.Center.Margin)

# crop type
dat_dive$Cropping.System <- as.factor(dat_dive$Cropping.System)
levels(dat_dive$Cropping.System)

# sampling method
dat_dive$Sampling.Method <- as.factor(dat_dive$Sampling.Method)
levels(dat_dive$Sampling.Method)

# lanscape complexity
dat_dive$Landscape.Complexity <- as.factor(dat_dive$Landscape.Complexity)
levels(dat_dive$Landscape.Complexity)



### give each study a unique identifier 
StudyIDs_abun <- data.frame(StudyID = unique(dat_abun$Reference.ID),
                            StudyNum = seq(1:length(unique(dat_abun$Reference.ID))))
dat_abun$StudyNum <- rep(NA, dim(dat_abun)[1])
for(i in 1:dim(dat_abun)[1]){
  dat_abun$StudyNum[i] <- StudyIDs_abun$StudyNum[which(StudyIDs_abun$StudyID==dat_abun$Reference.ID[i])]
}
write.csv(StudyIDs_abun, "0_StudyIDs_abun.csv")

StudyIDs_dive <- data.frame(StudyID = unique(dat_dive$Reference.ID),
                            StudyNum = seq(1:length(unique(dat_dive$Reference.ID))))
dat_dive$StudyNum <- rep(NA, dim(dat_dive)[1])
for(i in 1:dim(dat_dive)[1]){
  dat_dive$StudyNum[i] <- StudyIDs_dive$StudyNum[which(StudyIDs_dive$StudyID==dat_dive$Reference.ID[i])]
}
write.csv(StudyIDs_dive, "0_StudyIDs_dive.csv")



# give each observation its own ID from 1 to n observations 
dat_abun$Obs.ID <- seq(1, dim(dat_abun)[1])
dat_dive$Obs.ID <- seq(1, dim(dat_dive)[1])

##### calculate hedges g
# hedges g is equal to: 
# (x1 – x2) / √((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n2-2)
# where x's represent sample means
# n's represent sample size
# s's are standard deviations
# Hedges’ g is a useful effect size in ecological meta-analysis because it statistically corrects for variance that may be introduced when sample sizes are small (Hedges 1981).


SMD_effect_sizes <-
  escalc( # This is the function in metafor that allows us to calculate an effect size for each row in a database
    "SMD",
    # Specify the effect size we want to calculate. In this case SMD for Standardized mean difference
    m1i = Org..Mean,
    # mean richness at invaded sites
    n1i = Org..n,
    # invaded site sample size
    sd1i = Org..SD,
    # invaded site SD
    m2i = Con..Mean,
    # mean richness at control sites
    n2i = Con..n,
    # control site sample size
    sd2i = Con..SD,
    # control site SD
    data = dat_abun # This is where the escalc function can find all the data for our meta-analysis
  )
dat_abun$hedges <- SMD_effect_sizes$yi
dat_abun$variance <- SMD_effect_sizes$vi


SMD_effect_sizes <-
  escalc( # This is the function in metafor that allows us to calculate an effect size for each row in a database
    "SMD",
    # Specify the effect size we want to calculate. In this case SMD for Standardized mean difference
    m1i = Org..Mean,
    # mean richness at invaded sites
    n1i = Org..n,
    # invaded site sample size
    sd1i = Org..SD,
    # invaded site SD
    m2i = Con..Mean,
    # mean richness at control sites
    n2i = Con..n,
    # control site sample size
    sd2i = Con..SD,
    # control site SD
    data = dat_dive # This is where the escalc function can find all the data for our meta-analysis
  )

dat_dive$hedges <- SMD_effect_sizes$yi
dat_dive$variance <- SMD_effect_sizes$vi




# export dataset used for analysis 
write.csv(dat_abun, "dat_abun.csv")
write.csv(dat_dive, "dat_dive.csv")



#### overall models
#### testing whether within-study variation is significant
#### estimating distribution of variance across the three levels
# note that sigma^2.1 is the variance between effect sizes within studies
# and sigma^2.2 is the variance between studies

r <- 16
effectsdat_abun <- data.frame(response=rep(NA,r),
                              group=rep(NA,r),
                              mean=rep(NA,r),
                              lci=rep(NA,r),
                              uci=rep(NA,r),
                              nobs=rep(NA,r),
                              nstudy=rep(NA,r))
effectsdat_dive <- data.frame(responses=rep(NA,r),
                              group=rep(NA,r),
                              mean=rep(NA,r),
                              lci=rep(NA,r),
                              uci=rep(NA,r),
                              nobs=rep(NA,r),
                              nstudy=rep(NA,r))


# overall abundance model
mod <- rma.mv(hedges, variance, 
              random = list(~ 1 | Obs.ID, ~ 1 | StudyNum), 
              tdist=TRUE, data=dat_abun)
effectsdat_abun[1,] <- c("abundance","Overall", mod$b, mod$ci.lb, mod$ci.ub, mod$k,length(unique(dat_abun$Reference.ID)))
# test for significance of variation within- and between-studies
mod_2 <- rma.mv(hedges, variance, random = list(~ 1 | Obs.ID, ~ 1 | StudyNum), tdist=TRUE, data=dat_abun, sigma2=c(0,NA))
mod_3 <- rma.mv(hedges, variance, random = list(~ 1 | Obs.ID, ~ 1 | StudyNum), tdist=TRUE, data=dat_abun, sigma2=c(NA,0))
# variance partitioning
sum.inverse.variances <- sum(1 / (dat_abun$variance))
list.inverse.variances.square <- 1 / (dat_abun$variance^2)
numerator <- (length(dat_abun$variance) - 1) * sum.inverse.variances
denominator <- (sum.inverse.variances) ^ 2 - sum(list.inverse.variances.square)
estimated.sampling.variance <- numerator / denominator
I2_1 <- (estimated.sampling.variance) / (mod$sigma2[1] + mod$sigma2[2] + estimated.sampling.variance)
I2_2 <- (mod$sigma2[1]) / (mod$sigma2[1] + mod$sigma2[2] + estimated.sampling.variance)
I2_3 <- (mod$sigma2[2]) / (mod$sigma2[1]  + mod$sigma2[2] + estimated.sampling.variance)

# write files
sink(file = "Model-output1_Overall_abundance.txt"); summary(mod); sink(file = NULL)
sink(file = "Model-output2_Significance-of-variation_abundance_within-study.txt"); anova(mod, mod_2); sink(file = NULL) # if p<0.05, within-study variance is significant
sink(file = "Model-output2_Significance-of-variation_abundance_between-study.txt"); anova(mod, mod_3); sink(file = NULL) # if p<0.05, between-study variance is significant
sink(file = "Model-output3_Partitioning-of-variation_abundance.txt"); data.frame(samplingvar = I2_1 * 100, withinstudyvar = I2_2 * 100, betweenstudyvar = I2_3 * 100); sink(file = NULL)




# Overall diversity model
mod <- rma.mv(hedges, variance, 
              random = list(~ 1 | Obs.ID, ~ 1 | StudyNum), 
              tdist=TRUE, data=dat_dive)
effectsdat_dive[1,] <- c("diversity","Overall", mod$b, mod$ci.lb, mod$ci.ub, mod$k,length(unique(dat_dive$Reference.ID)))
# test for significance of variation within- and between-studies
mod_2 <- rma.mv(hedges, variance, random = list(~ 1 | Obs.ID, ~ 1 | StudyNum), tdist=TRUE, data=dat_dive, sigma2=c(0,NA))
mod_3 <- rma.mv(hedges, variance, random = list(~ 1 | Obs.ID, ~ 1 | StudyNum), tdist=TRUE, data=dat_dive, sigma2=c(NA,0))
# variance partitioning
sum.inverse.variances <- sum(1 / (dat_dive$variance))
list.inverse.variances.square <- 1 / (dat_dive$variance^2)
numerator <- (length(dat_dive$variance) - 1) * sum.inverse.variances
denominator <- (sum.inverse.variances) ^ 2 - sum(list.inverse.variances.square)
estimated.sampling.variance <- numerator / denominator
I2_1 <- (estimated.sampling.variance) / (mod$sigma2[1] + mod$sigma2[2] + estimated.sampling.variance)
I2_2 <- (mod$sigma2[1]) / (mod$sigma2[1] + mod$sigma2[2] + estimated.sampling.variance)
I2_3 <- (mod$sigma2[2]) / (mod$sigma2[1]  + mod$sigma2[2] + estimated.sampling.variance)
# write files
sink(file = "Model-output1_Overall_diversity.txt"); summary(mod); sink(file = NULL)
sink(file = "Model-output2_Significance-of-variation_diversity_within-study.txt"); anova(mod, mod_2); sink(file = NULL) # if p<0.05, within-study variance is significant
sink(file = "Model-output2_Significance-of-variation_diversity_between-study.txt"); anova(mod, mod_3); sink(file = NULL) # if p<0.05, between-study variance is significant
sink(file = "Model-output3_Partitioning-of-variation_diversity.txt"); data.frame(samplingvar = I2_1 * 100, withinstudyvar = I2_2 * 100, betweenstudyvar = I2_3 * 100); sink(file = NULL)







#### Testing for moderator effects
# Note that you need "mods= ~0 + moderators" for categorical variables
# Note also that order does not matter when adding dummy-coded variables to the model
# Note also that we exclude groups that are represented by less than three studies


### start by analyzing abundance responses

## Taxonomic.Category
levels(dat_abun$Taxonomic.Category)
taxlev <- c("Bees", "Bumblebees", "Butterflies", "Hoverflies", "Solitary Bees") ##categories represented by at least 3 studies 
mod <- rma.mv(hedges, variance, mods= ~0+Taxonomic.Category, 
              random = list(~ 1 | Obs.ID, ~ 1 | StudyNum), tdist=TRUE, 
              data=dat_abun[which(dat_abun$Taxonomic.Category %in% c("Bees", "Bumblebees", "Butterflies", "Hoverflies", "Solitary Bees")),])
sink(file = "Model-output4_Single-moderator_abundance_Taxonomic.Category.txt"); summary(mod); sink(file = NULL)
mod
effectsdat_abun[2,] <- c("abundance",rownames(mod$b)[1], mod$b[1,1], mod$ci.lb[1], mod$ci.ub[1], length(which(dat_abun$Taxonomic.Category==taxlev[1])), length(unique(dat_abun$Reference.ID[which(dat_abun$Taxonomic.Category==taxlev[1])])))
effectsdat_abun[3,] <- c("abundance",rownames(mod$b)[2], mod$b[2,1], mod$ci.lb[2], mod$ci.ub[2], length(which(dat_abun$Taxonomic.Category==taxlev[2])), length(unique(dat_abun$Reference.ID[which(dat_abun$Taxonomic.Category==taxlev[2])])))
effectsdat_abun[4,] <- c("abundance",rownames(mod$b)[3], mod$b[3,1], mod$ci.lb[3], mod$ci.ub[3], length(which(dat_abun$Taxonomic.Category==taxlev[3])), length(unique(dat_abun$Reference.ID[which(dat_abun$Taxonomic.Category==taxlev[3])])))
effectsdat_abun[5,] <- c("abundance",rownames(mod$b)[4], mod$b[4,1], mod$ci.lb[4], mod$ci.ub[4], length(which(dat_abun$Taxonomic.Category==taxlev[4])), length(unique(dat_abun$Reference.ID[which(dat_abun$Taxonomic.Category==taxlev[4])])))
effectsdat_abun[6,] <- c("abundance",rownames(mod$b)[5], mod$b[5,1], mod$ci.lb[5], mod$ci.ub[5], length(which(dat_abun$Taxonomic.Category==taxlev[5])), length(unique(dat_abun$Reference.ID[which(dat_abun$Taxonomic.Category==taxlev[5])])))

## Edge.Center.Margin
levels(dat_abun$Edge.Center.Margin)
taxlev <- c("Across System", "Center", "Edge", "Margin") ##categories represented by at least 3 studies
mod <- rma.mv(hedges, variance, mods= ~0+Edge.Center.Margin, 
              random = list(~ 1 | Obs.ID, ~ 1 | StudyNum), tdist=TRUE, data=dat_abun[which(dat_abun$Edge.Center.Margin %in% c("Across System", "Center", "Edge", "Margin")),])
sink(file = "Model-output4_Single-moderator_abundance_Edge.Center.Margin.txt"); summary(mod); sink(file = NULL)
mod
effectsdat_abun[7,] <- c("abundance",rownames(mod$b)[1], mod$b[1,1], mod$ci.lb[1], mod$ci.ub[1], length(which(dat_abun$Edge.Center.Margin==taxlev[1])), length(unique(dat_abun$Reference.ID[which(dat_abun$Edge.Center.Margin==taxlev[1])])))
effectsdat_abun[8,] <- c("abundance",rownames(mod$b)[2], mod$b[2,1], mod$ci.lb[2], mod$ci.ub[2], length(which(dat_abun$Edge.Center.Margin==taxlev[2])), length(unique(dat_abun$Reference.ID[which(dat_abun$Edge.Center.Margin==taxlev[2])])))
effectsdat_abun[9,] <- c("abundance",rownames(mod$b)[3], mod$b[3,1], mod$ci.lb[3], mod$ci.ub[3], length(which(dat_abun$Edge.Center.Margin==taxlev[3])), length(unique(dat_abun$Reference.ID[which(dat_abun$Edge.Center.Margin==taxlev[3])])))
effectsdat_abun[10,] <- c("abundance",rownames(mod$b)[4], mod$b[4,1], mod$ci.lb[4], mod$ci.ub[4], length(which(dat_abun$Edge.Center.Margin==taxlev[4])), length(unique(dat_abun$Reference.ID[which(dat_abun$Edge.Center.Margin==taxlev[4])])))


## Cropping.System
levels(dat_abun$Cropping.System)
taxlev <- c("Cereal", "Mixed", "Perennial") ##categories represented by at least 3 studies
mod <- rma.mv(hedges, variance, mods= ~0+Cropping.System, 
              random = list(~ 1 | Obs.ID, ~ 1 | StudyNum), tdist=TRUE, data=dat_abun[which(dat_abun$Cropping.System %in% c("Cereal", "Mixed", "Perennial")),])
sink(file = "Model-output4_Single-moderator_abundance_Cropping.System.txt"); summary(mod); sink(file = NULL)
mod
effectsdat_abun[11,] <- c("abundance",rownames(mod$b)[1], mod$b[1,1], mod$ci.lb[1], mod$ci.ub[1], length(which(dat_abun$Cropping.System==taxlev[1])), length(unique(dat_abun$Reference.ID[which(dat_abun$Cropping.System==taxlev[1])])))
effectsdat_abun[12,] <- c("abundance",rownames(mod$b)[2], mod$b[2,1], mod$ci.lb[2], mod$ci.ub[2], length(which(dat_abun$Cropping.System==taxlev[2])), length(unique(dat_abun$Reference.ID[which(dat_abun$Cropping.System==taxlev[2])])))
effectsdat_abun[13,] <- c("abundance",rownames(mod$b)[3], mod$b[3,1], mod$ci.lb[3], mod$ci.ub[3], length(which(dat_abun$Cropping.System==taxlev[3])), length(unique(dat_abun$Reference.ID[which(dat_abun$Cropping.System==taxlev[3])])))


## Sampling Method 
levels(dat_abun$Sampling.Method)
taxlev <- c("Combination", "Pan Traps", "Transects") ##categories represented by at least 3 studies
mod <- rma.mv(hedges, variance, mods= ~0+Sampling.Method, 
              random = list(~ 1 | Obs.ID, ~ 1 | StudyNum), tdist=TRUE, data=dat_abun[which(dat_abun$Sampling.Method %in% c("Combination", "Pan Traps", "Transects")),])
sink(file = "Model-output4_Single-moderator_abundance_Sampling.Method.txt"); summary(mod); sink(file = NULL)
mod
effectsdat_abun[14,] <- c("abundance",rownames(mod$b)[1], mod$b[1,1], mod$ci.lb[1], mod$ci.ub[1], length(which(dat_abun$Sampling.Method==taxlev[1])), length(unique(dat_abun$Reference.ID[which(dat_abun$Sampling.Method==taxlev[1])])))
effectsdat_abun[15,] <- c("abundance",rownames(mod$b)[2], mod$b[2,1], mod$ci.lb[2], mod$ci.ub[2], length(which(dat_abun$Sampling.Method==taxlev[2])), length(unique(dat_abun$Reference.ID[which(dat_abun$Sampling.Method==taxlev[2])])))
effectsdat_abun[16,] <- c("abundance",rownames(mod$b)[3], mod$b[3,1], mod$ci.lb[3], mod$ci.ub[3], length(which(dat_abun$Sampling.Method==taxlev[3])), length(unique(dat_abun$Reference.ID[which(dat_abun$Sampling.Method==taxlev[3])])))


## Landscape Complexity
levels(dat_abun$Landscape.Complexity)
taxlev <- c("Complex", "Gradient", "Simple") 
mod <- rma.mv(hedges, variance, mods= ~0+Landscape.Complexity, 
              random = list(~ 1 | Obs.ID, ~ 1 | StudyNum), tdist=TRUE, data=dat_abun[which(dat_abun$Landscape.Complexity %in% c("Complex", "Gradient", "Simple")),])
sink(file = "Model-output4_Single-moderator_abundance_Landscape.Complexity.txt"); summary(mod); sink(file = NULL)
mod
effectsdat_abun[17,] <- c("abundance",rownames(mod$b)[1], mod$b[1,1], mod$ci.lb[1], mod$ci.ub[1], length(which(dat_abun$Landscape.Complexity==taxlev[1])), length(unique(dat_abun$Reference.ID[which(dat_abun$Landscape.Complexity==taxlev[1])])))
effectsdat_abun[18,] <- c("abundance",rownames(mod$b)[2], mod$b[2,1], mod$ci.lb[2], mod$ci.ub[2], length(which(dat_abun$Landscape.Complexity==taxlev[2])), length(unique(dat_abun$Reference.ID[which(dat_abun$Landscape.Complexity==taxlev[2])])))
effectsdat_abun[19,] <- c("abundance",rownames(mod$b)[3], mod$b[3,1], mod$ci.lb[3], mod$ci.ub[3], length(which(dat_abun$Landscape.Complexity==taxlev[3])), length(unique(dat_abun$Reference.ID[which(dat_abun$Landscape.Complexity==taxlev[3])])))

### now analyze diversity responses 

## Taxonomic.Category
levels(dat_dive$Taxonomic.Category)
taxlev <- c("Bees", "Bumblebees", "Butterflies", "Hoverflies", "Moths", "Solitary Bees") ##categories represented by at least 3 studies, minus "Pollinating Insects" (repetitive/misleading)
mod <- rma.mv(hedges, variance, mods= ~0+Taxonomic.Category, 
              random = list(~ 1 | Obs.ID, ~ 1 | StudyNum), tdist=TRUE, 
              data=dat_dive[which(dat_dive$Taxonomic.Category %in% c("Bees", "Bumblebees", "Butterflies", "Hoverflies", "Moths", "Solitary Bees")),])
sink(file = "Model-output4_Single-moderator-diversity_Taxonomic.Category.txt"); summary(mod); sink(file = NULL)
mod
effectsdat_dive[2,] <- c("diversity",rownames(mod$b)[1], mod$b[1,1], mod$ci.lb[1], mod$ci.ub[1], length(which(dat_dive$Taxonomic.Category==taxlev[1])), length(unique(dat_dive$Reference.ID[which(dat_dive$Taxonomic.Category==taxlev[1])])))
effectsdat_dive[3,] <- c("diversity",rownames(mod$b)[2], mod$b[2,1], mod$ci.lb[2], mod$ci.ub[2], length(which(dat_dive$Taxonomic.Category==taxlev[2])), length(unique(dat_dive$Reference.ID[which(dat_dive$Taxonomic.Category==taxlev[2])])))
effectsdat_dive[4,] <- c("diversity",rownames(mod$b)[3], mod$b[3,1], mod$ci.lb[3], mod$ci.ub[3], length(which(dat_dive$Taxonomic.Category==taxlev[3])), length(unique(dat_dive$Reference.ID[which(dat_dive$Taxonomic.Category==taxlev[3])])))
effectsdat_dive[5,] <- c("diversity",rownames(mod$b)[4], mod$b[4,1], mod$ci.lb[4], mod$ci.ub[4], length(which(dat_dive$Taxonomic.Category==taxlev[4])), length(unique(dat_dive$Reference.ID[which(dat_dive$Taxonomic.Category==taxlev[4])])))
effectsdat_dive[6,] <- c("diversity",rownames(mod$b)[5], mod$b[5,1], mod$ci.lb[5], mod$ci.ub[5], length(which(dat_dive$Taxonomic.Category==taxlev[5])), length(unique(dat_dive$Reference.ID[which(dat_dive$Taxonomic.Category==taxlev[5])])))
effectsdat_dive[7,] <- c("diversity",rownames(mod$b)[6], mod$b[6,1], mod$ci.lb[6], mod$ci.ub[6], length(which(dat_dive$Taxonomic.Category==taxlev[6])), length(unique(dat_dive$Reference.ID[which(dat_dive$Taxonomic.Category==taxlev[6])])))



## Edge.Center.Margin
levels(dat_dive$Edge.Center.Margin)
taxlev <- c("Across System", "Center", "Edge", "Margin")
mod <- rma.mv(hedges, variance, mods= ~0+Edge.Center.Margin, 
              random = list(~ 1 | Obs.ID, ~ 1 | StudyNum), tdist=TRUE, 
              data=dat_dive[which(dat_dive$Edge.Center.Margin %in% c("Across System", "Center", "Edge", "Margin")),])
sink(file = "Model-output4_Single-moderator-diversity_Edge.Center.Margin.txt"); summary(mod); sink(file = NULL)
mod
effectsdat_dive[8,] <- c("diversity",rownames(mod$b)[1], mod$b[1,1], mod$ci.lb[1], mod$ci.ub[1], length(which(dat_dive$Edge.Center.Margin==taxlev[1])), length(unique(dat_dive$Reference.ID[which(dat_dive$Edge.Center.Margin==taxlev[1])])))
effectsdat_dive[9,] <- c("diversity",rownames(mod$b)[2], mod$b[2,1], mod$ci.lb[2], mod$ci.ub[2], length(which(dat_dive$Edge.Center.Margin==taxlev[2])), length(unique(dat_dive$Reference.ID[which(dat_dive$Edge.Center.Margin==taxlev[2])])))
effectsdat_dive[10,] <- c("diversity",rownames(mod$b)[3], mod$b[3,1], mod$ci.lb[3], mod$ci.ub[3], length(which(dat_dive$Edge.Center.Margin==taxlev[3])), length(unique(dat_dive$Reference.ID[which(dat_dive$Edge.Center.Margin==taxlev[3])])))
effectsdat_dive[11,] <- c("diversity",rownames(mod$b)[4], mod$b[4,1], mod$ci.lb[4], mod$ci.ub[4], length(which(dat_dive$Edge.Center.Margin==taxlev[4])), length(unique(dat_dive$Reference.ID[which(dat_dive$Edge.Center.Margin==taxlev[4])])))


## Cropping.System
levels(dat_dive$Cropping.System)
taxlev <- c("Cereal", "Mixed", "Pasture", "Perennial") ##categories represented by at least 3 studies
mod <- rma.mv(hedges, variance, mods= ~0+Cropping.System, 
              random = list(~ 1 | Obs.ID, ~ 1 | StudyNum), tdist=TRUE, 
              data=dat_dive[which(dat_dive$Cropping.System %in% c("Cereal", "Mixed", "Pasture", "Perennial")),])
sink(file = "Model-output4_Single-moderator-diversity_Cropping.System.txt"); summary(mod); sink(file = NULL)
mod
effectsdat_dive[12,] <- c("diversity",rownames(mod$b)[1], mod$b[1,1], mod$ci.lb[1], mod$ci.ub[1], length(which(dat_dive$Cropping.System==taxlev[1])), length(unique(dat_dive$Reference.ID[which(dat_dive$Cropping.System==taxlev[1])])))
effectsdat_dive[13,] <- c("diversity",rownames(mod$b)[2], mod$b[2,1], mod$ci.lb[2], mod$ci.ub[2], length(which(dat_dive$Cropping.System==taxlev[2])), length(unique(dat_dive$Reference.ID[which(dat_dive$Cropping.System==taxlev[2])])))
effectsdat_dive[14,] <- c("diversity",rownames(mod$b)[3], mod$b[3,1], mod$ci.lb[3], mod$ci.ub[3], length(which(dat_dive$Cropping.System==taxlev[3])), length(unique(dat_dive$Reference.ID[which(dat_dive$Cropping.System==taxlev[3])])))
effectsdat_dive[15,] <- c("diversity",rownames(mod$b)[4], mod$b[4,1], mod$ci.lb[4], mod$ci.ub[4], length(which(dat_dive$Cropping.System==taxlev[4])), length(unique(dat_dive$Reference.ID[which(dat_dive$Cropping.System==taxlev[4])])))


## Sampling Method 
levels(dat_dive$Sampling.Method)
taxlev <- c("Combination", "Other", "Pan Traps", "Transects", "Trap Nests") ##categories represented by at least 3 studies
mod <- rma.mv(hedges, variance, mods= ~0+Sampling.Method, 
              random = list(~ 1 | Obs.ID, ~ 1 | StudyNum), tdist=TRUE, data=dat_dive[which(dat_dive$Sampling.Method %in% c("Combination", "Pan Traps", "Transects", "Trap Nests", "Other")),])
sink(file = "Model-output4_Single-moderator_diversity_Sampling.Method.txt"); summary(mod); sink(file = NULL)
mod
effectsdat_dive[16,] <- c("diversity",rownames(mod$b)[1], mod$b[1,1], mod$ci.lb[1], mod$ci.ub[1], length(which(dat_dive$Sampling.Method==taxlev[1])), length(unique(dat_dive$Reference.ID[which(dat_dive$Sampling.Method==taxlev[1])])))
effectsdat_dive[17,] <- c("diversity",rownames(mod$b)[2], mod$b[2,1], mod$ci.lb[2], mod$ci.ub[2], length(which(dat_dive$Sampling.Method==taxlev[2])), length(unique(dat_dive$Reference.ID[which(dat_dive$Sampling.Method==taxlev[2])])))
effectsdat_dive[18,] <- c("diversity",rownames(mod$b)[3], mod$b[3,1], mod$ci.lb[3], mod$ci.ub[3], length(which(dat_dive$Sampling.Method==taxlev[3])), length(unique(dat_dive$Reference.ID[which(dat_dive$Sampling.Method==taxlev[3])])))
effectsdat_dive[19,] <- c("diversity",rownames(mod$b)[4], mod$b[4,1], mod$ci.lb[4], mod$ci.ub[4], length(which(dat_dive$Sampling.Method==taxlev[4])), length(unique(dat_dive$Reference.ID[which(dat_dive$Sampling.Method==taxlev[4])])))
effectsdat_dive[20,] <- c("diversity",rownames(mod$b)[5], mod$b[5,1], mod$ci.lb[5], mod$ci.ub[5], length(which(dat_dive$Sampling.Method==taxlev[5])), length(unique(dat_dive$Reference.ID[which(dat_dive$Sampling.Method==taxlev[5])])))


## Landscape Complexity 
levels(dat_dive$Landscape.Complexity)
taxlev <- c("Complex", "Gradient", "Simple")
mod <- rma.mv(hedges, variance, mods= ~0+Landscape.Complexity, 
              random = list(~ 1 | Obs.ID, ~ 1 | StudyNum), tdist=TRUE, data=dat_dive[which(dat_dive$Landscape.Complexity %in% c("Complex", "Gradient", "Simple")),])
sink(file = "Model-output4_Single-moderator_diversity_Landscape.Complexity.txt"); summary(mod); sink(file = NULL)
mod
effectsdat_dive[21,] <- c("diversity",rownames(mod$b)[1], mod$b[1,1], mod$ci.lb[1], mod$ci.ub[1], length(which(dat_dive$Landscape.Complexity==taxlev[1])), length(unique(dat_dive$Reference.ID[which(dat_dive$Landscape.Complexity==taxlev[1])])))
effectsdat_dive[22,] <- c("diversity",rownames(mod$b)[2], mod$b[2,1], mod$ci.lb[2], mod$ci.ub[2], length(which(dat_dive$Landscape.Complexity==taxlev[2])), length(unique(dat_dive$Reference.ID[which(dat_dive$Landscape.Complexity==taxlev[2])])))
effectsdat_dive[23,] <- c("diversity",rownames(mod$b)[3], mod$b[3,1], mod$ci.lb[3], mod$ci.ub[3], length(which(dat_dive$Landscape.Complexity==taxlev[3])), length(unique(dat_dive$Reference.ID[which(dat_dive$Landscape.Complexity==taxlev[3])])))




write.csv(effectsdat_abun, "Model-output4_Single-moderator_effects-by-category_abundance.csv")
write.csv(effectsdat_dive, "Model-output4_Single-moderator_effects-by-category_diversity.csv")
