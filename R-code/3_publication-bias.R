

library(metafor) # version 3.0-2




dat_abun <- read.csv("dat_abun.csv")
dat_dive <- read.csv("dat_dive.csv")


### Eggers test
# This test is to determine if there is significant publication bias. 
# An intercept which is significantly different from zero would indicate the overall 
# relationship between the precision and magnitude of effect sizes is asymmetrical, 
# and therefore biased. 

#  Abundance model
dat <- dat_abun
mod <- rma.mv(hedges, variance, mods = ~ variance, random = list(~ 1 | Obs.ID, ~ 1 | StudyNum), tdist=TRUE, data=dat)
sink(file = "Model-output6_Publication-bias-eggers_abundance.txt"); summary(mod); sink(file = NULL)

mod
# interpretation of the abundance eggers test: 
# the intercept is -0.4219 with confidence interval of -0.8951  to 0.0514
# because the intercept is not different from zero, we can conclude that 
# abundance effect sizes in our dataset are unbiased.

#  Diversity model
dat <- dat_dive
mod <- rma.mv(hedges, variance, mods = ~ variance, random = list(~ 1 | Obs.ID, ~ 1 | StudyNum), tdist=TRUE, data=dat)
sink(file = "Model-output6_Publication-bias-eggers_diversity.txt"); summary(mod); sink(file = NULL)

mod
# interpretation of the diversity eggers test:
# the intercept is 0.1496 with confidence interval of -0.2230 to 0.5222
# because the intercept is not different from zero, we can conclude that
# diversity effect sizes in our dataset are unbiased.


### Sensitivity analysis by jacknife procedure, 
### where each study is removed and overall mean recalculated
# If the mean hedges g does not change significantly after removing 
# each study in the dataset, then we can conclude that the overall
# response is robust against removal of any one study.


dat <- dat_abun
r <- length(unique(dat$Reference.ID))
ids <- unique(dat$Reference.ID)
effectsdat_abun <- data.frame(study=rep(NA,r),
                             mean=rep(NA,r),
                             lci=rep(NA,r),
                             uci=rep(NA,r),
                             nobs=rep(NA,r))
for(i in 1:r){
  study <- ids[i]
  effectsdat_abun$study[i] <- study
  mod <- rma.mv(hedges, variance,random = list(~ 1 | Obs.ID, ~ 1 | StudyNum), tdist=TRUE, data=dat[-which(dat$Reference.ID==study),])
  effectsdat_abun$mean[i] <- mod$b
  effectsdat_abun$lci[i] <- mod$ci.lb
  effectsdat_abun$uci[i] <- mod$ci.ub
  effectsdat_abun$nobs[i] <- mod$k
}

write.csv(effectsdat_abun, "Model-output6_Publication-bias-sensitivity-jacknife_abundance.csv")

effectsdat_abun

# interpretation of the abundance jacknife test: 
# inspect the 'mean' column of effectsdat_abun
# the mean hedges g does not change significantly after removal of each study
# compared to the overall hedges g in "Model-output1_Overall_abundance.txt"
# therefore, the overall abundance response is robust against removal of any one study



dat <- dat_dive
r <- length(unique(dat$Reference.ID))
ids <- unique(dat$Reference.ID)
effectsdat_dive <- data.frame(study=rep(NA,r),
                             mean=rep(NA,r),
                             lci=rep(NA,r),
                             uci=rep(NA,r),
                             nobs=rep(NA,r))
for(i in 1:r){
  study <- ids[i]
  effectsdat_dive$study[i] <- study
  mod <- rma.mv(hedges, variance,random = list(~ 1 | Obs.ID, ~ 1 | StudyNum), tdist=TRUE, data=dat[-which(dat$Reference.ID==study),])
  effectsdat_dive$mean[i] <- mod$b
  effectsdat_dive$lci[i] <- mod$ci.lb
  effectsdat_dive$uci[i] <- mod$ci.ub
  effectsdat_dive$nobs[i] <- mod$k
}

write.csv(effectsdat_dive, "Model-output6_Publication-bias-sensitivity-jacknife_diversity.csv")


effectsdat_dive

# interpretation of the diversity jacknife test: 
# inspect the 'mean' column of effectsdat_dive
# the mean hedges g does not change significantly after removal of each study
# compared to the overall hedges g in "Model-output1_Overall_diversity.txt"
# therefore, the overall diversity response is robust against removal of any one study


mdat <- merge(effectsdat_abun, effectsdat_dive, by="study", all=TRUE, sort = TRUE)


write.csv(mdat, "Model-output6_Publication-bias-sensitivity-jacknife_merged.csv")




