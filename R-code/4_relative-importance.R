
# Purpose of this code: 
# Use model selection analysis (R package glmulti) to determine the relative 
# importance of our five moderators (sampling location, cropping system, sampling 
# method, landscape complexity, and taxonomic group) for abundance and diversity 
# effect sizes. This method explores all possible models given a set of moderators, 
# then calculates the relative importance of moderators using the sum of Akaike 
# weights for all models containing each respective moderator. 
# Because of our limited dataset, we weight each model using the the Akaike Information  
# Criterioncorrected for small sample sizes.
# We limit possible models to those that include between one and five main effects but 
# no interaction effects due to small sample size and thus low statistical power to 
# detect such effects. We considered moderators with relative importance values above 
# 0.8 to be important.

# In this model selection approach, 
# the importance value for a particular predictor is equal to the 
# sum of the weights/probabilities for the models in which the variable 
# appears. So, a variable that shows up in lots of models with large 
# weights will receive a high importance value. In that sense, these 
# values can be regarded as the overall support for each variable across 
# all models in the candidate set. In the plot below, vertical dashed line is drawn at 0.8, 
# which is often used as a cutoff to differentiate between important 
# and not so important variables, but this is again a more or less 
# arbitrary division (and a cutoff of 0.5 has also been at times used/suggested).


library(metafor)
library(glmulti)
library(rJava)

dat_abun <- read.csv("dat_abun.csv")
dat_dive <- read.csv("dat_dive.csv")

## for each dataset, 
# response variable: Hedges' g -> "hedges"
# variance: "variance"
# predictor variables: Taxonomic.Category + Cropping.System + Sampling.Method + Sampling.Location + Landscape.Complexity 


# subset our datasets to include only those effect sizes which include all relevant information:
dat_abun1 <- dat_abun[!apply(dat_abun[,c("hedges", "variance", 
                                      "Taxonomic.Category" , 
                                      "Cropping.System" , 
                                      "Edge.Center.Margin" ,
                                      "Sampling.Method" , 
                                      "Landscape.Complexity")], 1, anyNA),]
dat_dive1 <- dat_dive[!apply(dat_dive[,c("hedges", "variance", 
                                         "Taxonomic.Category" , 
                                         "Cropping.System" , 
                                         "Edge.Center.Margin" ,
                                         "Sampling.Method" , 
                                         "Landscape.Complexity")], 1, anyNA),]



### abundance

# model selection
rma.glmulti <- function(formula, data, ...)
  rma(formula, variance, data=data, method="ML", ...)

res_abun <- glmulti(hedges ~ Taxonomic.Category + Cropping.System + Edge.Center.Margin + Sampling.Method + Landscape.Complexity, data=dat_abun1,
               level=1, fitfunction=rma.glmulti, crit="aicc", confsetsize=128)


# plot variable (i.e. terms) importances
ww = exp(-(res_abun@crits - res_abun@crits[1])/2)
ww=ww/sum(ww)
# handle synonymies for interactions
# this translates to unique notations (e.g. x:y and y:x are the same)
clartou=function(res_abun) {
  sort(strsplit(res_abun, ":")[[1]])-> pieces
  if (length(pieces)>1) paste(pieces[1],":",pieces[2],sep="")
  else res_abun
}
# list terms in models
tet = lapply(res_abun@formulas, function(res_abun) sapply(attr(delete.response(terms(res_abun)),"term.labels"), clartou))
# all unique terms
unique(unlist(tet))-> allt_abun
# importances
sapply(allt_abun, function(res_abun) sum(ww[sapply(tet, function(t) res_abun%in%t)]))-> imp_abun
allt_abun2 <- c("Sampling Location", "Landscape Complexity", "Crop Type", 
           "Sampling Method", "Pollinator Group")



### diversity

# model selection
rma.glmulti <- function(formula, data, ...)
  rma(formula, variance, data=data, method="ML", ...)

res_dive <- glmulti(hedges ~ Taxonomic.Category + Cropping.System + Edge.Center.Margin + Sampling.Method + Landscape.Complexity, data=dat_dive1,
                    level=1, fitfunction=rma.glmulti, crit="aicc", confsetsize=128)


# plot variable (i.e. terms) importances
ww = exp(-(res_dive@crits - res_dive@crits[1])/2)
ww=ww/sum(ww)
# handle synonymies for interactions
# this translates to unique notations (e.g. x:y and y:x are the same)
clartou=function(res_dive) {
  sort(strsplit(res_dive, ":")[[1]])-> pieces
  if (length(pieces)>1) paste(pieces[1],":",pieces[2],sep="")
  else res_dive
}
# list terms in models
tet = lapply(res_dive@formulas, function(res_dive) sapply(attr(delete.response(terms(res_dive)),"term.labels"), clartou))
# all unique terms
unique(unlist(tet))-> allt_dive
# importances
sapply(allt_dive, function(res_dive) sum(ww[sapply(tet, function(t) res_dive%in%t)]))-> imp_dive
allt_dive2 <- c("Pollinator Group", "Crop Type", "Landscape Complexity",  "Sampling Method",
                "Sampling Location")


# create a dataframe which includes relative importances for abundance and diversity
abundat <- data.frame(name=allt_abun2, 
                     val.abun=imp_abun)
divedat <- data.frame(name=allt_dive2, 
                      val.dive=imp_dive)

mdat <- merge(abundat, divedat, by="name", all=TRUE)

mdat <- mdat[order(mdat$val.abun),]

# Plot

jpeg("Figures-4_Variable-importance.jpeg", width = 4000, height = 2400, res=600)

layout(matrix(c(1, 2,  # First, second
                     3, 3), # and third plot
                   nrow = 2,
                   ncol = 2,
                   byrow = TRUE), widths=c(1,1), heights=c(0.9,0.05))
par(oma=c(0,12,0,1), mar=c(3,1,3,1))

# get colors 
cols <- natparks.pals("SmokyMtns", 14)

barplot(mdat$val.abun,
        xlab="",xlim=c(0,1), 
        ylab="",horiz=T,las=2, 
        names.arg=mdat$name,
        main="A) Abundance", font.main=1, las=1, mgp=c(1.75,0.5,0), tck=-0.05, xaxp=c(0,1,5),
        col=cols[6])
box(lty = 'solid', col = 'black')
abline(v=0.8, col="black", lty=2, cex=2)

barplot(mdat$val.dive,
        xlab="",xlim=c(0,1), 
        ylab="",horiz=T,las=2, 
        names.arg=NA,
        main="B) Diversity", font.main=1, las=1, mgp=c(1.75,0.5,0), tck=-0.05, xaxp=c(0,1,5),
        col=cols[2])
box(lty = 'solid', col = 'black')
abline(v=0.8, col="black", lty=2, cex=2)

par(mar=c(0,0,0,0))
plot(2,2,xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
text(0.5,0.5, "Relative importance", cex=1)

dev.off()

# compare to fig to double check
imp_abun
imp_dive



# interpretation: 
# Model selection analysis to determine relative importance of moderators show that 
# sampling location is the most important moderators for the response of abundance to 
# organic agriculture. However, taxonomic group is the most important moderator for 
# diversity response to organic agriculture. 



