#E cell optimization for 160830 (8-30-16)
#libraries
library(ggplot2)
library(reshape2)
library(plyr)
library(rsm)

#import data
df <- read.csv("/Volumes/mainland-1/Projects/TAARs/E\ cell\ Optimization/Results/E Cell Optimization_160830.csv")
#df  <- read.csv("/Users/mkamarck/Documents/School\ and\ Lab/Mainland\ Lab/E-Cell-Optimization/Results/E Cell Optimization_160830.csv")
#df <- df[1:96,1:14]
#code the data
dfforrsm <- subset(df, UseforRSM == "y")
#dfforrsm <- dfforrsm[,-c(4,5)]
df.coded <- coded.data(dfforrsm, x1 ~(Receptor-10)/5, x2 ~ (SV40 - 1)/.5, x3 ~ (CRE - 25)/10)
#normalized data
normData <- df.coded[,c(1:3,8,11,14)]
#######################################################
#look at some rsm data for normalzied values
norm_30 <- rsm(TMA_30_norm ~ FO(x1,x2,x3) , data = normData)
summary(norm_30)
persp (norm_30, ~ x1 +x2 + x3, col = rainbow(50), contours = "colors") 
#SO
norm_30 <- rsm(TMA_30_norm ~ SO(x1,x2,x3) , data = normData)
summary(norm_30)
persp (norm_30, ~ x1 +x2 + x3, col = rainbow(50), contours = "colors") 

norm_100 <- rsm(TMA_100_norm ~ FO(x1,x2,x3), data = normData)
summary(norm_100)
persp (norm_100, ~ x1+ x2 + x3, col = rainbow(50), contours = "colors") 
norm_100 <- rsm(TMA_100_norm ~ SO(x1,x2,x3), data = normData)
summary(norm_100)
persp (norm_100, ~ x1+ x2 + x3, col = rainbow(50), contours = "colors") 
#300
norm_300 <- rsm(TMA_300_norm ~ FO(x1,x2,x3), data = normData)
summary(norm_300)
persp (norm_300, ~ x1+ x2 + x3, col = rainbow(50), contours = "colors") 
norm_300 <- rsm(TMA_300_norm ~ SO(x1,x2,x3), data = normData)
summary(norm_300)
persp (norm_300, ~ x1+ x2 + x3, col = rainbow(50), contours = "colors") 

#################################
#rsm for anova
normData <- dfforrsm[,c(1:3,8,11,14)]
normData.melt <- melt(normData, c("Receptor", "SV40", "CRE"))

normData.aov <- ddply(subset(normData.melt, variable %in% c("TMA_30_norm", "TMA_100_norm", "TMA_300_norm")), .variables = c("Receptor", "SV40", "CRE"), function(x) test = summary(aov(value~variable, x))[[1]][["Pr(>F)"]][[1]])
#code it.
normData.aov <- coded.data(normData.aov, x1 ~(Receptor-10)/5, x2 ~ (SV40 - 1)/.5, x3 ~ (CRE - 25)/10)

anova <- rsm(V1 ~ FO(x1,x2,x3), data = normData.aov)
summary(anova)
persp (anova, ~ x1 + x2 + x3, col = rainbow(50), contours = "colors")

#rsm 2nd order
anova <- rsm(V1 ~ SO(x1,x2,x3), data = normData.aov)
summary(anova)
persp (anova, ~  x1+x2 + x3, col = rainbow(50), contours = "colors")

#These all look terrible.
#########################################################################################################################
#Take the top 10 p values and graph these to see what the difference is 
#topTen <- normData.aov[which(normData.aov$V1<0.001),]
normData.aov <- ddply(subset(normData.melt, variable %in% c("TMA_30_norm", "TMA_100_norm", "TMA_300_norm")), .variables = c("Receptor", "SV40", "CRE"), function(x) test = summary(aov(value~variable, x))[[1]][["Pr(>F)"]][[1]])

#put the p values with the original 
normData.aov <- decode.data(normData.aov)
normData.aov.all <- merge(normData, normData.aov)
#get the top ten
normData.sort <- normData.aov.all[order(normData.aov.all$V1),]
normData.topTen <- normData.sort[c(1:30),]

#Graph all of these
#how can I even do that.
#try one
normData.topTen.melt <- melt(normData.topTen, c("Receptor", "SV40", "CRE", "V1"))
ggplot(subset(normData.topTen.melt, Receptor == 3.75 & SV40 == 0.75 & CRE == 30),aes(x = variable, y = value)) +
  geom_boxplot()
#plot all of them
ggplot(normData.topTen.melt, aes(x = variable, y = value)) +
  geom_boxplot() +
  facet_grid(V1 ~Receptor)
#cool so these are all going in the same way. 


#now I want to see if the norm values and luc values that are the highest are the same ones that have the lowest p values, or if there is someoverlap
#next I want to check the signal to noise and see if those with the lowest signal to noise are also those with highest p values. 
#calculate signal to noise
normData.signalToNoise <- ddply(normData, .variables = c("Receptor", "SV40", "CRE"), .fun = summarize, mean.TMA_30 = mean(TMA_30_norm), mean.TMA_100 = mean(TMA_100_norm), mean.TMA_300 = mean(TMA_300_norm), sd.TMA_30 = sd(TMA_30_norm), sd.TMA_100 = sd(TMA_100_norm), sd.TMA_300 = sd(TMA_300_norm))
normData.signalToNoise$SN_30 <- normData.signalToNoise$mean.TMA_30/normData.signalToNoise$sd.TMA_30
normData.signalToNoise$SN_100 <- normData.signalToNoise$mean.TMA_100/normData.signalToNoise$sd.TMA_100
normData.signalToNoise$SN_300 <- normData.signalToNoise$mean.TMA_300/normData.signalToNoise$sd.TMA_300
normData.signalToNoise$SN_avg <- (normData.signalToNoise$SN_30 + normData.signalToNoise$SN_100 + normData.signalToNoise$SN_300)/3

normData.signalToNoise_sort <- normData.signalToNoise[order(normData.signalToNoise$SN_avg, decreasing = TRUE),]
normData.signalToNoise_topTen <- normData.signalToNoise_sort[c(1:10),c(1:3)]

#get top ten with the average largest signal to noise ratio over the three TMA concentrations
#mark each top ten with a note for which it is then compile them
#define the dupsBetweenGroups function
dupsBetweenGroups <- function (df, idcol) {
  # df: the data frame
  # idcol: the column which identifies the group each row belongs to
  # Get the data columns to use for finding matches
  datacols <- setdiff(names(df), idcol)
  # Sort by idcol, then datacols. Save order so we can undo the sorting later.
  sortorder <- do.call(order, df)
  df <- df[sortorder,]
  # Find duplicates within each id group (first copy not marked)
  dupWithin <- duplicated(df)
  # With duplicates within each group filtered out, find duplicates between groups. 
  # Need to scan up and down with duplicated() because first copy is not marked.
  dupBetween = rep(NA, nrow(df))
  dupBetween[!dupWithin] <- duplicated(df[!dupWithin,datacols])
  dupBetween[!dupWithin] <- duplicated(df[!dupWithin,datacols], fromLast=TRUE) | dupBetween[!dupWithin]
  # ============= Replace NA's with previous non-NA value ==============
  # This is why we sorted earlier - it was necessary to do this part efficiently
  # Get indexes of non-NA's
  goodIdx <- !is.na(dupBetween)
  # These are the non-NA values from x only
  # Add a leading NA for later use when we index into this vector
  goodVals <- c(NA, dupBetween[goodIdx])
  # Fill the indices of the output vector with the indices pulled from
  # these offsets of goodVals. Add 1 to avoid indexing to zero.
  fillIdx <- cumsum(goodIdx)+1
  # The original vector, now with gaps filled
  dupBetween <- goodVals[fillIdx]
  # Undo the original sort
  dupBetween[sortorder] <- dupBetween
  # Return the vector of which entries are duplicated across groups
  return(dupBetween)
}
#############################

library(DataCombine)
normData.signalToNoise_topTen$df <- "signalToNoise"
normData.pvalue_topTen <- normData.topTen[, c(1:3)]
normData.pvalue_topTen$df <- "pValue"

compare <- rbind(normData.pvalue_topTen, normData.signalToNoise_topTen)
# Find the rows which have duplicates in a different group.
compare$dupRows <- dupsBetweenGroups(compare, "df")
trueForBoth <- compare[which(compare$dupRows == "TRUE" & compare$df == "signalToNoise"),]
#there doesn't seem to be a pattern to these. 
###########
#Lets see what luc data looks like, find the highest value and then find the noise for this
lucData <- df.coded[,c(1:3,6,9,12)]
lucData.sort_30 <- lucData[order(lucData$TMA_30_luc, decreasing = TRUE),]
lucData.sort_30_top20 <- decode.data(lucData.sort_30[c(1:20),c(1:3)])
lucData.sort_30_top20$df <- "TMA_30"
lucData.sort_100 <- lucData[order(lucData$TMA_100_luc, decreasing = TRUE),]
lucData.sort_100_top20 <- decode.data(lucData.sort_100[c(1:20),c(1:3)])
lucData.sort_100_top20$df <- "TMA_100"
lucData.sort_300 <- lucData[order(lucData$TMA_300_luc, decreasing = TRUE),]
lucData.sort_300_top20 <- decode.data(lucData.sort_300[c(1:20),c(1:3)])
lucData.sort_300_top20$df <- "TMA_300"
#combine
luc.combo <- rbind(lucData.sort_30_top20, lucData.sort_100_top20 ,lucData.sort_300_top20)
luc.combo$dupRows <- dupsBetweenGroups(luc.combo, "df")
luc.combo.same <- luc.combo[which(luc.combo$dupRows == "TRUE"),]
#only 4 that weren't duplicated, which means in general that for one value of TMA, all three combos are high. 
#lets now look at overlap between the highest luciferase value and the other top CRE
#get top 10 from TMA at 100
lucData.top10_TMA100 <- decode.data(lucData[order(lucData$TMA_100_luc, decreasing = TRUE),])
lucData.top10_TMA100 <- lucData.top10_TMA100[c(1:10), c(1:3)]
lucData.top10_TMA100$df <- "lucDataTMA100"
#go back and also look at noise for this...

pvalue.topten <- unique(normData.pvalue_topTen)
compare_norm.p.luc <- rbind(pvalue.topten, normData.signalToNoise_topTen, lucData.top10_TMA100)
compare_norm.p.luc$dupRows <- dupsBetweenGroups(compare_norm.p.luc, "df")
trueForBoth <- compare_norm.p.luc[which(compare_norm.p.luc$dupRows == "TRUE"),]

matches.order <- trueForBoth[order(trueForBoth$Receptor, trueForBoth$SV40, trueForBoth$CRE),]


#add in noise to my weird thing
luc.signalToNoise <- ddply(decode.data(lucData), .variables = c("Receptor", "SV40", "CRE"), .fun = summarize, mean.TMA_30 = mean(TMA_30_luc), mean.TMA_100 = mean(TMA_100_luc), mean.TMA_300 = mean(TMA_300_luc), sd.TMA_30 = sd(TMA_30_luc), sd.TMA_100 = sd(TMA_100_luc), sd.TMA_300 = sd(TMA_300_luc))
luc.signalToNoise$SN_30 <- luc.signalToNoise$mean.TMA_30/luc.signalToNoise$sd.TMA_30
luc.signalToNoise$SN_100 <- luc.signalToNoise$mean.TMA_100/luc.signalToNoise$sd.TMA_100
luc.signalToNoise$SN_300 <- luc.signalToNoise$mean.TMA_300/luc.signalToNoise$sd.TMA_300
luc.signalToNoise$SN_avg <- (luc.signalToNoise$SN_30 + luc.signalToNoise$SN_100 + luc.signalToNoise$SN_300)/3

luc.STN.sort <- luc.signalToNoise[order(luc.signalToNoise$SN_avg, decreasing = TRUE),]
luc.STN.topTen <- luc.STN.sort[c(1:10), c(1:3)]
luc.STN.topTen$df <- "lucSignalToNoise"
#combine again
compare_norm.p.luc <- rbind(pvalue.topten, normData.signalToNoise_topTen, lucData.top10_TMA100, luc.STN.topTen)
compare_norm.p.luc$dupRows <- dupsBetweenGroups(compare_norm.p.luc, "df")
trueForBoth <- compare_norm.p.luc[which(compare_norm.p.luc$dupRows == "TRUE"),]
matches.order <- trueForBoth[order(trueForBoth$Receptor, trueForBoth$SV40, trueForBoth$CRE),]
#write.csv(matches.order, "/Volumes/mainland-1/Projects/TAARs/E\ cell\ Optimization/Results/ECellOptimization_160830_topHits.csv")
# receptor 5, SV40 .75, and CRE 30 seems to be the best on all of these fronts. 
#so definitely lower receptor 

#now try modeling norm signal to noise using rsm
normData.signalToNoise.code <- coded.data(normData.signalToNoise, x1 ~(Receptor-10)/5, x2 ~ (SV40 - 1)/.5, x3 ~ (CRE - 25)/10)
rsm.norm <- rsm(SN_avg ~ FO(x1,x2,x3), data = normData.signalToNoise.code)
summary(rsm.norm)
persp (rsm.norm, ~  x1+x2 + x3, col = rainbow(50), contours = "colors")

rsm.norm <- rsm(SN_avg ~ SO(x1,x2,x3), data = normData.signalToNoise.code)
summary(rsm.norm)
#pdf("/Volumes/mainland/Projects/TAARs/E\ cell\ Optimization/Results/ECellAnalysis_RSM_SigToNoise_160811.pdf")
persp (rsm.norm, ~  x1+x2 + x3, col = rainbow(50), contours = "colors")
#dev.off()
#so this confirms what our other analyses are saying which is that low receptor is good along with low SV40 and high CRE
#Do other signal to noise calculations look similar?

######################################################
#new analysis
#compare differences between the new concentrations and the old ones
#1. Receptor = 3.75, SV40 = 0.75, CRE = 30
#2. Receptor = 5, SV40 = 0.5, CRE = 25
#3. Receptor = 3.75, SV40 = 0.5, CRE = 25 (2 and 3 are similar) but 3 is much more noisy at 300uM
receptor  <- 3.75
sv40  <- .5
cre  <- 25
df_compareOrig  <- subset(df, Receptor_Type == "hTAAR5")
df_compareOrig  <- df_compareOrig[,c(1:3,8,11,14)]
#take out some values that came up in the matches list and compare them to the original values
df_compareOrig.1  <- subset(df_compareOrig, (Receptor == 5 & SV40 == 5 & CRE == 10) | (Receptor == receptor & SV40 == sv40 & CRE == cre))
df_compareOrig.melt.1  <- melt(df_compareOrig.1, c("Receptor", "SV40", "CRE"))
ggplot(df_compareOrig.melt.1, aes(x = factor(variable), y = value, color = factor(CRE))) +
    geom_boxplot()
#well, the graph definitely looks like the new one is better, but we should still compare signal to noise between the two

#I need to come back and plot more of these

##################################
#Examine the wells with the experiments with Rho...
#want to plot these side by side with the corresponding hTAAR5 experiments. 
rhoComp  <- subset(df, Receptor == 5)
rhoComp.norm  <- rhoComp[,c(1:3,8,11,14, 16)]
rhoComp.normMelt  <- melt(rhoComp.norm, c("Receptor", "SV40", "CRE", "Receptor_Type"))
ggplot(rhoComp.normMelt, aes(x = variable, y = value, colour = Receptor_Type)) +
  geom_boxplot() +
  facet_grid(SV40~CRE)
##################################
#regular ploting of all the other ones by receptor type
plotReceptor2.5  <- subset(df, Receptor == 2.5 & Receptor_Type == "hTAAR5")
plotReceptor2.5  <- plotReceptor2.5[,c(1:3,8,11,14)]
plotReceptor2.5_melt  <- melt(plotReceptor2.5, c("Receptor", "SV40", "CRE"))
ggplot(plotReceptor2.5_melt, aes(x = variable, y = value, colour = factor(SV40))) +
  geom_boxplot() +
  facet_grid(SV40~CRE)
#receptor = 3.75
plotReceptor3.75  <- subset(df, Receptor == 3.75 & Receptor_Type == "hTAAR5")
plotReceptor3.75  <- plotReceptor3.75[,c(1:3,8,11,14)]
plotReceptor3.75_melt  <- melt(plotReceptor3.75, c("Receptor", "SV40", "CRE"))
ggplot(plotReceptor3.75_melt, aes(x = variable, y = value, colour = factor(SV40))) +
  geom_boxplot() +
  facet_grid(SV40~CRE)
################################
#ploting of all of these by SV40
#I think either 0.5 or 0.75 will have the best outcomes
#receptor = 3.75
plotSV40  <- subset(dfforrsm, SV40 == .5)
plotSV40  <- plotSV40[,c(1:3,8,11,14)]
plotSV40.melt  <- melt(plotSV40, c("Receptor", "SV40", "CRE"))
ggplot(plotSV40.melt, aes(x = variable, y = value, colour = factor(Receptor))) +
  geom_boxplot() +
  facet_grid(Receptor~CRE)
#for SV40 at 0.5 I like receptor at 3.75 or at 5, 2.5 doesn't give a big enough response
#for SV40 at .75 it looks really noisy in general - receptor at 3.75 and CRE at 30 was on the top list
#if .75 looks weird across the board it could be that I didn't mix the beginning component well enough - just a thought. 
#1 is also really noisy
################################
#Trying to figure out what other combinations I will try for the dose response runs.
#top values for TMA 300 normalized
normData.sort2  <- normData[order(normData$TMA_30_norm, decreasing = FALSE),]
#I ended up using most of the combnations from the combination database
################################
library(drc)
#plot the data from the h3A cells. 
#graph the h3a results.
h3A <- read.csv("/Users/mkamarck/Documents/School\ and\ Lab/Mainland\ Lab/E-Cell-Optimization/Results/ECellOptimization_H3AControl_160830.csv")
h3A_avgs <- ddply(h3A, .variables = c("TMA.Concentration", "Plate_number", "Receptor"), .fun = summarize, average = mean(norm), standdev = sd(norm))
h3A_avgs$logTMA <- log10(h3A_avgs$TMA.Concentration)
for(i in 1:length(h3A_avgs$TMA.Concentration)){
  if(h3A_avgs$TMA.Concentration[i] == 0){
    h3A_avgs$logTMA[i] = -10
  }
}
#graph both plates together without the fit line
ggplot(h3A_avgs, aes(x = logTMA, y = average, colour = factor(Receptor))) +
  geom_point() + 
  geom_errorbar(aes(ymin = average - (standdev/sqrt(3)), ymax = average + (standdev/sqrt(3)))) +
  facet_grid(.~Plate_number)

#graph them with the fit line
h3A.sub <- subset(h3A, Plate_number == 7 & Receptor == 830)
h3A.sub$logTMA <- log10(h3A.sub$TMA.Concentration)
for(i in 1:length(h3A.sub$TMA.Concentration)){
  if(h3A.sub$TMA.Concentration[i] == 0){
    h3A.sub$logTMA[i] = -10
  }
}
h3A.sub <- na.omit(h3A.sub)
h3A.sub[which(h3A.sub$logTMA==Inf)] = NA

#the key is that we need to use the non-log form of TMA becuse the equation turns it into the log form on its own
MOD1 <- drm(formula = norm ~ TMA.Concentration, data = h3A.sub[complete.cases(h3A.sub),], lowerl=c(-10,NA,0.0000000000001),upperl=c(NA,10,0.01),fct = LL.4(fixed=c(1,NA,NA,NA),names=(c("Slope", "Top", "Bottom", "ED"))))
summary(MOD1)
#and do the plate 7 for rho
h3A.sub <- subset(h3A, Plate_number == 7 & Receptor == 999)
h3A.sub$logTMA <- log10(h3A.sub$TMA.Concentration)
h3A.sub[which(h3A.sub$logTMA==Inf)] = 10
MOD2 <- drm(formula = norm ~ TMA.Concentration, data = h3A.sub[complete.cases(h3A.sub),], fct = LL.4(fixed=c(1,NA,NA,NA),names=(c("Slope", "Top", "Bottom", "ED"))))
summary(MOD2)

# MOD1<-drm(formula = NormalizedLuc ~ concentration, data = x[complete.cases(x),], fct = LL.4(fixed=c(1,NA,NA,NA),names=(c("Slope", "Top", "Bottom", "ED"))))
# plotDR(temp,plot=TRUE)
eq1 = function(x){summary(MOD1)$coefficients[2,1] + (summary(MOD1)$coefficients[1,1]-summary(MOD1)$coefficients[2,1])/(1+10^((log10(summary(MOD1)$coefficients[3,1]) - x)*1))}
eq2 = function(x){summary(MOD2)$coefficients[2,1] + (summary(MOD2)$coefficients[1,1]-summary(MOD2)$coefficients[2,1])/(1+10^((log10(summary(MOD2)$coefficients[3,1]) - x)*1))}

#the first problem was that i wasn't using the log10 for TMA concentrations so the model didn't fit what I was actually graphing
#pdf("/Users/mkamarck/Documents/School\ and\ Lab/Mainland\ Lab/E-Cell-Optimization/Results/160830_h3ADR_plate7.pdf")
ggplot(subset(h3A_avgs, Plate_number == 7), aes(x = logTMA, y = average, colour = factor(Receptor))) +
  geom_point() +
  geom_errorbar(aes(ymin = average - (standdev/sqrt(6)), ymax = average + (standdev/sqrt(6)))) +
  ylab ("Normalized Luciferase Value") +
  stat_function(fun = eq1, color="#F8766D") +
  stat_function(fun = eq2) 
#dev.off()


#I will eventually figure out how to calculate p values for these things.   


