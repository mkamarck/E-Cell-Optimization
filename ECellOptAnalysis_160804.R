#Analysis from the next run 160804 (8-4-16)
#libraries
library(ggplot2)
library(reshape2)
library(plyr)
library(rsm)

#I'm not sure how to look at the receptor for this - I might not want to include it as a factor, and then just compare to see whether the 15 and 17.5 give the same results
#start by analyzing without receptor type
#import data
df <- read.csv("/Volumes/mainland/Projects/TAARs/E\ cell\ Optimization/Results/ECellOptimization_160804.csv")
#df <- df[1:96,1:14]
#code the data
dfforrsm <- subset(df, UseforRSM == "y")
df.coded <- coded.data(dfforrsm, x1 ~(Receptor-16.25)/1.25, x2 ~ (SV40 - 2)/1, x3 ~ (CRE - 23.75)/6.25)

#make dataframes for each of the three kinds of data
lucData <- df.coded[,c(1:3,6,9,12)]
RLData <- df.coded[,c(1:3,7,10,13)]
normData <- df.coded[,c(1:3,8,11,14)]

#rsm for the normalized data
#rsm analysis
normData.Low <- subset(normData, x1 == -1)
normData.High <- subset(normData, x1 == 1)
#normData_coded <- coded.data(normData, x2 ~ (SV40 - 2)/1, x3 ~ (CRE - 23.75)/6.25)

norm_30 <- rsm(TMA_30_norm ~ FO(x1,x2,x3) + TWI(x2,x3), data = normData)
summary(norm_30)
persp (norm_30, ~ x1 +x2 + x3, col = rainbow(50), contours = "colors") 
#I can't use the receptor in the SO equation becuas eof the way I designe the experiment.

norm_30 <- rsm(TMA_30_norm ~ SO(x2,x3), data = normData.Low)
summary(norm_30)
persp (norm_30, ~ x2 + x3, col = rainbow(50), contours = "colors")
norm_30 <- rsm(TMA_30_norm ~ SO(x2,x3), data = normData.High)
summary(norm_30)
persp (norm_30, ~ x2 + x3, col = rainbow(50), contours = "colors") 
#This looks about the same with both receptor concentrations

#try with norm_100
#the lack of fit was significant, so I should try second order maybe?
norm_100 <- rsm(TMA_100_norm ~ SO(x2,x3), data = normData.Low)
summary(norm_100)
persp (norm_100, ~ x2 + x3, col = rainbow(50), contours = "colors") 
norm_100 <- rsm(TMA_100_norm ~ SO(x2,x3), data = normData.High)
summary(norm_100)
persp (norm_100, ~ x2 + x3, col = rainbow(50), contours = "colors") 

#Norm_300
norm_300 <- rsm(TMA_300_norm ~ SO(x2,x3), data = normData.Low)
summary(norm_300)
persp (norm_300, ~ x2 + x3, col = rainbow(50), contours = "colors") 
norm_300 <- rsm(TMA_300_norm ~ SO(x2,x3), data = normData.High)
summary(norm_300)
persp (norm_300, ~ x2 + x3, col = rainbow(50), contours = "colors")

###################################
#Lets look at the p values
#start with uncoded values
normData <- dfforrsm[,c(1:3,8,11,14)]
normData.melt <- melt(normData, c("Receptor", "SV40", "CRE"))
normDatattest_30to100 <- ddply(subset(normData.melt, variable %in% c("TMA_30_norm", "TMA_100_norm")), .variables = c("Receptor", "SV40", "CRE"), function(x) test = t.test(value~variable, x)$p.value)
normDatattest_30to100 <- coded.data(normDatattest_30to100, x1 ~(Receptor-16.25)/1.25, x2 ~ (SV40 - 2)/1, x3 ~ (CRE - 23.75)/6.25)

# exampledf <- subset(normData.melt, variable %in% c("TMA_30_norm", "TMA_100_norm") & Receptor == 10 & CRE == 10 & SV40 == 2.5)
# test <- t.test(value~variable, exampledf)$p.value
#visualize the ttests- copy joel
ggplot(normDatattest_30to100, aes(x = x1, y = V1)) +
  geom_point()

ggplot(normDatattest_30to100, aes(x = x2, y = V1)) +
  geom_point()

ggplot(normDatattest_30to100, aes(x = x3, y = V1)) +
  geom_point()
#can't see much from these tests, except maybe there are some outliers 

#now lets run some RSM stuff
#FO test
ttest_30to100 <- rsm(V1 ~ FO(x1,x2,x3), data = normDatattest_30to100)
summary(ttest_30to100)
persp (ttest_30to100, ~ x1 + x2 + x3, col = rainbow(50), contours = "colors")

#Separate by receptor
normDatattest_30to100.Low <- subset(normDatattest_30to100, x1 == -1) 
normDatattest_30to100.High <- subset(normDatattest_30to100, x1 == 1) 

#rsm
ttest_30to100 <- rsm(V1 ~ SO(x2,x3), data = normDatattest_30to100.Low)
summary(ttest_30to100)
persp (ttest_30to100, ~  x2 + x3, col = rainbow(50), contours = "colors")
ttest_30to100 <- rsm(V1 ~ SO(x2,x3), data = normDatattest_30to100.High)
summary(ttest_30to100)
persp (ttest_30to100, ~  x2 + x3, col = rainbow(50), contours = "colors")


normDatattest_100to300 <- ddply(subset(normData.melt, variable %in% c("TMA_100_norm", "TMA_300_norm")), .variables = c("Receptor", "SV40", "CRE"), function(x) test = t.test(value~variable, x)$p.value)
normDatattest_100to300 <- coded.data(normDatattest_100to300, x1 ~(Receptor-16.25)/1.25, x2 ~ (SV40 - 2)/1, x3 ~ (CRE - 23.75)/6.25)
normDatattest_100to300.Low <- subset(normDatattest_100to300, x1 == -1) 
normDatattest_100to300.High <- subset(normDatattest_100to300, x1 == 1) 
#rsm
ttest_100to300 <- rsm(V1 ~ SO(x2,x3), data = normDatattest_100to300.Low)
summary(ttest_100to300)
persp (ttest_100to300, ~  x2 + x3, col = rainbow(50), contours = "colors")
ttest_100to300 <- rsm(V1 ~ SO(x2,x3), data = normDatattest_100to300.High)
summary(ttest_100to300)
persp (ttest_100to300, ~  x2 + x3, col = rainbow(50), contours = "colors")
#this one looks weird, I think the receptor value needs to be lower, in general. the 15 looks much more like we want...

##################
#try an anova
normData.aov <- ddply(subset(normData.melt, variable %in% c("TMA_30_norm", "TMA_100_norm", "TMA_300_norm")), .variables = c("Receptor", "SV40", "CRE"), function(x) test = summary(aov(value~variable, x))[[1]][["Pr(>F)"]][[1]])
#code it.
normData.aov <- coded.data(normData.aov, x1 ~(Receptor-16.25)/1.25, x2 ~ (SV40 - 2)/1, x3 ~ (CRE - 23.75)/6.25)

anova <- rsm(V1 ~ FO(x1,x2,x3), data = normData.aov)
summary(anova)
persp (anova, ~ x1 + x2 + x3, col = rainbow(50), contours = "colors")
#divide it by receptor
normData.aov.Low <- subset(normData.aov, x1 == -1) 
normData.aov.High <- subset(normData.aov, x1 == 1) 
#rsm 2nd order
anova <- rsm(V1 ~ SO(x2,x3), data = normData.aov.Low)
summary(anova)
persp (anova, ~  x2 + x3, col = rainbow(50), contours = "colors")
anova <- rsm(V1 ~ SO(x2,x3), data = normData.aov.High)
summary(anova)
persp (anova, ~  x2 + x3, col = rainbow(50), contours = "colors")
#again high receptor looks totally different from low receptor. 

#Now I need to look to see what the p values are like for the originals
#then eliminate any combos that have p values higher than the originals. 
df.orig <- subset(df, UseforRSM == "n")
df.orig.norm <- df.orig[,c(1:3,8,11,14)]
df.norm.melt <- melt(df.orig.norm, c("Receptor", "SV40", "CRE"))
test.orig <- aov(value~variable, df.norm.melt)
summary(test.orig)
p.orig <- summary(test.orig)[[1]][["Pr(>F)"]][[1]]
p.orig

#eliminate rows with p values less than that with the original summary
normData.low.p <- normData.aov[which(normData.aov$V1<p.orig),]
#I should have been able to predict this, but none of the things have p values lower than the original becuae there are so many more replications of the original p-value
#try taking just the top three (not actually the top one cause that's messed up) from each\
df.orig.norm_short <- df.orig[c(3:5),c(1:3,8,11,14)]
df.norm.melt <- melt(df.orig.norm_short, c("Receptor", "SV40", "CRE"))
test.orig <- aov(value~variable, df.norm.melt)
summary(test.orig)
p.orig <- summary(test.orig)[[1]][["Pr(>F)"]][[1]]
p.orig
normData.low.p <- normData.aov[which(normData.aov$V1<p.orig),]
#hmmm this also has the lowest p value by a lot... I wonder why that is.  Maybe I need to go back to a really low value of receptor with my improved values of SV40 and CRE


                            