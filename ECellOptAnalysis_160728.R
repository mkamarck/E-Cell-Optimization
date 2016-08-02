#Analysis from the next run 160728 (7-28-16)
#libraries
library(ggplot2)
library(reshape2)
library(plyr)
library(rsm)

#import data
df <- read.csv("/Volumes/mainland/Projects/TAARs/E\ cell\ Optimization/Results/ECellOptimization_160728.csv")
df <- df[1:96,1:14]
#make dataframes for each of the three kinds of data
lucData <- df[,c(1:6,9,12)]
RLData <- df[,c(1:5,7,10,13)]
normData <- df[,c(1:5,8,11,14)]

#######################################################################
#First lets run an RSM analysis on the normalized data to see how it looks
#code the normalized data
#coded.data(ChemReact1, x1 ~ (Time - 85)/5, x2 ~ (Temp - 175)/5)
#this coding will be different from the last analysis
#receptor - middle is 17.5 and difference from middle is 7.5
normData_coded <- coded.data(normData, x1 ~(Receptor-17.5)/7.5, x2 ~ (SV40 - 8.75)/6.25, x3 ~ (CRE - 15)/5)
#rsm analysis
norm_30 <- rsm(TMA_30_norm ~ FO(x1,x2,x3), data = normData_coded)
summary(norm_30)
contour(norm_30, ~x1+x2+x3, image = TRUE)
#the lack of fit was significant, so I should try second order maybe?
norm_30 <- rsm(TMA_30_norm ~ SO(x1,x2,x3), data = normData_coded)
summary(norm_30)
contour(norm_30, ~x1+x2+x3, image = TRUE)
#lack of fit was still significant, try other combos
norm_30 <- rsm(TMA_30_norm ~ FO(x1,x2,x3) + TWI(x1,x2,x3), data = normData_coded)
summary(norm_30)
contour(norm_30, ~x1+x2+x3, image = TRUE)

norm_30 <- rsm(TMA_30_norm ~ SO(x1,x2,x3) + TWI(x1,x2,x3), data = normData_coded)
summary(norm_30)
contour(norm_30, ~x1+x2+x3, image = TRUE)
persp (norm_30, ~ x1 + x2 + x3, col = rainbow(50), contours = "colors") #These three graphs are actually quite useful.

norm_30 <- rsm(TMA_30_norm ~ SO(x2,x3) + TWI(x2,x3) + PQ(x2,x3), data = normData_coded)
summary(norm_30)
contour(norm_30, ~x2+x3, image = TRUE)

#I'm clearly not doing the analysis right becuase to my eyes it looks like the top values of the norm are at high cre and low SV40, but here they are showing the s
#stationary points at the opposite
#this means that I'm either not coding the variables correctly
#or I'm not running the analysis correclty
#or I don't understand the output summary correctly

#lets try the p value thing now that we have replications of everything
# delete.na <- function(DF, n=0) {
#   log <- apply(df, 2, is.na)
#   logindex <- apply(log, 1, function(x) sum(x) <= n)
#   df[logindex, ]
# }
#######################################################################
normData.melt <- melt(normData, c("Row", "Column", "Receptor", "SV40", "CRE"))
normDatattest_30to100 <- ddply(subset(normData.melt, variable %in% c("TMA_30_norm", "TMA_100_norm")), .variables = c("Receptor", "SV40", "CRE"), function(x) test = t.test(value~variable, x)$p.value)

# exampledf <- subset(normData.melt, variable %in% c("TMA_30_norm", "TMA_100_norm") & Receptor == 10 & CRE == 10 & SV40 == 2.5)
# test <- t.test(value~variable, exampledf)$p.value
#visualize the ttests- copy joel
ggplot(normDatattest_30to100, aes(x = Receptor, y = V1)) +
  geom_point()

ggplot(normDatattest_30to100, aes(x = SV40, y = V1)) +
  geom_point()

ggplot(normDatattest_30to100, aes(x = CRE, y = V1)) +
  geom_point()
#nothing obvious from these plots except that maybe CRE works better at higher values


#rsm analysis of the ttests
#code the ttest data
ttest_30to100.coded<- coded.data(normDatattest_30to100, x1 ~(Receptor-17.5)/7.5, x2 ~ (SV40 - 8.75)/6.25, x3 ~ (CRE - 15)/5)

#FO test
ttest_30to100 <- rsm(V1 ~ FO(x1,x2,x3), data = ttest_30to100.coded)
summary(ttest_30to100)
persp (ttest_30to100, ~ x1 + x2 + x3, col = rainbow(50), contours = "colors")

ttest_30to100 <- rsm(V1 ~ SO(x1,x2,x3), data = ttest_30to100.coded)
summary(ttest_30to100)
contour(ttest_30to100, ~x1+x2+x3, image = TRUE)
persp (ttest_30to100, ~ x1 + x2 + x3, col = rainbow(50), contours = "colors")
#this definitely isn't fitting right

ttest_30to100 <- rsm(V1 ~ SO(x1,x2,x3) + TWI(x1,x2,x3), data = ttest_30to100.coded)
summary(ttest_30to100)
persp (ttest_30to100, ~ x1 + x2 + x3, col = rainbow(50), contours = "colors")

ttest_30to100 <- rsm(V1 ~ FO(x1,x2,x3) + TWI(x1,x2,x3), data = ttest_30to100.coded)
summary(ttest_30to100)
persp (ttest_30to100, ~ x1 + x2 + x3, col = rainbow(50), contours = "colors")

#It really looks like receptor is not having a lot of effect compared to SV40 and CRE balance

#######################################################################
#Lets try the ttest between 100 and 300 - even though the repeats of 300 are questionable so this data may not be helpful.
normDatattest_100to300 <- ddply(subset(normData.melt, variable %in% c("TMA_300_norm", "TMA_100_norm")), .variables = c("Receptor", "SV40", "CRE"), function(x) test = t.test(value~variable, x)$p.value)
#visualize
ggplot(normDatattest_100to300, aes(x = Receptor, y = V1)) +
  geom_point()
#looks here like receptor is peaking around 20

ggplot(normDatattest_100to300, aes(x = SV40, y = V1)) +
  geom_point()
#can't really tell if this is showing us anything

ggplot(normDatattest_100to300, aes(x = CRE, y = V1)) +
  geom_point()
#Cre is getting better as it gets larger

#code data
ttest_100to300.coded<- coded.data(normDatattest_100to300, x1 ~(Receptor-17.5)/7.5, x2 ~ (SV40 - 8.75)/6.25, x3 ~ (CRE - 15)/5)

#rsm analysis
ttest_100to300 <- rsm(V1 ~ FO(x1,x2,x3), data = ttest_100to300.coded)
summary(ttest_100to300)
persp (ttest_100to300, ~ x1 + x2 + x3, col = rainbow(50), contours = "colors") #this is useless


ttest_100to300 <- rsm(V1 ~ SO(x1,x2,x3), data = ttest_100to300.coded)
summary(ttest_100to300)
persp (ttest_100to300, ~ x1 + x2 + x3, col = rainbow(50), contours = "colors")

ttest_100to300 <- rsm(V1 ~ SO(x1,x2,x3) + TWI(x1,x2,x3), data = ttest_100to300.coded)
summary(ttest_100to300)
persp (ttest_100to300, ~ x1 + x2 + x3, col = rainbow(50), contours = "colors")
#these contour graphs don't make sense with what the actual results are which means that the fit is really bad

#the stationary points do make sense though becuase x1 receptor is negative, which means that its at a maximum x2 SV40 is negative which means it has a maximum and CRE is positive which means it only has a minimum
#this means for the next experiment we need to use low values of SV40, similar values of receptor, and high values of CRE

#figure out how to write this up in your lab notebook... also see if this makes sense with the previous data.

