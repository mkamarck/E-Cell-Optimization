#Analyze E cell optimization - first run 160719 (7-19-16)

#libraries
library(ggplot2)
library(reshape2)
library(plyr)
library(rsm)

#import data
df <- read.csv("/Volumes/mainland/Projects/TAARs/E\ cell\ Optimization/Results/ECellOptimization_160719.csv")
df <- df[,1:14]

#Goal: conduct a bunch of t-tests to see how different the different concentraitons of TMA look with different concentrations of the vectors
#We want to find the combination of vectors that is the best at telling apart the different concentrations of TMA
#make dataframes for each of the three kinds of data

lucData <- df[,1:8]
RLData <- df[,c(1:5, 9:11)]
normData <- df[,c(1:5, 12:14)]
########################
#start by analyzing normalized data - none of this is working...
#I need to get p values for how they are different at the different concentrations, but I can't get it to work, so I'm going to come back
 normData.melt <- melt(normData, c("Row", "Column", "Receptor", "SV40", "CRE"))
 normDatattest <- ddply(normData.melt, .variables = c("Receptor", "SV40", "CRE"), function(x) tests = summary(aov(x$value~x$variable))[[1]]$'Pr(>F)'[1])
# 
# 
# test =summary(aov(normData.melt,value~variable))[[1]][["Pr(>F)"]]
# test = summary(aov(normData.melt,value~variable))[[1]]$'Pr(>F)'
# 
# est <- aov(value~variable, data = normData.melt)
# summary(est)
# summary(est)[[1]]$'Pr(>F)'[1]
# test <- summary(aov(value~variable, data = normData.melt))[[1]]$'Pr(>F)'[1]
# 
# normData$test30_100 <- t.test(TMA_30_norm ~TMA_100_norm, data = normData)
##################################################################

#Start by analyzing using rsm package to maximize the normData response. I will start with TMA_100
normData_TMA100 <- subset(normData, select = c("Receptor", "SV40", "CRE", "TMA_100_norm"))
#code the receptor, sv40 and CRE into equally spaced variables, example 
#coded.data(ChemReact1, x1 ~ (Time - 85)/5, x2 ~ (Temp - 175)/5)
normTMA_100_coded <- coded.data(normData_TMA100, x1 ~(Receptor-15)/10, x2 ~ (SV40 - 8.75)/6.25, x3 ~ (CRE - 11.25)/8.75)
#try and fit models
rsm_testFO <- rsm(TMA_100_norm ~FO(x1,x2,x3), data = normTMA_100_coded)
summary(rsm_testFO)
#lack of fit is non-significant, but I would think that a second order equation would be better? I don't know how you decide which one to choose
rsm_testSO <- rsm(TMA_100_norm ~SO(x1,x2,x3), data = normTMA_100_coded)
summary(rsm_testSO)
#two way interaction adds onto the normal FO rsm test, which seems to fit fine - I can also try this with SO
rsm_testTWI <- rsm(TMA_100_norm ~ FO(x1,x2,x3) + TWI(x1,x2,x3), data = normTMA_100_coded)
summary(rsm_testTWI)
#SO and TWI
rsm_testTWISO <- rsm(TMA_100_norm ~ SO(x1,x2,x3) + TWI(x1,x2,x3), data = normTMA_100_coded)
summary(rsm_testTWISO)
#how do I know whether to use a first or second order - I'm pretty sure the two way interaction makes sense with all of them, but why then does the receptor value not seem to matter
#we use second order when we think the optimal response is somewhere in our response surface
#with FO we use steepest ascent to decide in what direction to move with our response surface
#lets try some visualization
#example
#contour(heli.rsm, ~ x1 + x2 + x3 + x4, image = TRUE,
 #       +   at = summary(heli.rsm)$canonical$xs)
contour(rsm_testTWISO, ~x1+x2+x3, image = TRUE)
#the hard part about this normalization part is that we have extra hidden factors as we truely have 2 response outputs rather than 1 that we are combining into one
#we may want to analyze the response surface differently so as to include luc and RL values

####################################################################
#lets try just looking at the luc values with TMA_100
lucData_TMA100 <- subset(lucData, select = c("Receptor", "SV40", "CRE", "TMA_100_luc"))
#code it
lucTMA_100_coded <- coded.data(lucData_TMA100, x1 ~(Receptor-15)/10, x2 ~ (SV40 - 8.75)/6.25, x3 ~ (CRE - 11.25)/8.75)
#FO, look at steepest ascent
rsm_testFO <- rsm(TMA_100_luc ~ FO(x1,x2,x3), data = lucTMA_100_coded)
summary(rsm_testFO)
#Fo TWI
rsm_testTWI <- rsm(TMA_100_luc ~ FO(x1,x2,x3) + TWI(x1,x2,x3), data = lucTMA_100_coded)
summary(rsm_testTWI)
#this test makes no sense because it seems to be giving me a minimum instead of a maximum
#SO TWI
rsm_testTWI <- rsm(TMA_100_luc ~ SO(x1,x2,x3) + TWI(x1,x2,x3), data = lucTMA_100_coded)
summary(rsm_testTWI)
#SO TEST IS not appropriate here

#visualization
contour(rsm_testTWI, ~x1+x3, image = TRUE)
#maybe the norm luc is right and the receptor concentration really doesn't matter...
#try this with different TMA concetrations
lucData_TMA30 <- subset(lucData, select = c("Receptor", "SV40", "CRE", "TMA_30_luc"))
lucTMA_30_coded <- coded.data(lucData_TMA30, x1 ~(Receptor-15)/10, x2 ~ (SV40 - 8.75)/6.25, x3 ~ (CRE - 11.25)/8.75)
#FO test
rsm_testFO <- rsm(TMA_30_luc ~ FO(x1,x2,x3), data = lucTMA_30_coded)
summary(rsm_testFO)
#Fo TWI
rsm_testTWI <- rsm(TMA_30_luc ~ FO(x1,x2,x3) + TWI(x1,x2,x3), data = lucTMA_30_coded)
summary(rsm_testTWI)
contour(rsm_testTWI, ~x1+x3, image = TRUE)

#this is starting to make sense - the FO gives direction of steepest ascent which will bring us to a maximum point - we cannot really get a maximum simply by looking at luc
#because as much as you increase CRE and receptor you will get higher values. 
#so what we really want to maximize is difference between different concentraitons, so we need to try a whole different way of analysis
##################################
#I wanted to try doing this analysis on p-values, but I can't get R to behave, or I can't figure out how to do the appropriate t-tests. 
normData.melt <- melt(normData, c("Row", "Column", "Receptor", "SV40", "CRE"))
normData.melt <- normData.melt[,3:7]
normDatattest_30to100 <- ddply(subset(normData.melt, variable %in% c("TMA_30_norm", "TMA_100_norm")), .variables = c("Receptor", "SV40", "CRE"), function(x) test = t.test(x$value~x$variable)$p.value)

exampledf <- subset(normData.melt, variable %in% c("TMA_30_norm", "TMA_100_norm"))
test <- t.test(value~variable, data = exampledf)$p.value
summary(test)
##################################
#lets just make a database with differences and then run the RSM analysis on it.
#we'll start with normalization data
#went back and added a PQ term - and now it seems to make more sense
normData$diffTMA_30to100 <- normData$TMA_100_norm - normData$TMA_30_norm
normData$diffTMA_100to300 <- normData$TMA_300_norm - normData$TMA_100_norm
normData$diffTMA_30to300 <- normData$TMA_300_norm - normData$TMA_30_norm
#code normData
normData.coded <- coded.data(normData, x1 ~(Receptor-15)/10, x2 ~ (SV40 - 8.75)/6.25, x3 ~ (CRE - 11.25)/8.75)
#FO, look at steepest ascent
diff_30to300_FO <- rsm(diffTMA_30to300 ~ FO(x1,x2,x3), data = normData.coded)
summary(diff_30to300_FO)
contour(diff_30to300_FO, ~x1+x2+x3, image = TRUE)

diff_30to100_FO <- rsm(diffTMA_30to100 ~ FO(x1,x2,x3), data = normData.coded)
summary(diff_30to100_FO)
contour(diff_30to100_FO, ~x1+x2, image = TRUE)

diff_100to300_FO <- rsm(diffTMA_100to300 ~ FO(x1,x2,x3), data = normData.coded)
summary(diff_100to300_FO)
contour(diff_100to300_FO, ~x2+x3, image = TRUE)
#try first order with interaction
diff_30to300_TWI <- rsm(diffTMA_30to300 ~ FO(x1,x2,x3) + TWI(x1,x2,x3), data = normData.coded)
summary(diff_30to300_TWI)
contour(diff_30to300_TWI, ~x1+x2+x3, image = TRUE)

diff_30to100_TWI <- rsm(diffTMA_30to100 ~ FO(x1,x2,x3) + TWI(x1,x2,x3), data = normData.coded)
summary(diff_30to100_TWI)
contour(diff_30to100_TWI, ~x1+x2+x3, image = TRUE)

diff_100to300_TWI <- rsm(diffTMA_100to300 ~ FO(x1,x2,x3) + TWI(x1,x2,x3), data = normData.coded)
summary(diff_100to300_TWI)
contour(diff_100to300_TWI, ~x1+x2+x3, image = TRUE)
#from these last two it seems like it might be that CRE and SV40 are converging around 10-15ng but its not clear whether the receptor is having any effect at all.

diff_100to300_TWI <- rsm(diffTMA_100to300 ~ FO(x1,x2,x3) + TWI(x1,x2,x3) + PQ(x1,x2,x3), data = normData.coded)
summary(diff_100to300_TWI)
contour(diff_100to300_TWI, ~x1+x2+x3, image = TRUE)

diff_30to100_TWI <- rsm(diffTMA_30to100 ~ FO(x1,x2,x3) + TWI(x1,x2,x3)+ PQ(x1,x2,x3), data = normData.coded)
summary(diff_30to100_TWI)
contour(diff_30to100_TWI, ~x1+x2+x3, image = TRUE)

diff_30to300_TWI <- rsm(diffTMA_30to300 ~ FO(x1,x2,x3) + TWI(x1,x2,x3)+ PQ(x1,x2,x3), data = normData.coded)
summary(diff_30to300_TWI)
contour(diff_30to300_TWI, ~x1+x2+x3, image = TRUE)
#when I add the PQ into the model it seems to make more sense and be converging around receptor at 15, SV40 at 10 and CRE at 12 - these are close to our original values but with the receptor higher
#of course for this the eigenanalysis is not a maximum for all three points, so I'm not sure whether it is optimal or we still need to work on values that are far away.

####################################
#I am going to do one last analysis with the luc data to see if the answers are different
lucData$diffTMA_30to100 <- lucData$TMA_100_luc - lucData$TMA_30_luc
lucData$diffTMA_100to300 <- lucData$TMA_300_luc - lucData$TMA_100_luc
lucData$diffTMA_30to300 <- lucData$TMA_300_luc - lucData$TMA_30_luc
lucData.coded <- coded.data(lucData, x1 ~(Receptor-15)/10, x2 ~ (SV40 - 8.75)/6.25, x3 ~ (CRE - 11.25)/8.75)

diff_30to300_FO <- rsm(diffTMA_30to300 ~ FO(x1,x3), data = lucData.coded)
summary(diff_30to300_FO)
contour(diff_30to300_FO, ~x1+x3, image = TRUE)

diff_30to100_FO <- rsm(diffTMA_30to100 ~ FO(x1,x3), data = lucData.coded)
summary(diff_30to100_FO)
contour(diff_30to100_FO, ~x1+x3, image = TRUE)

diff_100to300_FO <- rsm(diffTMA_100to300 ~ FO(x1,x3), data = lucData.coded)
summary(diff_100to300_FO)
contour(diff_100to300_FO, ~x1+x3, image = TRUE)

#try first order with interaction
diff_30to300_TWI <- rsm(diffTMA_30to300 ~ FO(x1,x3) + TWI(x1,x3), data = lucData.coded)
summary(diff_30to300_TWI)
contour(diff_30to300_TWI, ~x1+x2+x3, image = TRUE)
#these aren't different enough between the different groups to get any useful analysis from it

diff_30to100_TWI <- rsm(diffTMA_30to100 ~ FO(x1,x3) + TWI(x1,x3), data = lucData.coded)
summary(diff_30to100_TWI)
contour(diff_30to100_TWI, ~x1+x3, image = TRUE)

diff_100to300_TWI <- rsm(diffTMA_100to300 ~ FO(x1,x3) + TWI(x1,x3), data = lucData.coded)
summary(diff_100to300_TWI)
contour(diff_100to300_TWI, ~x1+x3, image = TRUE)

diff_100to300_TWI <- rsm(diffTMA_100to300 ~ FO(x1,x3) + TWI(x1,x3) + PQ(x1,x3), data = lucData.coded)
summary(diff_100to300_TWI)
contour(diff_100to300_TWI, ~x1+x3, image = TRUE)


#This analysis wasn't very helpful becuase its really only useful for doing x1 and x3 - try that.
#when I try that with first order the point of steepest ascent is right around where we normally have our values

#I think to really do this analysis right I still need to figure out how to conduct it with multiple outcomes. 
#It could also be that because I only ran 1 sample of each, if I did p values that would help
#the other thing that might be useful is to redesign the experiment - i didn't design the experiment to sample the response space very well, some parts of the space are sampled better than others and that could be skewing the results. 

#I could also try including the concentration of TMA as a variable and seeing if it comes out correctly with the model - this way it may almost work as a control becuase I know where the top response of TMA concentration should be

