#analyze the data from 9-30-16
#this was an antagonist screen comparing results between H3A and E cells, and to make sure that the same top hits came up as the last time we ran the screen. 
#import
library(tidyverse)
library(drc)
#import data
plate1 <- read.csv("/Volumes/mainland/Projects/TAARs/E\ cell\ Optimization/AntagScreenRepeat_mini/160930_plate1.csv", stringsAsFactors = FALSE)[,1:11]
plate2 <- read.csv("/Volumes/mainland/Projects/TAARs/E\ cell\ Optimization/AntagScreenRepeat_mini/160930_plate2.csv", stringsAsFactors = FALSE)[,1:11]
plate3 <- read.csv("/Volumes/mainland/Projects/TAARs/E\ cell\ Optimization/AntagScreenRepeat_mini/160930_plate3.csv", stringsAsFactors = FALSE)[,1:11]
plate4 <- read.csv("/Volumes/mainland/Projects/TAARs/E\ cell\ Optimization/AntagScreenRepeat_mini/160930_plate4.csv", stringsAsFactors = FALSE)[,1:11]
#combine per cell type
eCells <- rbind(plate1, plate2)
H3A <- rbind(plate3, plate4)
#calculate the normalized response
eCells$norm <- eCells$Luc/eCells$RL
H3A$norm <- H3A$Luc/H3A$RL
###########DOSE RESPONSE###############################
#becuase of all the TMA in 384 well plate, we ran a dose response to see if we are actually applying the EC80 to cells.  
#In other words, if there is enough bleedover between wells, we might be seeing more like 600 or more uM TMA applied to cells than 500uM
#therefore we ran a dose response in the middle of the plate to see if the curve looked the same as we were expecting. 

#E-CELLS
DRTMA.eCells <- subset(eCells, experimentType == "DR")
DRTMA.eCells$logTMA <- log10(DRTMA.eCells$concentration_ag) #take logTMA for graphing (not for fitting the function)
DRTMA.eCells$concentration_ag[which(DRTMA.eCells$logTMA==-10)] = 0
#create the fit line - TAAR5
TAAR5.eCells <- subset(DRTMA.eCells, CloneKey == 830)
#the key is that we need to use the non-log form of TMA becuse the equation turns it into the log form on its own
modelTAAR5 <- drm(formula = norm ~ concentration_ag, data = TAAR5.eCells, fct = LL.4(fixed=c(1,NA,NA,NA),names=(c("Slope", "Top", "Bottom", "ED"))))
summary(modelTAAR5)
TAAR5.eq = function(x){summary(modelTAAR5)$coefficients[2,1] + (summary(modelTAAR5)$coefficients[1,1]-summary(modelTAAR5)$coefficients[2,1])/(1+10^((log10(summary(modelTAAR5)$coefficients[3,1]) - x)*1))}
#fit for rho
Rho.eCells <- subset(DRTMA.eCells, CloneKey == 999)
modelRho <- drm(formula = norm ~ concentration_ag, data = Rho.eCells, fct = LL.4(fixed=c(1,NA,NA,NA),names=(c("Slope", "Top", "Bottom", "ED"))))
summary(modelRho)
Rho.eq = function(x){summary(modelRho)$coefficients[2,1] + (summary(modelRho)$coefficients[1,1]-summary(modelRho)$coefficients[2,1])/(1+10^((log10(summary(modelRho)$coefficients[3,1]) - x)*1))}
#plot TMA DR
ggplot(DRTMA.eCells, aes(x = logTMA, y = norm, colour = factor(CloneKey))) +
  geom_point() +
  #geom_errorbar(aes(ymin = average - (standdev/sqrt(6)), ymax = average + (standdev/sqrt(6)))) +
  ylab ("Normalized Luciferase Value") +
  stat_function(fun = TAAR5.eq, colour = "#F8766D") +  
  stat_function(fun = Rho.eq, colour = "#00BFC4")
#that's better, but the function still isn't fitting well.  - maybe I need to try putting this in prism
library(scales)
show_col(hue_pal()(2))

#H3A
DR.H3A <- subset(H3A, experimentType == "DR")
DR.H3A$logTMA <- log10(DR.H3A$concentration_ag) #take logTMA for graphing (not for fitting the function)
DR.H3A$concentration_ag[which(DR.H3A$logTMA==-10)] = 0
TAAR5.H3A <- subset(DR.H3A, CloneKey == 830)
Rho.H3A <- subset(DR.H3A, CloneKey == 999)
#model
modelTAAR5 <- drm(formula = norm ~ concentration_ag, data = TAAR5.H3A, fct = LL.4(fixed=c(1,NA,NA,NA),names=(c("Slope", "Top", "Bottom", "ED"))))
summary(modelTAAR5)
TAAR5.eq = function(x){summary(modelTAAR5)$coefficients[2,1] + (summary(modelTAAR5)$coefficients[1,1]-summary(modelTAAR5)$coefficients[2,1])/(1+10^((log10(summary(modelTAAR5)$coefficients[3,1]) - x)*1))}
#fit for rho
modelRho <- drm(formula = norm ~ concentration_ag, data = Rho.H3A, fct = LL.4(fixed=c(1,NA,NA,NA),names=(c("Slope", "Top", "Bottom", "ED"))))
summary(modelRho)
Rho.eq = function(x){summary(modelRho)$coefficients[2,1] + (summary(modelRho)$coefficients[1,1]-summary(modelRho)$coefficients[2,1])/(1+10^((log10(summary(modelRho)$coefficients[3,1]) - x)*1))}
#plot TMA DR
ggplot(DR.H3A, aes(x = logTMA, y = norm, colour = factor(CloneKey))) +
  geom_point() +
  #geom_errorbar(aes(ymin = average - (standdev/sqrt(6)), ymax = average + (standdev/sqrt(6)))) +
  ylab ("Normalized Luciferase Value") +
  stat_function(fun = TAAR5.eq, colour = "#F8766D") +  
  stat_function(fun = Rho.eq, colour = "#00BFC4")

#This second one is weird - why is the H3A response so low? 
#why is the forskolin response at baseline so high?
#why dn't the forskolin and the TMA responses match?
#how then should I analyze this? 
