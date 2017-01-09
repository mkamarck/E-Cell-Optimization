#analyze the data from 9-30-16
#this was an antagonist screen comparing results between H3A and E cells, and to make sure that the same top hits came up as the last time we ran the screen. 
#import
library(tidyverse)
library(drc)
#import data
data <- read.csv("/Volumes/mainland/Projects/TAARs/E\ cell\ Optimization/AntagScreenRepeat_mini/161104_Results.csv", stringsAsFactors = FALSE)
#get rid of unnecessary columns
#gives all data
data <- subset(data, select = c("TransfectionDate", "PlateNum", "CloneKey", "Odor1K", "Odor2K","Odor1C", "Odor2C", "PlateLocation", "LucData", "RLData", "Notes","LabNoteBook", "TransfectionType", "CellType"))
#cuts down more columns to the useful ones for this analysis
data.sub <- subset(data, select = c("PlateNum", "CloneKey", "Odor1K", "Odor2K","Odor1C", "Odor2C", "PlateLocation", "LucData", "RLData", "CellType"))

#calculate normalized luc response
data.sub$norm <- data.sub$LucData/data.sub$RLData
#combine per cell type
eCells <- subset(data.sub, PlateNum == 1 | PlateNum == 2)
H3A <- subset(data.sub, PlateNum == 3 | PlateNum == 4)

###########DOSE RESPONSE###############################
#becuase of all the TMA in 384 well plate, we ran a dose response to see if we are actually applying the EC80 to cells.  
#In other words, if there is enough bleedover between wells, we might be seeing more like 600 or more uM TMA applied to cells than 500uM
#therefore we ran a dose response in the middle of the plate to see if the curve looked the same as we were expecting. 

#concentration of TMA alone
TMA_conc <- log10(500e-6)
for_conc <- log10(1e-6)

#E-CELLS
DRTMA.eCells <- subset(eCells, (CloneKey == 830 & (Odor1C !=500e-6 & Odor1C!=1)) | (CloneKey == 999 & (Odor1C != 1e-6 & Odor1C !=1)))
DRTMA.eCells$logTMA <- log10(DRTMA.eCells$Odor1C) #take logTMA for graphing (not for fitting the function)
DRTMA.eCells$logTMA[which(DRTMA.eCells$Odor1C == 0)] = -10
#create the fit line - TAAR5
TAAR5.eCells <- subset(DRTMA.eCells, CloneKey == 830)
#the key is that we need to use the non-log form of TMA becuse the equation turns it into the log form on its own
modelTAAR5 <- drm(formula = norm ~ Odor1C, data = TAAR5.eCells, fct = LL.4(fixed=c(1,NA,NA,NA),names=(c("Slope", "Top", "Bottom", "ED"))))
summary(modelTAAR5)
TAAR5.eq = function(x){summary(modelTAAR5)$coefficients[2,1] + (summary(modelTAAR5)$coefficients[1,1]-summary(modelTAAR5)$coefficients[2,1])/(1+10^((log10(summary(modelTAAR5)$coefficients[3,1]) - x)*1))}
####fit for rho on graph#####
#fit for rho
Rho.eCells <- subset(DRTMA.eCells, CloneKey == 999)
modelRho <- drm(formula = norm ~ Odor1C, data = Rho.eCells, fct = LL.4(fixed=c(1,NA,NA,NA),names=(c("Slope", "Top", "Bottom", "ED"))))
summary(modelRho)
Rho.eq = function(x){summary(modelRho)$coefficients[2,1] + (summary(modelRho)$coefficients[1,1]-summary(modelRho)$coefficients[2,1])/(1+10^((log10(summary(modelRho)$coefficients[3,1]) - x)*1))}
#plot TMA DR
  #Agx1 is from E2 - is TMA with CD293 (1ul) so it should be about equal to the TMA on the line, but its lower - if anything I was expecting it to be higher
  #code to put forskolin by itself on the ecell graph
#   Rho <- subset(eCells, PlateNum == 1 & Odor1K == 77)
#   LowRL<-subset(Rho, Rho$RL < 1000) #how many have low RL?
  #none
  #RhoAlone<-subset(Rho, antagInStockKey == 811) ##Agonist only
######graph###########
#pdf("/Volumes/mainland-1/Projects/TAARs/E\ cell\ Optimization/AntagScreenRepeat_mini/TMA_DR_ECells.pdf")
ggplot(DRTMA.eCells, aes(x = logTMA, y = norm, colour = factor(CloneKey))) +
  geom_point() +
  #geom_point(data = Agx1, aes(x = TMA_conc, y = norm, colour = "blue")) +
  #geom_point(data = RhoAlone, aes(x = for_conc, y = norm, colour = "red")) +
  #geom_errorbar(aes(ymin = average - (standdev/sqrt(6)), ymax = average + (standdev/sqrt(6)))) +
  ylab ("Normalized Luciferase Value") +
  stat_function(fun = TAAR5.eq, colour = "#F8766D") +  
  stat_function(fun = Rho.eq, colour = "#00BA38") 
#dev.off()

#For 11-4-16; I think I didnt do a dR but just did TMA by itself in all of the wells that didn't say CD293 - I'm not sure but this certainly didn't work here
#I'm guessing I also only put CD293 in the well that has the lowest values but claims to have some amount of TMA. 


# library(scales)
# show_col(hue_pal()(4))
#  ggplot(TAAR5.eCells, aes(x = logTMA, y = norm)) +
#    geom_point() +
#    stat_function(fun = TAAR5.eq)

#H3A
DR.H3A <- subset(H3A, (CloneKey == 830 & (Odor1C !=500e-6 & Odor1C!=1)) | (CloneKey == 999 & (Odor1C != 1e-6 & Odor1C !=1)))
DR.H3A$logTMA <- log10(DR.H3A$Odor1C) #take logTMA for graphing (not for fitting the function)
DR.H3A$logTMA[which(DRTMA.eCells$Odor1C == 0)] = -10
TAAR5.H3A <- subset(DR.H3A, CloneKey == 830)
Rho.H3A <- subset(DR.H3A, CloneKey == 999)
#model
modelTAAR5 <- drm(formula = norm ~ Odor1C, data = TAAR5.H3A, fct = LL.4(fixed=c(1,NA,NA,NA),names=(c("Slope", "Top", "Bottom", "ED"))))
summary(modelTAAR5)
TAAR5.eq = function(x){summary(modelTAAR5)$coefficients[2,1] + (summary(modelTAAR5)$coefficients[1,1]-summary(modelTAAR5)$coefficients[2,1])/(1+10^((log10(summary(modelTAAR5)$coefficients[3,1]) - x)*1))}
#fit for rho
# modelRho <- drm(formula = norm ~ Odor1C, data = Rho.H3A, fct = LL.4(fixed=c(1,NA,NA,NA),names=(c("Slope", "Top", "Bottom", "ED"))))
# summary(modelRho)
# Rho.eq = function(x){summary(modelRho)$coefficients[2,1] + (summary(modelRho)$coefficients[1,1]-summary(modelRho)$coefficients[2,1])/(1+10^((log10(summary(modelRho)$coefficients[3,1]) - x)*1))}
# #calculate EC50
# ED(modelRho, c(50, 80), interval = "delta")


#plot TMA DR
#pdf("/Volumes/mainland-1/Projects/TAARs/E\ cell\ Optimization/AntagScreenRepeat_mini/TMA_DR_H3A.pdf")
ggplot(DR.H3A, aes(x = logTMA, y = norm, colour = factor(CloneKey))) +
  geom_point() +
  #geom_errorbar(aes(ymin = average - (standdev/sqrt(6)), ymax = average + (standdev/sqrt(6)))) +
  ylab ("Normalized Luciferase Value") +
  stat_function(fun = TAAR5.eq, colour = "#F8766D") #+  
  #stat_function(fun = Rho.eq, colour = "#00BFC4")
#dev.off()

#This one is weird becuase it kindof looks like a dose response - is it supposed to be ?

#I don't know what's going on here...
