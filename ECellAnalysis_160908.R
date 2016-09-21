#Analysis of the dose response curves done with our top hits 9-8-16
#libraries
library(ggplot2)
library(reshape2)
library(plyr)
library(drc)

#import data
df <- read.csv("/Volumes/mainland/Projects/TAARs/E\ cell\ Optimization/Results/ECellOptimization_160908.csv")
df <- df[,1:11] #get rid of empty columns from the csv file
#create values for normalized luciferase
df$norm <- df$Luc/df$RL
df$logTMA <- log10(df$TMAConcentration)
df$logTMA[which(df$logTMA==-Inf)] = -10
#####HERE I SHOULD ADD SOME CODE TO RUN CHECKS ON THE DATA SO THAT THE RL VALUES ARENT TOO LOW, ETC.

#plot the dose response of the original concentrations - that's from masterMix 12
#plot the dose response curves with the fit line
df.mm.orig <- subset(df, masterMix == 12 & receptorType == "TAAR5")
#the key is that we need to use the non-log form of TMA becuse the equation turns it into the log form on its own
modelTAAR5.orig <- drm(formula = norm ~ TMAConcentration, data = df.mm.orig[complete.cases(df.mm.orig),], fct = LL.4(fixed=c(1,NA,NA,NA),names=(c("Slope", "Top", "Bottom", "ED"))))
summary(modelTAAR5.orig)
df.mm.orig <- subset(df, masterMix == 12 & receptorType == "Rho")
modelRho.orig <- drm(formula = norm ~ TMAConcentration, data = df.mm.orig[complete.cases(df.mm.orig),], fct = LL.4(fixed=c(1,NA,NA,NA),names=(c("Slope", "Top", "Bottom", "ED"))))
summary(modelRho.orig)

# MOD1<-drm(formula = NormalizedLuc ~ concentration, data = x[complete.cases(x),], fct = LL.4(fixed=c(1,NA,NA,NA),names=(c("Slope", "Top", "Bottom", "ED"))))
# plotDR(temp,plot=TRUE)
orig1 = function(x){summary(modelTAAR5.orig)$coefficients[2,1] + (summary(modelTAAR5.orig)$coefficients[1,1]-summary(modelTAAR5.orig)$coefficients[2,1])/(1+10^((log10(summary(modelTAAR5.orig)$coefficients[3,1]) - x)*1))}
orig2 = function(x){summary(modelRho.orig)$coefficients[2,1] + (summary(modelRho.orig)$coefficients[1,1]-summary(modelRho.orig)$coefficients[2,1])/(1+10^((log10(summary(modelRho.orig)$coefficients[3,1]) - x)*1))}

#pdf("/Users/mkamarck/Documents/School\ and\ Lab/Mainland\ Lab/E-Cell-Optimization/Results/160830_h3ADR_plate7.pdf")
ggplot(subset(df, masterMix == 12), aes(x = logTMA, y = norm, colour = factor(receptorType))) +
  geom_point() +
  #geom_errorbar(aes(ymin = average - (standdev/sqrt(6)), ymax = average + (standdev/sqrt(6)))) +
  ylab ("Normalized Luciferase Value") +
  stat_function(fun = orig1) +  
  stat_function(fun = orig2, color="#F8766D") +
  ggtitle (paste("Dose Response to TMA: Receptor=", df$Receptor[1], " SV40=", df$SV40[1], " CRE=", df$CRE[1]))
#dev.off()


#########################
#plot 1
#plot the dose response curves with the fit line
df.mm1 <- subset(df, masterMix == 1 & receptorType == "TAAR5")
#the key is that we need to use the non-log form of TMA becuse the equation turns it into the log form on its own
modelTAAR5 <- drm(formula = norm ~ TMAConcentration, data = df.mm1[complete.cases(df.mm1),], fct = LL.4(fixed=c(1,NA,NA,NA),names=(c("Slope", "Top", "Bottom", "ED"))))
summary(modelTAAR5)
df.mm1 <- subset(df, masterMix == 1 & receptorType == "Rho")
modelRho <- drm(formula = norm ~ TMAConcentration, data = df.mm1[complete.cases(df.mm1),], fct = LL.4(fixed=c(1,NA,NA,NA),names=(c("Slope", "Top", "Bottom", "ED"))))
summary(modelRho)

# MOD1<-drm(formula = NormalizedLuc ~ concentration, data = x[complete.cases(x),], fct = LL.4(fixed=c(1,NA,NA,NA),names=(c("Slope", "Top", "Bottom", "ED"))))
# plotDR(temp,plot=TRUE)
eq1 = function(x){summary(modelTAAR5)$coefficients[2,1] + (summary(modelTAAR5)$coefficients[1,1]-summary(modelTAAR5)$coefficients[2,1])/(1+10^((log10(summary(modelTAAR5)$coefficients[3,1]) - x)*1))}
eq2 = function(x){summary(modelRho)$coefficients[2,1] + (summary(modelRho)$coefficients[1,1]-summary(modelRho)$coefficients[2,1])/(1+10^((log10(summary(modelRho)$coefficients[3,1]) - x)*1))}

#pdf("/Users/mkamarck/Documents/School\ and\ Lab/Mainland\ Lab/E-Cell-Optimization/Results/160830_h3ADR_plate7.pdf")
ggplot(subset(df, masterMix == 1), aes(x = logTMA, y = norm, colour = factor(receptorType))) +
  geom_point() +
  #geom_errorbar(aes(ymin = average - (standdev/sqrt(6)), ymax = average + (standdev/sqrt(6)))) +
  ylab ("Normalized Luciferase Value") +
  stat_function(fun = eq1) +  #, color="#F8766D"
  stat_function(fun = eq2) +
  stat_function(fun = orig1, color = "slateblue") +
  stat_function(fun = orig2, color = "seagreen") +
  ggtitle (paste("Dose Response to TMA: Receptor=", df$Receptor[1], " SV40=", df$SV40[1], " CRE=", df$CRE[1]))
#dev.off()
#########################

#plot loop
#plot the dose response curves with the fit line
#pdf("/Volumes/mainland/Projects/TAARs/E\ cell\ Optimization/Results/160908_DRgraphs.pdf", width = 6, height = 4)
for(i in 1:11){ #11 is for the first 11 master mixes, we are leaving out 12 because it is graphed in each one as a comparison
  rm(df.mm1)
  rm(df.mm2)
  rm(modelTAAR5)
  rm(modelRho)
  rm(eq1)
  rm(eq2)
  print(i)
  df.mm1 <- subset(df, masterMix == i & receptorType == "TAAR5")
  #the key is that we need to use the non-log form of TMA becuse the equation turns it into the log form on its own
  #make the models
  tryCatch({
    modelTAAR5 <- drm(formula = norm ~ TMAConcentration, data = df.mm1[complete.cases(df.mm1),], fct = LL.4(fixed=c(1,NA,NA,NA),names=(c("Slope", "Top", "Bottom", "ED"))))
    summary(modelTAAR5)
  }, error=function(e){cat("ERROR:", conditionMessage(e), "\n")})
  df.mm2 <- subset(df, masterMix == i & receptorType == "Rho")
  tryCatch({
    modelRho <- drm(formula = norm ~ TMAConcentration, data = df.mm2[complete.cases(df.mm2),], fct = LL.4(fixed=c(1,NA,NA,NA),names=(c("Slope", "Top", "Bottom", "ED"))))
    #summary(modelRho) 
  }, error=function(e){cat("ERROR:", conditionMessage(e), "\n")})
  
  #set the equations equal to a variable so we can graph them later
  eq1 = function(x){summary(modelTAAR5)$coefficients[2,1] + (summary(modelTAAR5)$coefficients[1,1]-summary(modelTAAR5)$coefficients[2,1])/(1+10^((log10(summary(modelTAAR5)$coefficients[3,1]) - x)*1))}
  eq2 = function(x){summary(modelRho)$coefficients[2,1] + (summary(modelRho)$coefficients[1,1]-summary(modelRho)$coefficients[2,1])/(1+10^((log10(summary(modelRho)$coefficients[3,1]) - x)*1))}
  #graph the data
  print(ggplot(subset(df, masterMix == i), aes(x = logTMA, y = norm, colour = factor(receptorType))) +
    geom_point() +
    #geom_errorbar(aes(ymin = average - (standdev/sqrt(6)), ymax = average + (standdev/sqrt(6)))) +
    ylab ("Normalized Luciferase Value") +
    stat_function(fun = eq1) +  
    stat_function(fun = eq2, color="#F8766D" ) +
    stat_function(fun = orig1, color = "slateblue") +
    stat_function(fun = orig2, color = "seagreen") +
    ggtitle (paste("Dose Response to TMA: Master Mix=", i, " Receptor=", subset(df, masterMix == i)$Receptor[1], " SV40=", subset(df, masterMix == i)$SV40[1], " CRE=", subset(df, masterMix == i)$CRE[1]))
  )
    
}
#dev.off()
#master mix 3 couldn't fit a function for TAAR5
##################################################
#plotting the same stuff but with luciferase values instead of normalized luciferase values
#luc models for the original
df.mm.orig <- subset(df, masterMix == 12 & receptorType == "TAAR5")
#the key is that we need to use the non-log form of TMA becuse the equation turns it into the log form on its own
modelTAAR5.orig.luc <- drm(formula = Luc ~ TMAConcentration, data = df.mm.orig[complete.cases(df.mm.orig),], fct = LL.4(fixed=c(1,NA,NA,NA),names=(c("Slope", "Top", "Bottom", "ED"))))
summary(modelTAAR5.orig.luc)
df.mm.orig <- subset(df, masterMix == 12 & receptorType == "Rho")
modelRho.orig.luc <- drm(formula = Luc ~ TMAConcentration, data = df.mm.orig[complete.cases(df.mm.orig),], fct = LL.4(fixed=c(1,NA,NA,NA),names=(c("Slope", "Top", "Bottom", "ED"))))
summary(modelRho.orig.luc)

# MOD1<-drm(formula = NormalizedLuc ~ concentration, data = x[complete.cases(x),], fct = LL.4(fixed=c(1,NA,NA,NA),names=(c("Slope", "Top", "Bottom", "ED"))))
# plotDR(temp,plot=TRUE)
orig1.luc = function(x){summary(modelTAAR5.orig.luc)$coefficients[2,1] + (summary(modelTAAR5.orig.luc)$coefficients[1,1]-summary(modelTAAR5.orig.luc)$coefficients[2,1])/(1+10^((log10(summary(modelTAAR5.orig.luc)$coefficients[3,1]) - x)*1))}
orig2.luc = function(x){summary(modelRho.orig.luc)$coefficients[2,1] + (summary(modelRho.orig.luc)$coefficients[1,1]-summary(modelRho.orig.luc)$coefficients[2,1])/(1+10^((log10(summary(modelRho.orig.luc)$coefficients[3,1]) - x)*1))}

#pdf("/Volumes/mainland/Projects/TAARs/E\ cell\ Optimization/Results/160908_DRgraphs_luc.pdf", width = 6, height = 4)
for(i in 1:11){ #11 is for the first 11 master mixes, we are leaving out 12 because it is graphed in each one as a comparison
  rm(df.mm1)
  rm(df.mm2)
  rm(modelTAAR5)
  rm(modelRho)
  rm(eq1)
  rm(eq2)
  print(i)
  df.mm1 <- subset(df, masterMix == i & receptorType == "TAAR5")
  #the key is that we need to use the non-log form of TMA becuse the equation turns it into the log form on its own
  #make the models
  tryCatch({
    modelTAAR5 <- drm(formula = Luc ~ TMAConcentration, data = df.mm1[complete.cases(df.mm1),], fct = LL.4(fixed=c(1,NA,NA,NA),names=(c("Slope", "Top", "Bottom", "ED"))))
    summary(modelTAAR5)
  }, error=function(e){cat("ERROR:", conditionMessage(e), "\n")})
  df.mm2 <- subset(df, masterMix == i & receptorType == "Rho")
  tryCatch({
    modelRho <- drm(formula = Luc ~ TMAConcentration, data = df.mm2[complete.cases(df.mm2),], fct = LL.4(fixed=c(1,NA,NA,NA),names=(c("Slope", "Top", "Bottom", "ED"))))
    #summary(modelRho) 
  }, error=function(e){cat("ERROR:", conditionMessage(e), "\n")})
  
  #set the equations equal to a variable so we can graph them later
  eq1 = function(x){summary(modelTAAR5)$coefficients[2,1] + (summary(modelTAAR5)$coefficients[1,1]-summary(modelTAAR5)$coefficients[2,1])/(1+10^((log10(summary(modelTAAR5)$coefficients[3,1]) - x)*1))}
  eq2 = function(x){summary(modelRho)$coefficients[2,1] + (summary(modelRho)$coefficients[1,1]-summary(modelRho)$coefficients[2,1])/(1+10^((log10(summary(modelRho)$coefficients[3,1]) - x)*1))}
  #graph the data
  print(ggplot(subset(df, masterMix == i), aes(x = logTMA, y = Luc, colour = factor(receptorType))) +
          geom_point() +
          #geom_errorbar(aes(ymin = average - (standdev/sqrt(6)), ymax = average + (standdev/sqrt(6)))) +
          ylab ("Normalized Luciferase Value") +
          stat_function(fun = eq1) +  
          stat_function(fun = eq2, color="#F8766D" ) +
          stat_function(fun = orig1.luc, color = "slateblue") + 
          stat_function(fun = orig2.luc, color = "seagreen") +
          ggtitle (paste("Dose Response to TMA: Master Mix=", i, " Receptor=", subset(df, masterMix == i)$Receptor[1], " SV40=", subset(df, masterMix == i)$SV40[1], " CRE=", subset(df, masterMix == i)$CRE[1]))
  )
}
#dev.off()
#this doesn't give us that much extra information, but it helps visualize what may be the best concentrations to use. 
###############################################
#signal to noise comparisons
#calculate the signal to noise by dividing the mean by the standard deviation
df.STN <- ddply(.data = subset(df, receptorType == "TAAR5"), .variables = c("masterMix", "receptorType", "TMAConcentration"), .fun = summarize, signal = mean(norm), noise = sd(norm))
df.STN$STN <- df.STN$signal/df.STN$noise
#make a variable with the TMA Concentrations
TMAConcentration <- unique(df$TMAConcentration)
masterMixes = 13 #defines the length of our matrix
#create two dataframes to compare the signal to noise values
order.mm <- matrix(ncol=length(TMAConcentration), nrow=masterMixes)
order.STN <- matrix(ncol=length(TMAConcentration), nrow=masterMixes)
for(i in 1:length(TMAConcentration)){
  print(TMAConcentration[i])
  df.STN.TMA <- subset(df.STN, df.STN$TMAConcentration == TMAConcentration[i]) #creates a subset of the database with a specific TMA concentration
  df.STN.sort <- df.STN.TMA[order(df.STN.TMA$STN, decreasing = TRUE),] #orders the STN at that concentration from highest to lowest
  order.mm[,i] <- df.STN.TMA[order(df.STN.TMA$STN, decreasing = TRUE),]$masterMix #write the order of the master mix numbers to the df; had to write it like this, otherwise it changed the order of the mastermix and I can't figure out why
  order.STN[,i] <- df.STN.sort$STN #write the order of the master mix numbers to the df
}

order.mm <- as.data.frame(order.mm)
order.STN <- as.data.frame(order.STN)
names(order.mm) <- TMAConcentration
names(order.STN) <- TMAConcentration

write.csv(c(order.mm, order.STN), file = "/Volumes/mainland/Projects/TAARs/E\ cell\ Optimization/Results/STNResults_160908.csv")
###############################################
#I can do a similar analysis of the signal to noise with p values, once I figure out how to get the p values from the models

###############################################

#GOALS
#1. graph this in comparison to the original - check
#2. figure out how to do statistics between the different models
#3. write this graph and statistics into a loop so I can do it for all the graphs - check
#4. Look at the H3A data and try to eliminate wells that didn't work
#5. compare signal to noise between the different options; look at average signal to noise and also at signal to noise at each of the various points



