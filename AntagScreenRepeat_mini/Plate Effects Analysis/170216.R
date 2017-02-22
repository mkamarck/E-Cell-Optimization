#Analyzing for plate effects due to degredation of signal and timing of adding luc and RL 
library(platetools)
library(ggplot2)
library(viridis)
library(plyr)
library(reshape2)
library(pipeR)
#read in data
df <- read.csv("170216_plateEffects.csv") %>>% 
  subset(select = c(PlateNum, CloneKey, Odor1K, Odor1C, PlateLocation, RLData, LucData, Notes, RunNum, RunType, expecting.Lum))

#plates 1 and 2 have additional information, add this information about the other plates#####

df$RunNum <- as.factor(df$RunNum)
df$RunType <- as.character(df$RunType)
for (i in 1:length(df$RunNum)){
  if(df$RunNum[i] %in% c(1,3,7,9)) {
    df$RunType[i] = "Luc"
  } 
  if(df$RunNum[i] %in% c(2,5,8,11)){
   df$RunType[i] = "RL1" 
  } 

  if(df$RunNum[i] %in% c(4,6,10,13)){
    df$RunType[i] = "RL2"
  }
  if(df$RunNum[i] %in% c(12, 15)){
    df$RunType[i] = "RL3"
  }
  if(df$RunNum[i] %in% c(14, 16)){
    df$RunType[i] = "RL4"
  }
}
for (i in 1:length(df$RunNum)){
  if(df$RunNum[i] %in% c(1,3,7:16)){
    df$expecting.Lum[i] = "y"
  }
  else
    df$expecting.Lum[i] = "n"
}


for (i in 1:length(df$RunNum)){
  if(df$RunNum[i] %in% c(1,2,4)){
    df$PlateNum[i] = 1
  }
  if(df$RunNum[i] %in% c(3,5,6)){
    df$PlateNum[i] = 2
  }
  if(df$RunNum[i] %in% c(7,8,10,12,14)){
    df$PlateNum[i] = 3
  }
  if(df$RunNum[i] %in% c(9,11,13,15,16)){
    df$PlateNum[i] = 4
  }
}

#put all data in a column called lumData
for(i in 1:length(df$RunNum)){
  if(df$RunType[i] == "Luc"){
    df$lumData[i] = df$LucData[i] 
  }
  else{
    df$lumData[i] = df$RLData[i]
  }
}

############analysis#####
for (i in 1:16){
  lum.data <- subset(df, RunNum == i)
  plotRL <- raw_map(data = lum.data$lumData, 
                    well = lum.data$PlateLocation,
                    plate = 384) +
    scale_fill_viridis()+
    theme_dark() +
    ggtitle(paste("plate ", i))
  print(plotRL)
 
}

#run the linear analysis
df$row <- df$PlateLocation
for (i in c(10, 20,1:9,11:19,21:24)){
  df$row<-gsub(paste0(i),"",df$row)
}
df$rowNum <- rep(1:16, 24)#make the rows into numbers for the linear model
#model all of them
for (i in 1:16){
  #model plate 1 RL
  print(paste("plate", i))
  print(paste("RL, plate", i))
  RL.data <- subset(df, RunNum == i)
  df.model.RL <- lm(lumData~rowNum, data = RL.data) #define the linear model
  df.eqrl = function(x){summary(df.model.RL)$coefficients[2,1]*x + (summary(df.model.RL)$coefficients[1,1])}
  print(summary(df.model.RL))
  
  #plot it
  #RL.data.melt <- melt(RL.data, c("PlateLocation", "PlateNum", "row", "rowNum", "CloneKey", "Odor1K", "Odor1C", "Notes")) 
  #graph!
  print(
    ggplot(RL.data , aes(x = row, y = lumData)) +
      geom_point() +
      geom_smooth(method = lm) +
      stat_function(fun = df.eqrl, colour = "#F8766D") +
      annotate("text", x=2, y = max(RL.data$lumData), label = paste("p=",summary(df.model.RL)$coefficients[2,4]), colour = "#F8766D")+
      annotate("text", x=5, y = max(RL.data$lumData), label = paste("R^2=",summary(df.model.RL)$r.squared), colour = "#F8766D")+
      ggtitle(paste("plate", i))
  )
}
##
ggplot(data = subset(df, PlateNum == 2 & (RunType %in% c("RL1", "RL2"))) , aes(x = row, y = lumData, color = RunType)) +
  geom_point()



ggplot(data = subset(df, PlateNum == 4 & (RunType %in% c("RL1", "RL2", "RL3", "RL4"))) , aes(x = row, y = lumData, color = RunType)) +
  geom_point()

#make linear models and then run ANOVA (not aov) command
#plate 3
RL1.3 <- lm(lumData~rowNum, data = subset(df, PlateNum == 3 & RunType == "RL1")) #define the linear model
funRL1.3 = function(x){summary(RL1.3)$coefficients[2,1]*x + (summary(RL1.3)$coefficients[1,1])}
RL2.3 <- lm(lumData~rowNum, data = subset(df, PlateNum == 3 & RunType == "RL2"))
funRL2.3 = function(x){summary(RL2.3)$coefficients[2,1]*x + (summary(RL2.3)$coefficients[1,1])}
RL3.3 <- lm(lumData~rowNum, data = subset(df, PlateNum == 3 & RunType == "RL3"))
funRL3.3 = function(x){summary(RL3.3)$coefficients[2,1]*x + (summary(RL3.3)$coefficients[1,1])}
RL4.3 <- lm(lumData~rowNum, data = subset(df, PlateNum == 3 & RunType == "RL4"))
funRL4.3 = function(x){summary(RL4.3)$coefficients[2,1]*x + (summary(RL4.3)$coefficients[1,1])}

#pdf("170216_plate3.pdf")
ggplot(data = subset(df, PlateNum == 3 & (RunType %in% c("RL1", "RL2", "RL3", "RL4"))) , aes(x = row, y = lumData, color = RunType)) +
  geom_point() +
  stat_function(fun = funRL1.3, color = "#F8766D")+
  annotate("text", x=8, y = 7000, label = paste("p=",summary(RL1.3)$coefficients[2,4]), colour = "#F8766D")+
  annotate("text", x=15, y = 7000, label = paste("R^2=",summary(RL1.3)$r.squared), colour = "#F8766D")+
  stat_function(fun = funRL2.3, color = "#7CAE00")+
  annotate("text", x=8, y = 6800, label = paste("p=",summary(RL2.3)$coefficients[2,4]), colour = "#7CAE00")+
  annotate("text", x=15, y = 6800, label = paste("R^2=",summary(RL2.3)$r.squared), colour = "#7CAE00")+
  stat_function(fun = funRL3.3, color = "#00BFC4")+
  annotate("text", x=8, y = 6600, label = paste("p=",summary(RL3.3)$coefficients[2,4]), colour = "#00BFC4")+
  annotate("text", x=15, y = 6600, label = paste("R^2=",summary(RL3.3)$r.squared), colour = "#00BFC4")+
  stat_function(fun = funRL4.3, color = "#C77CFF") +
  annotate("text", x=8, y = 6400, label = paste("p=",summary(RL4.3)$coefficients[2,4]), colour = "#C77CFF")+
  annotate("text", x=15, y = 6400, label = paste("R^2=",summary(RL4.3)$r.squared), colour = "#C77CFF")
#dev.off()
  
#plate 4
RL1.3 <- lm(lumData~rowNum, data = subset(df, PlateNum == 4 & RunType == "RL1")) #define the linear model
funRL1.3 = function(x){summary(RL1.3)$coefficients[2,1]*x + (summary(RL1.3)$coefficients[1,1])}
RL2.3 <- lm(lumData~rowNum, data = subset(df, PlateNum == 4 & RunType == "RL2"))
funRL2.3 = function(x){summary(RL2.3)$coefficients[2,1]*x + (summary(RL2.3)$coefficients[1,1])}
RL3.3 <- lm(lumData~rowNum, data = subset(df, PlateNum == 4 & RunType == "RL3"))
funRL3.3 = function(x){summary(RL3.3)$coefficients[2,1]*x + (summary(RL3.3)$coefficients[1,1])}
RL4.3 <- lm(lumData~rowNum, data = subset(df, PlateNum == 4 & RunType == "RL4"))
funRL4.3 = function(x){summary(RL4.3)$coefficients[2,1]*x + (summary(RL4.3)$coefficients[1,1])}

#pdf("170216_plate4.pdf")
ggplot(data = subset(df, PlateNum == 4 & (RunType %in% c("RL1", "RL2", "RL3", "RL4"))) , aes(x = row, y = lumData, color = RunType)) +
  geom_point() +
  stat_function(fun = funRL1.3, color = "#F8766D")+
  annotate("text", x=8, y = 6300, label = paste("p=",summary(RL1.3)$coefficients[2,4]), colour = "#F8766D")+
  annotate("text", x=15, y = 6300, label = paste("R^2=",summary(RL1.3)$r.squared), colour = "#F8766D")+
  stat_function(fun = funRL2.3, color = "#7CAE00")+
  annotate("text", x=8, y = 6100, label = paste("p=",summary(RL2.3)$coefficients[2,4]), colour = "#7CAE00")+
  annotate("text", x=15, y = 6100, label = paste("R^2=",summary(RL2.3)$r.squared), colour = "#7CAE00")+
  stat_function(fun = funRL3.3, color = "#00BFC4")+
  annotate("text", x=8, y = 5900, label = paste("p=",summary(RL3.3)$coefficients[2,4]), colour = "#00BFC4")+
  annotate("text", x=15, y = 5900, label = paste("R^2=",summary(RL3.3)$r.squared), colour = "#00BFC4")+
  stat_function(fun = funRL4.3, color = "#C77CFF") +
  annotate("text", x=8, y = 5700, label = paste("p=",summary(RL4.3)$coefficients[2,4]), colour = "#C77CFF")+
  annotate("text", x=15, y = 5700, label = paste("R^2=",summary(RL4.3)$r.squared), colour = "#C77CFF")
#dev.off()

  