#Analyzing the plate that was run with hTAAR5 and TMA (for whole plate) for plate effects
library(platetools)
library(ggplot2)
library(viridis)
library(plyr)
library(reshape2)
library(pipeR)

#####1-18-17#####
#import data
df <- read.csv("170207_plateEffects.csv") %>>% subset(select = c(PlateNum, CloneKey, Odor1K, Odor1C, PlateLocation, RLData, LucData, Notes))

#graph RL of all plates
for (i in 1:6){
  RL.data <- subset(df, PlateNum == i)
  plotRL <- raw_map(data = RL.data$RLData, 
          well = RL.data$PlateLocation,
          plate = 384) +
    scale_fill_viridis()+
    theme_dark() +
    ggtitle(paste("plate ", i))
  print(plotRL)
  
}
RL.data.5a <- subset(df, PlateNum == 6)
raw_map(data = RL.data.5a$LucData, 
        well = RL.data.5a$PlateLocation,
        plate = 384) +
  scale_fill_viridis()+
  theme_dark() +
  ggtitle("plate 5 second run through")

#graph luc data
for (i in 1:5){
  RL.data <- subset(df, PlateNum == i)
  plotRL <- raw_map(data = RL.data$LucData, 
                    well = RL.data$PlateLocation,
                    plate = 384) +
    scale_fill_viridis()+
    theme_dark() +
    ggtitle(paste("plate ", i))
  print(plotRL)
  
}

#####model the data for effects#####
df$row <- df$PlateLocation
for (i in c(10, 20,1:9,11:19,21:24)){
  df$row<-gsub(paste0(i),"",df$row)
}
df$rowNum <- rep(1:16, 24)#make the rows into numbers for the linear model
#model all of them
for (i in 1:6){
  #model plate 1 RL
  print(paste("plate", i))
  print(paste("RL, plate", i))
  RL.data <- subset(df, PlateNum == i)
  df.model.RL <- lm(RLData~rowNum, data = RL.data) #define the linear model
  df.eqrl = function(x){summary(df.model.RL)$coefficients[2,1]*x + (summary(df.model.RL)$coefficients[1,1])}
  print(summary(df.model.RL))
  print(paste("luc, plate", i))
  #put data in df
  #   p.r.data[i+1,4] <- summary(df.model.RL)$r.squared
  #   p.r.data[i+1,3] <- summary(df.model.RL)$coefficients[2,4]
  df.model <- lm(LucData~rowNum, data = RL.data) #define the linear model
  df.eqluc = function(x){summary(df.model)$coefficients[2,1]*x + (summary(df.model)$coefficients[1,1])}
  print(summary(df.model))
  #   p.r.data[i+1,6] <- summary(df.model)$r.squared
  #   p.r.data[i+1,5] <- summary(df.model)$coefficients[2,4]
  #plot it
  RL.data.melt <- melt(RL.data, c("PlateLocation", "PlateNum", "row", "rowNum", "CloneKey", "Odor1K", "Odor1C", "Notes")) 
  #graph!
  print(
    ggplot(data = subset(RL.data.melt, variable %in% c("RLData", "LucData")) , aes(x = row, y = value, colour = variable)) +
      geom_point() +
      geom_smooth(method = lm) +
      stat_function(fun = df.eqrl, colour = "#F8766D") +
      annotate("text", x=2, y = max(RL.data$LucData), label = paste("p=",summary(df.model.RL)$coefficients[2,4]), colour = "#F8766D")+
      annotate("text", x=5, y = max(RL.data$LucData), label = paste("R^2=",summary(df.model.RL)$r.squared), colour = "#F8766D")+
      stat_function(fun = df.eqluc, colour = "#00BFC4") +
      annotate("text", x=2, y = max(RL.data$LucData) + 200, label = paste("p=",summary(df.model)$coefficients[2,4]), colour = "#00BFC4")+
      annotate("text", x=5, y = max(RL.data$LucData) + 200, label = paste("R^2=",summary(df.model)$r.squared), colour = "#00BFC4")+
      ggtitle(paste("plate", i))
  )
}


#########old code to pull from#####
#make a column that just has the row information from the plate
df$row <- df$PlateLocation
for (i in c(10, 20,1:9,11:19,21:24)){
  df$row<-gsub(paste0(i),"",df$row)
}
df$rowNum <- rep(1:16, 24)#make the rows into numbers for the linear model
df.model <- lm(RLData~rowNum, data = df) #define the linear model
df.eq = function(x){summary(df.model)$coefficients[2,1]*x + (summary(df.model)$coefficients[1,1])}
#plot this with the model
ggplot(df, aes(x = row, y = RLData)) +
  geom_point() +
  #geom_smooth(method = lm) +
  stat_function(fun = df.eq, colour = "#00BFC4") 
#add stuff to p and r value data frame
#   p.r.data <-  as.data.frame(c("1-17-17", "1-26-17", "1-26-17", "1-26-17", "1-26-17"))
#   p.r.data$plateNum <- c(1,1,2,3,4)
#   names(p.r.data) <- c("plateDate", "plateNum")
#   p.r.data$RL_pValue <- rep(0,5)
#   p.r.data$RL_R <- rep(0,5)
#   p.r.data$Luc_pValue <- rep(0,5)
#   p.r.data$Luc_R <- rep(0,5)
#   p.r.data[1,4] <- summary(df.model)$r.squared
#   p.r.data[1,3] <- summary(df.model)$coefficients[2,4]
#Lets check it for luciferase too, to be sure that there isn't a gradient
#I should probably check for gradient in the other direction too...
df.model <- lm(LucData~rowNum, data = df) #define the linear model
df.eq = function(x){summary(df.model)$coefficients[2,1]*x + (summary(df.model)$coefficients[1,1])}
ggplot(df, aes(x = row, y = LucData)) +
  geom_point() +
  geom_smooth(method = lm) +
  stat_function(fun = df.eq, colour = "#00BFC4") 
#wow, luciferase does have a gradient - I think I saw this


#fill in this database
#   p.r.data[1,6] <- summary(df.model)$r.squared
#   p.r.data[1,5] <- summary(df.model)$coefficients[2,4]
#check for horizontal effects on the plates
df$column <- df$PlateLocation
columns <- list("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P")
for(i in 1:length(columns)){
  df$column <- gsub(paste0(columns[i]), "", df$column)
}
df$column <- as.numeric(df$column)

df.model <- lm(RLData~column, data = df) #define the linear model
df.eq = function(x){summary(df.model)$coefficients[2,1]*x + (summary(df.model)$coefficients[1,1])}
ggplot(df, aes(x = column, y = RLData)) +
  geom_point() +
  #geom_smooth(method = lm) +
  stat_function(fun = df.eq, colour = "#00BFC4") +
  annotate("text", x=10, y = 7000, label = paste("p=",summary(df.model)$coefficients[2,4]), colour = "#00BFC4") +
  ggtitle("RL across columns")

df.model <- lm(LucData~column, data = df) #define the linear model
df.eq = function(x){summary(df.model)$coefficients[2,1]*x + (summary(df.model)$coefficients[1,1])}
ggplot(df, aes(x = column, y = LucData)) +
  geom_point() +
  #geom_smooth(method = lm) +
  stat_function(fun = df.eq, colour = "#00BFC4") +
  annotate("text", x=10, y = 5000, label = paste("p=",summary(df.model)$coefficients[2,4]), colour = "#00BFC4") +
  ggtitle("Luc across columns")
#put luc data in database


#########1-26-17##### 
#import data from 1-26-17
plateTest <- read.csv("170127-PlateEffects.csv")

plateTest$norm <- plateTest$LucData/plateTest$RLData

RLmax <- mean(plateTest$RLData) + 2*sd(plateTest$RLData)
plateTest <- ddply(plateTest, .variables = c("PlateNum"), function(x) (subset(x, x$RLData<(mean(x$RLData) + 3*sd(plateTest$RLData)))))

#plot the maps
raw_map(data = plateTest[which(plateTest$PlateNum == 1),]$RLData, 
        well = plateTest[which(plateTest$PlateNum == 1),]$PlateLocation, 
        plate = 384) +
  scale_fill_viridis() +
  theme_dark()


raw_map(data = plateTest[which(plateTest$PlateNum ==2),]$RLData, 
        well = plateTest[which(plateTest$PlateNum == 2),]$PlateLocation, 
        plate = 384) +
  scale_fill_viridis() +
  theme_dark()



raw_map(data = plateTest[which(plateTest$PlateNum ==3),]$RLData, 
        well = plateTest[which(plateTest$PlateNum == 3),]$PlateLocation, 
        plate = 384) +
  scale_fill_viridis() +
  theme_dark()

raw_map(data = plateTest[which(plateTest$PlateNum ==4),]$RLData, 
        well = plateTest[which(plateTest$PlateNum == 4),]$PlateLocation, 
        plate = 384) +
  scale_fill_viridis()+
  theme_dark()

#now lets look at the models of RL
#plate 1
plateTest <- read.csv("170127-PlateEffects.csv")

plateTest$norm <- plateTest$LucData/plateTest$RLData
plateTest$row <- df$PlateLocation
for (i in c(10, 20,1:9,11:19,21:24)){
 plateTest$row<-gsub(paste0(i),"",plateTest$row)
}
plateTest$rowNum <- rep(1:16, 24*4)#make the rows into numbers for the linear model

for (i in 1:4){
  #model plate 1 RL
  print(paste("plate", i))
  print(paste("RL, plate", i))
  df.model.RL <- lm(RLData~rowNum, data = plateTest[which(plateTest$PlateNum ==i),]) #define the linear model
  df.eqrl = function(x){summary(df.model.RL)$coefficients[2,1]*x + (summary(df.model.RL)$coefficients[1,1])}
  print(summary(df.model.RL))
  print(paste("luc, plate", i))
  #put data in df
#   p.r.data[i+1,4] <- summary(df.model.RL)$r.squared
#   p.r.data[i+1,3] <- summary(df.model.RL)$coefficients[2,4]
  df.model <- lm(LucData~rowNum, data = plateTest[which(plateTest$PlateNum ==i),]) #define the linear model
  df.eqluc = function(x){summary(df.model)$coefficients[2,1]*x + (summary(df.model)$coefficients[1,1])}
  print(summary(df.model))
#   p.r.data[i+1,6] <- summary(df.model)$r.squared
#   p.r.data[i+1,5] <- summary(df.model)$coefficients[2,4]
  #plot it
  plateTest2 <- subset(plateTest, select = c("PlateLocation", "RLData", "LucData", "norm", "PlateNum", "row", "rowNum"))
  plateTest2 <- melt(plateTest2, c("PlateLocation", "PlateNum", "row", "rowNum") )
#   ggplot(plateTest[which(plateTest$PlateNum ==i),], aes(x = row, y = RLData)) +
#     geom_point() +
#     geom_smooth(method = lm) +
#     stat_function(fun = df.eq, colour = "#00BFC4") #+
#     #ggtitle(paste("plate", i, "RL"))
#   

  print(
    ggplot(data = subset(plateTest2, PlateNum == i & variable %in% c("RLData", "LucData")) , aes(x = row, y = value, colour = variable)) +
      geom_point() +
      geom_smooth(method = lm) +
      stat_function(fun = df.eqrl, colour = "#F8766D") +
      annotate("text", x=2, y = 10000, label = paste("p=",summary(df.model.RL)$coefficients[2,4]), colour = "#F8766D")+
      stat_function(fun = df.eqluc, colour = "#00BFC4") +
      annotate("text", x=2, y = 9800, label = paste("p=",summary(df.model)$coefficients[2,4]), colour = "#00BFC4")+
      ggtitle(paste("plate", i))
    )
  }
#export the chart of p and r values
#write.csv(p.r.data, "P and R values for verticle plate effects.csv")

#horizontal effects
plateTest$column <- plateTest$PlateLocation
columns <- list("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P")
for(i in 1:length(columns)){
  plateTest$column <- gsub(paste0(columns[i]), "", plateTest$column)
}
plateTest$column <- as.numeric(plateTest$column)
for (i in 1:4){
  #model plate 1 RL
  print(paste("plate", i))
  print(paste("RL, plate", i))
  df.model.RL <- lm(RLData~column, data = plateTest[which(plateTest$PlateNum ==i),]) #define the linear model
  df.eqrl = function(x){summary(df.model.RL)$coefficients[2,1]*x + (summary(df.model.RL)$coefficients[1,1])}
  print(summary(df.model.RL))
  print(paste("luc, plate", i))
  df.model <- lm(LucData~column, data = plateTest[which(plateTest$PlateNum ==i),]) #define the linear model
  df.eqluc = function(x){summary(df.model)$coefficients[2,1]*x + (summary(df.model)$coefficients[1,1])}
  print(summary(df.model))
  p.r.data[i+1,4] <- summary(df.model.RL)$r.squared
  p.r.data[i+1,3] <- summary(df.model.RL)$coefficients[2,4]
  p.r.data[i+1,6] <- summary(df.model)$r.squared
  p.r.data[i+1,5] <- summary(df.model)$coefficients[2,4]
  #plot it
  plateTest2 <- subset(plateTest, select = c("PlateLocation", "RLData", "LucData", "norm", "PlateNum", "row", "rowNum", "column"))
  plateTest2 <- melt(plateTest2, c("PlateLocation", "PlateNum", "row", "rowNum", "column") )

  print(
    ggplot(data = subset(plateTest2, PlateNum == i & variable %in% c("RLData", "LucData")) , aes(x = column, y = value, colour = variable)) +
      geom_point() +
      stat_function(fun = df.eqrl, colour = "#F8766D") +
      annotate("text", x=5, y = 10000, label = paste("p=",summary(df.model.RL)$coefficients[2,4]), colour = "#F8766D")+
      stat_function(fun = df.eqluc, colour = "#00BFC4") +
      annotate("text", x=5, y = 9800, label = paste("p=",summary(df.model)$coefficients[2,4]), colour = "#00BFC4")+
      ggtitle(paste("plate", i))
  )
}

#write.csv(p.r.data, "P and R values for horizontal plate effects.csv")
#####Is this comperable across plates?!?!?!#####
#if it is we can normalize it



