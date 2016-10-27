##############analyze the data from 9-30-16######
#this was an antagonist screen comparing results between H3A and E cells, and to make sure that the same top hits came up as the last time we ran the screen. 
#import
library(tidyverse)
library(drc)
#import data
plate1 <- read.csv("/Volumes/mainland-1/Projects/TAARs/E\ cell\ Optimization/AntagScreenRepeat_mini/160930_plate1.csv", stringsAsFactors = FALSE)[,1:11]
plate2 <- read.csv("/Volumes/mainland-1/Projects/TAARs/E\ cell\ Optimization/AntagScreenRepeat_mini/160930_plate2.csv", stringsAsFactors = FALSE)[,1:11]
plate3 <- read.csv("/Volumes/mainland-1/Projects/TAARs/E\ cell\ Optimization/AntagScreenRepeat_mini/160930_plate3.csv", stringsAsFactors = FALSE)[,1:11]
plate4 <- read.csv("/Volumes/mainland-1/Projects/TAARs/E\ cell\ Optimization/AntagScreenRepeat_mini/160930_plate4.csv", stringsAsFactors = FALSE)[,1:11]
#combine per cell type
eCells <- rbind(plate1, plate2)
H3A <- rbind(plate3, plate4)
#calculate the normalized response
eCells$norm <- eCells$Luc/eCells$RL
H3A$norm <- H3A$Luc/H3A$RL

#####ANALYZE THE ANTAG RESPONSE####
#E Cells
#plate 2
x = 817 #agonist #
y = "TMA" #agonist name
z = "hTAAR5" #receptor

Screen<-subset(eCells, Plate == 2 & agonistInStockKey == x & experimentType %in% c("antagScreen", "control"))
LowRL<-subset(Screen, Screen$RL < 1000) #how many have low RL?
Screen<-subset(Screen, Screen$RL > 1000) #subset out low RL values

Screen<-Screen[,c( "plateLocation", "antagInStockKey", "norm")] #columns we need

#Antag<-subset(Screen, Odor1InstockKey == x) ##subset fof agonist
Agx2<-subset(Screen, antagInStockKey == x) ##Agonist + Agonist in screen
Agx1<-subset(Screen, antagInStockKey == 811) ##Agonist only
Screen["PerAg"]<-Screen$norm/(mean(Agx1$norm))*100 ##% of avg of Agonist only
Screen["PerInhib"]<-100-Screen$PerAg ##%inhib
Antag <- Screen[order(Screen$PerAg) , ]
#add rankings to antag
Antag$rank <- c(1:length(Antag$norm))
#pull out the "top hits"
OldTopHits.E <- Antag[Antag$antagInStockKey %in% c(678, 486, 1221, 527, 821, 827, 817),]

#############Look at forskolin data ####
x = 77
Screen<-subset(eCells, Plate == 1 & agonistInStockKey == x & experimentType %in% c("antagScreen", "control"))
LowRL<-subset(Screen, Screen$RL < 1000) #how many have low RL?
Screen<-subset(Screen, Screen$RL > 1000) #subset out low RL values

Screen<-Screen[,c( "plateLocation", "antagInStockKey", "norm")] #columns we need

#Antag<-subset(Screen, Odor1InstockKey == x) ##subset fof agonist
Agx2<-subset(Screen, antagInStockKey == x) ##Agonist + Agonist in screen
Agx1<-subset(Screen, antagInStockKey == 811) ##Agonist only
Screen["PerAg"]<-Screen$norm/(mean(Agx1$norm))*100 ##% of avg of Agonist only
Screen["PerInhib"]<-100-Screen$PerAg ##%inhib
Antag <- Screen[order(Screen$PerAg) , ]
#add rankings to antag
Antag$rank <- c(1:length(Antag$norm))
OldTopHits.E.rho <- Antag[Antag$antagInStockKey %in% c(678, 486, 1221, 527, 821, 827, 817),]
#definitly something has gone wrong because 486 and 1221 are supposed to both be PEB

AntagResponse(eCells, 1, 77)
# getAgResponse <- function(Screen, agKey){
#   subScreen <- subset
# }

