#analyze the data from 1-11-17
#This mini screen made a new odor dilution block and added lots of extra TMA alone variables
#import
library(tidyverse)
library(drc)
#import data
data <- read.csv("/Volumes/mainland/Projects/TAARs/E\ cell\ Optimization/AntagScreenRepeat_mini/170111_Results.csv", stringsAsFactors = FALSE)
#get rid of unnecessary columns
#gives all data
data <- subset(data, select = c("TransfectionDate", "PlateNum", "CloneKey", "Odor1K", "Odor2K","Odor1C", "Odor2C", "PlateLocation", "LucData", "RLData", "Notes","LabNoteBook", "TransfectionType", "CellType"))
#cuts down more columns to the useful ones for this analysis
data.sub <- subset(data, select = c("PlateNum", "CloneKey", "Odor1K", "Odor2K","Odor1C", "Odor2C", "PlateLocation", "LucData", "RLData", "CellType", "TransfectionType"))

#calculate normalized luc response
data.sub$norm <- data.sub$LucData/data.sub$RLData
#combine per cell type
eCells <- subset(data.sub, CellType == "E Cells")
H3A <- subset(data.sub, CellType == "H3A")

#####ANALYZE THE ANTAG RESPONSE: E Cells####
#E Cells
#plate 2
x = 817 #agonist #
y = "TMA" #agonist name
z = "hTAAR5" #receptor

Screen<-subset(eCells, Odor1K == x ) # get subset for the plate we are interested in
LowRL<-subset(Screen, Screen$RLData < 1000) #how many have low RL? (for H3A cells)
Screen<-subset(Screen, Screen$RLData > 1000) #subset out low RL values

#for E Cells, there are no RL values that low
#find standard deviation of the RL values for E Cells
meanRL <- mean(Screen$RLData)
sdRL <- sd(Screen$RLData)
#set min and max for more than 2 sd's outside of the mean RL value
minRL <- meanRL - 2*sdRL
maxRL <- meanRL + 2*sdRL
LowRL <- subset(Screen, Screen$RLData < minRL)
#no low RL for this set of E Cells
HighRL <- subset(Screen, Screen$RLData > maxRL)
#there are some values for highRL, but we are going to keep them for now. 

#Antag<-subset(Screen, Odor1InstockKey == x) ##subset fof agonist
Agx2 <- subset(Screen, Odor2K == x) ##Agonist + Agonist in screen
Agx1 <- subset(Screen, Odor2K == 811) ##Agonist only
Screen$PerAg <- Screen$norm/(mean(Agx1$norm))*100 ##% of avg of Agonist only
Screen$PerInhib <- 100-Screen$PerAg ##%inhib
Antag <- Screen[order(Screen$PerAg) , ]
#add rankings to antag
Antag$rankInhib <- c(1:length(Antag$norm))
#pull out the "top hits"
OldTopHits.E <- Antag[Antag$Odor2K %in% c(678, 486, 1221, 527, 821, 827, 817),]

#####Graph E Cell results#####
ggplot(data = Antag, aes(x = rankInhib, y = PerAg)) +
  geom_point() + 
  geom_point(data = subset(Antag, Odor2K == 811), aes(colour = "TMA Alone")) +
  ylab("% Response of mean 'agonist alone' response") +
  xlab("Rank of % inhibition") +
  ggtitle("E Cells with hTAAR5 Antagonist Response") +
  theme(legend.position = c(.1,.9),
        legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))



#####make heatmaps of raw values - E Cells #####
#make a heatmap of the plate - specifically of the wells that I'm interested in
library(platetools)
library(ggthemes)
library(viridis)


#for aesthetics change the source code of plate tools
#trace(plt384, edit = TRUE)

#my attempt
#pdf("/Volumes/mainland/Projects/TAARs/E\ cell\ Optimization/AntagScreenRepeat_mini/Plate Effects Analysis/170111_ECell_TAAR5_Luc.pdf")
raw_map(data = Screen$LucData,
      well = Screen$PlateLocation, 
      plate =384) +
  theme_dark() +
  scale_fill_viridis() +
  ggtitle("Luc values of E-Cells and hTAAR5")
#dev.off()
#to try to figure out 
#1. look at other 3 plates
#2. look at luciferase by itself
#3. look at Rl by itself
#4. look at b-scores? try other things they suggested on the website
#5. figure out what the hell is going on in here!

# panel.spacing.x
# raw_map(data = eCells$norm,
#         well = eCells$PlateLocation, 
#         plate =384, 
#         plate_id = eCells$PlateNum) +
#   theme_dark() +
#   scale_fill_viridis() +
#   ggtitle("normalized luciferase values of E Cells")
#some error message about spacing that I don't understand, but this doesn't illuminate the differences as well becuase the normalized values are so different. 
#pdf("/Volumes/mainland/Projects/TAARs/E\ cell\ Optimization/AntagScreenRepeat_mini/Plate Effects Analysis/170111_ECell_Rho_Luc.pdf")
eCellsRho <- subset(eCells, Odor1K == 77)
raw_map(data = eCellsRho$LucData,
        well = eCellsRho$PlateLocation, 
        plate =384) +
  theme_dark() +
  scale_fill_viridis() +
  ggtitle("Luc values of Rho and E Cells")
#dev.off()
#also we can look at H3A cells


#############Look at H3A data####
x = 817
Screen<-subset(H3A, Odor1K == 817)
LowRL<-subset(Screen, Screen$RL < 1000) #how many have low RL?
#Screen<-subset(Screen, Screen$RL > 1000) #subset out low RL values

#Antag<-subset(Screen, Odor1InstockKey == x) ##subset fof agonist
Agx2 <- subset(Screen, Odor2K == x) ##Agonist + Agonist in screen
Agx1 <- subset(Screen, Odor2K == 811) ##Agonist only
Screen$PerAg <- Screen$norm/(mean(Agx1$norm))*100 ##% of avg of Agonist only
Screen$PerInhib <- 100-Screen$PerAg ##%inhib
Antag <- Screen[order(Screen$PerAg) , ]
#add rankings to antag
Antag$rank <- c(1:length(Antag$norm))
#pull out the "top hits"
OldTopHits.E <- Antag[Antag$Odor2K %in% c(678, 486, 1221, 527, 821, 827, 817),]

#AntagResponse(eCells, 1, 77)
# getAgResponse <- function(Screen, agKey){
#   subScreen <- subset
# }
#############look at plate maps H3A####
#pdf("/Volumes/mainland/Projects/TAARs/E\ cell\ Optimization/AntagScreenRepeat_mini/Plate Effects Analysis/170111_H3A_TAAR5_rl.pdf")
raw_map(data = Screen$RLData,
        well = Screen$PlateLocation, 
        plate =384) +
  theme_dark() +
  scale_fill_viridis() +
  ggtitle("RL values of hTAAR5 with H3A cells")
#dev.off()
#there is still a gradient - definitely a gradient, but its not quite as dramatic as with the E cells. 

#look at H3A with Rho
H3ARho <- subset(H3A, Odor1K == 77)
pdf("/Volumes/mainland/Projects/TAARs/E\ cell\ Optimization/AntagScreenRepeat_mini/Plate Effects Analysis/170111_H3A_Rho_RL.pdf")
raw_map(data = H3ARho$RLData,
        well = H3ARho$PlateLocation, 
        plate =384) +
  #geom_point(type= ,size=3, fill=) +
  theme_dark() +
  scale_fill_viridis() +
  ggtitle("RL values of Rho with H3A cells")
dev.off()
