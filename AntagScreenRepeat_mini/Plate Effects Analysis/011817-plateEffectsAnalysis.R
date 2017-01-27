#Analyzing the plate that was run with hTAAR5 and TMA (for whole plate) for plate effects
library(platetools)
library(ggplot2)
library(viridis)

#import data
df <- read.csv("011817-PlateEffectsAnalysis.csv")
#make map
raw_map(data = df$RLData, 
        well = df$PlateLocation, 
        plate = 384) +
  scale_fill_viridis()
  
raw_map(data = df$LucData, 
        well = df$PlateLocation, 
        plate = 384) +
  scale_fill_viridis()

df$norm <- df$LucData/df$RLData

raw_map(data = df$norm, 
        well = df$PlateLocation, 
        plate = 384) +
  scale_fill_viridis()