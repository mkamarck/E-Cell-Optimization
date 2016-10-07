setwd("S:/mainland/Projects/TAARs/Antag screen") #for PC
setwd("/Volumes/mainland/Projects/TAARs/Antag screen") #for MAc

library(plyr)
library(reshape2)
library(ggplot2)

x = 817 #agonist #
y = "TMA" #agonist name
z = "hTAAR5" #receptor

Screen<-read.table("ScreenData140716_2.csv",header=T,sep=",") ##read in data
LowRL<-subset(Screen, Screen$RLData < 1000) #how many have low RL?
Screen<-subset(Screen, Screen$RLData > 1000) #subset out low RL values
Screen$NormLuc=Screen$Luc/Screen$RL #normalize by RL
Screen<-Screen[,c("CloneKey", "Odor1InstockKey", "Odor2InstockKey", "NormLuc")] #columns we need

##Antag List
Antag<-subset(Screen, Odor1InstockKey == x) ##subset fof agonist
Agx2<-subset(Antag, Odor2InstockKey == x) ##Agonist + Agonist in screen
Agx1<-subset(Antag, Odor2InstockKey == 811) ##Agonist only
Antag["PerAg"]<-Antag$NormLuc/(mean(Agx1$NormLuc))*100 ##% of avg of Agonist only
Antag["PerInhib"]<-100-Antag$PerAg ##%inhib
Antag <- Antag[order(Antag$PerAg) , ]
#Agx1<-subset(Antag, Odor2InstockKey == 811)
#write.table(subset(Antag, select = c("Odor2InstockKey", "PerInhib")), file = "/Volumes/mainland/Projects/TAARs/Chemical Similarity/AntagOutput.txt", sep = "\t", row.names = F) #edit by MK to get stuff for chemically similar analysis

write.table(Antag, file = paste0(z, ".txt"), sep="\t", row.names = F)

Antag <- subset(Antag, Odor2InstockKey != 917)
temp<-subset(Antag, Odor2InstockKey == 811)

theme_set(theme_gray(base_size = 18))
ggplot(Antag,aes(x=Odor2InstockKey,y=PerInhib)) + geom_point(colour="black", size=1, shape=20) +
  geom_hline(y=85,color="red") +
  #geom_point(data = subset(Antag, Odor2InstockKey == x), aes(x = Odor2InstockKey, y = PerInhib), colour = "#9900CC") +
  #geom_text(data = subset(Antag, Odor2InstockKey == x), aes(x= Odor2InstockKey, y = PerInhib, label = paste(y, "+", y, sep=" ")), colour = "#9900CC", size = 4, vjust=2) + 
  geom_point(data = subset(Antag, Odor1InstockKey == x & Odor2InstockKey == 811), aes(x = Odor2InstockKey, y = PerInhib), colour = "blue") + 
  geom_text(data = subset(Antag, Odor1InstockKey == x & Odor2InstockKey == 811)[1,], aes(x= Odor2InstockKey, y = PerInhib, label = "phenylethylamine"), colour = "blue", size = 4, vjust=-1) +
  #geom_point(data = subset(Antag, Odor1InstockKey == x & Odor2InstockKey == 265), aes(x = Odor2InstockKey, y = PerInhib), colour = "orange") + 
  #geom_text(data = subset(Antag, Odor1InstockKey == x & Odor2InstockKey == 265)[1,], aes(x= Odor2InstockKey, y = PerInhib, label = paste(y, "+", "DMSO", sep=" ")), colour = "orange", size = 4, vjust=-1) + 
  ylab("Inhibition(%)") +
  xlab("Antagonist") + 
  theme(panel.background = element_rect(fill="grey95"))

ggplot(Antag, aes(x = NormLuc)) + 
  geom_freqpoly(binwidth=0.005) + 
  scale_y_log10()+scale_x_continuous(breaks=seq(-3,10,1)) + 
  xlab("normalized luciferase value") +
  ylab("log(count)")


##forskolin
For<-subset(Screen, Odor1InstockKey == "77")
ForOnly<-subset(For, Odor2InstockKey == "77")
For["PerFor"]<-(For$NormLuc/mean(ForOnly$NormLuc))*100
For["PerInhib"]<-100-For$PerFor
For <- For[order(For$PerFor) , ]
write.table(For, file = "forAntag.txt", sep="\t", row.names = F)

ggplot(Antag, aes(x = NormLuc)) + 
  geom_freqpoly(binwidth=0.001) + 
  scale_y_log10()+scale_x_continuous(breaks=seq(-3,10,1)) + 
  xlab("normalized luciferase value") +
  ylab("log(count)")




#AgOnly["PerAg"]<-AgOnly$NormLuc/2.32*100
#AgOnly["PerInhib"]<-100-AgOnly$PerAg
#AgOnly <- AgOnly[order(AgOnly$PerAg) , ]

#geom_text(data = AgOnly, aes(x=487, y = -18.72293186, label = "phenethylamine"), colour = "blue", size = 4, vjust=2) +
#geom_point(data=inhibAND,aes(x=id,y=inhib_AND),colour="blue", size=1, shape=20, label=inhibAND$odor2name)+
#geom_text(data=inhibAND,aes(x=id,y=inhib_AND,label=odor2name),colour="blue",size=3,hjust=1, vjust=-1, face="bold")+
#theme_bw() +
#theme(axis.text.x = element_blank(),
#panel.grid.major = element_blank(),
#panel.grid.minor = element_blank(),
#axis.ticks.x = element_blank()) +