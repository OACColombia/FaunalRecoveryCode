##~~~~~~~~~~~~~~ Map of the studies included ~~~~~~~~~~~~~~#
#Clear memory
rm(list=ls())
#working in my laptop, set the directory
setwd("C:/Users/Flaco/Dropbox/Tesis MSc/Capitulo 1 - Review/Articulos base/Results/GEB/Submitted/Resubmission")
#Call the file MainData.csv
mydata <- read.csv("C:/Users/Flaco/Dropbox/Tesis MSc/Capitulo 1 - Review/Articulos base/Results/GEB/Submitted/Resubmission/MainData.csv")
dataAll=mydata
#Map
#First install the packages "maps", "PBSmapping", and "data.table"
install.packages("maps")
install.packages("PBSmapping")
install.packages("data.table")
#Then call the packages, and additional "ggplot2"
library("PBSmapping")
library("ggplot2")
library("maps")
library("data.table")
#set worldmap data
worldmap = map_data("world")
setnames(worldmap, c("X","Y","PID","POS","region","subregion"))
worldmap = clipPolys(worldmap, xlim=xlim,ylim=ylim, keepExtra=TRUE)
#graph
data=subset(dataAll,Forest.cov!="OGF")
data$Taxa <- factor(data$Taxa,levels=c('Amphibians','Reptiles','Birds','Mammals'))
data$Forest.cov<-factor(data$Forest.cov,levels=c('ES','YSF','MSF','OSF'))
Map <- ggplot(data,aes(x=Mean_Longitude,y=Mean_Latitude))+
  borders("world", fill="lightgray")+
  geom_hline(yintercept=0, linetype="dashed", color="gray")+
  geom_point(data=data,shape=21, stroke=0.5, colour="white", fill="black", size = 8)+
  theme_bw()+
  labs(y="", x="")+
  scale_x_continuous(limits = c(-130,180), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-60,35), expand = c(0, 0)) +
  annotate("text",x=-100, y=-10, label="Neotropic\n n = 46\n(M = 39, D = 6, S = 1)", size=8)+
  annotate("text",x=-5, y=-7, label="Afrotropic\n n = 8\n(M = 8)", size=8)+
  annotate("text",x=80, y=-3, label="Indo Malay\n n = 10\n(M = 10)", size=8)+
  annotate("text",x=150, y=8, label="Australasia\n n = 8\n(M = 4, S = 4)", size=8)+
  theme(legend.position = c(0.065, 0.5), 
        legend.title = element_blank( ),legend.background = element_rect( ), 
        legend.text = element_text(colour=rgb(0,0,0), size = 15),
        legend.key.size=unit(.8, "cm"))+
  theme(axis.text.x = element_text(hjust = 1), 
        axis.text=element_text(size=0),axis.title=element_text(size=0,face="bold"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())
Map#the map of the studies
