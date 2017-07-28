#~~~~~~~~~~~~   Forest succession by age  ~~~~~#
#Clear memory
rm(list=ls())
#working in my laptop, set the directory
setwd("C:/Users/Flaco/Dropbox/Vertebrate Recovery - Nature/data")
#Call the file MainData.csv
mydata <- read.csv("C:/Users/Flaco/Dropbox/Vertebrate Recovery - Nature/data/MainData.csv")
data=mydata
names(data)
attach(data)
#Kruskal-Wallis and Bonferroni multi comparison of medians, age~succession(Habitat in table)
#Working with ggplot2
library(plyr)
library(Rmisc)
library(ggplot2)
kruskal.test(Mean.Age~SuccStage)#there are differences
pairwise.t.test(Mean.Age,SuccStage, p.adj="bonf")#The three categories are differents
data$SuccStage <- factor(data$SuccStage,levels=c('ES','YSF','MSF','OSF'))
#IQR for each succession category
ES<-subset(data,SuccStage=="ES")
summary(ES$Mean.Age)
YSF<-subset(data,SuccStage=="YSF")
summary(YSF$Mean.Age)
MSF<-subset(data,SuccStage=="MSF")
summary(MSF$Mean.Age)
OSF<-subset(data,SuccStage=="OSF")
summary(OSF$Mean.Age)
nrow(ES)
nrow(YSF)
nrow(MSF)
nrow(OSF)
a=ggplot(data = data, aes(x=SuccStage, y=Mean.Age))+
  geom_boxplot()+theme_bw()+
  annotate("text", x=1, y=15, label="0 (0,0)", size=5)+
  annotate("text", x=2, y=25, label="7 (5,13)", size=5)+
  annotate("text", x=3, y=50, label="22 (15,30)", size=5)+
  annotate("text", x=4, y=105, label="41 (27,65)", size=5)+
  labs(y="Years after abandonment\n", x="\n Successional stages")+  
  scale_x_discrete(labels = c("\n Early \n succession \n (n=57)", 
                              "\n Young \n secondary forest \n (n=107)",
                              "\n Mid-successional \n secondary forest \n (n=59)",
                              "\n Old \n secondary forest \n (n=12)"))+
  theme(axis.text=element_text(size=22), axis.title=element_text(size=25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
a  #Figure age by successional stages