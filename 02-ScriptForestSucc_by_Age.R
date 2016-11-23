#~~~~~~~~~~~~   Forest succession by age  ~~~~~#
#Clear memory
rm(list=ls())
#working in my laptop, set the directory
setwd("C:/Users/Flaco/Dropbox/Tesis MSc/Capitulo 1 - Review/Articulos base/Results/GEB/Submitted/Resubmission")
#Call the file MainData.csv
mydata <- read.csv("C:/Users/Flaco/Dropbox/Tesis MSc/Capitulo 1 - Review/Articulos base/Results/GEB/Submitted/Resubmission/MainData.csv")
dataAll=mydata
#Include just the "secondary forests" but old-growth forests (reference)
data=subset(dataAll,Forest.cov!="OGF")
names(data)
attach(data)
#Kruskal-Wallis and Bonferroni multi comparison of medians, age~succession(Habitat in table)
#Working with ggplot2
library(plyr)
library(Rmisc)
library(ggplot2)
kruskal.test(Age~Forest.cov)#there are differences
pairwise.t.test(Age,Forest.cov, p.adj="bonf")#The three categories are differents
data$Forest.cov <- factor(data$Forest.cov,levels=c('ES','YSF','MSF','OSF'))
#IQR for each succession category
ES<-subset(data,Forest.cov=="ES")
summary(ES$Age)
YSF<-subset(data,Forest.cov=="YSF")
summary(YSF$Age)
MSF<-subset(data,Forest.cov=="MSF")
summary(MSF$Age)
OSF<-subset(data,Forest.cov=="OSF")
summary(OSF$Age)
a=ggplot(data = data, aes(x=Forest.cov, y=Age))+
  geom_boxplot()+theme_bw()+
  annotate("text", x=1, y=10, label="0 (0,0)", size=5)+
  annotate("text", x=2, y=25, label="8 (5,13)", size=5)+
  annotate("text", x=3, y=50, label="25 (20,35)", size=5)+
  annotate("text", x=4, y=82, label="60 (55,65)", size=5)+
  labs(y="Years after abandonment\n", x="\n Successional stages")+  
  scale_x_discrete(labels = c("\n Early \n succession \n (n=27)", 
                              "\n Young \n secondary forest \n (n=119)",
                              "\n Mid-successional \n secondary forest \n (n=49)",
                              "\n Old \n secondary forest \n (n=4)"))+
  theme(axis.text=element_text(size=22), axis.title=element_text(size=25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
a  #Figure age by successional stages
