#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Fit single model, age as single predictor for faunal recovery
#Clear memory
rm(list=ls())
setwd("c:/Users/Flaco/Dropbox/Tesis MSc/Capitulo 1 - Review/Articulos base/Results/GEB/Submitted/Resubmission")
#Call the file MainData.csv
mydata <- read.csv("C:/Users/Flaco/Dropbox/Tesis MSc/Capitulo 1 - Review/Articulos base/Results/GEB/Submitted/Resubmission/MainData.csv")
dataAll=mydata
data=subset(dataAll,Forest.cov!="OGF")
summary(data)
attach(data)
head(data)
nrow(data)#199 observations
#SS
data <- data[!apply(data[,c("Sim_S_toOGF", 
                            "Age", 
                            "Taxa", 
                            "Biome")], 1, anyNA),]
nrow(data)#170, I loss rows of studies without list
attach(data)
#~~~Simple model Amphibians
A<-subset(data, Taxa=="Amphibians")
nrow(A)#47 observations
length(unique(A$Nstudy))#14 studies
SSA<-lm(A$Sim_S_toOGF~1+A$Age)
summary(SSA)
A <- data.frame(A, predict(SSA, A, interval="confidence"))
head(A)#fit = model; lwr=lower confidence; upr=upper confidence
A <- data.frame(A, predict(SSA, A, interval="predict"))
head(A)#fit.1 = model; lwr.1=lower prediction; upr.1=upper prediction
upPIA<-lm(A$upr.1~1+A$Age)#upper interval line
lwPIA<-lm(A$lwr.1~1+A$Age)#lower interval line
#~~~~~Reptiles
R<-subset(data, Taxa=="Reptiles")
nrow(R)#37 observations
length(unique(R$Nstudy))#12 studies
SSR<-lm(R$Sim_S_toOGF~1+R$Age)
summary(SSR)
R <- data.frame(R, predict(SSR, R, interval="confidence"))
head(R)#fit = model; lwr=lower confidence; upr=upper confidence
R <- data.frame(R, predict(SSR, R, interval="predict"))
head(R)#fit.1 = model; lwr.1=lower prediction; upr.1=upper prediction
upPIR<-lm(R$upr.1~1+R$Age)#upper interval line
lwPIR<-lm(R$lwr.1~1+R$Age)#lower interval line
#~~~~~Birds
B<-subset(data, Taxa=="Birds")
nrow(B)#43 observations
length(unique(B$Nstudy))#21 studies
SSB<-lm(B$Sim_S_toOGF~1+B$Age)
summary(SSB)
B <- data.frame(B, predict(SSB, B, interval="confidence"))
head(R)#fit = model; lwr=lower confidence; upr=upper confidence
B <- data.frame(B, predict(SSB, B, interval="predict"))
head(B)#fit.1 = model; lwr.1=lower prediction; upr.1=upper prediction
upPIB<-lm(B$upr.1~1+B$Age)#upper interval line
lwPIB<-lm(B$lwr.1~1+B$Age)#lower interval line
#~~~~~~Mammals
M<-subset(data, Taxa=="Mammals")
nrow(M)#43 observations
length(unique(M$Nstudy))#18 studies
SSM<-lm(M$Sim_S_toOGF~1+M$Age)
summary(SSM)
M <- data.frame(M, predict(SSM, M, interval="confidence"))
head(M)#fit = model; lwr=lower confidence; upr=upper confidence
M <- data.frame(M, predict(SSM, M, interval="predict"))
head(M)#fit.1 = model; lwr.1=lower prediction; upr.1=upper prediction
upPIM<-lm(M$upr.1~1+M$Age)#upper interval line
lwPIM<-lm(M$lwr.1~1+M$Age)#lower interval line

#~~~Figure fit null model Sorensen Similarity and prediction time of recovery
library(ggplot2)
data$Taxa <- factor(data$Taxa,levels=c('Birds', 'Mammals','Reptiles','Amphibians'))

##changing the graph
#1~ points
p<-ggplot(data,aes(x=data$Age,y=data$Sim_S_toOGF, linetype=Taxa,
                   group=Taxa, fill=Taxa, colour=Taxa,
                   shape=Nbiome))+
  geom_point(data=A, aes(A$Age, A$Sim_S_toOGF, fill=Taxa), 
             size = 5, stroke = 1, alpha=1, colour="black", fill="gray")+
  geom_point(data=R, aes(R$Age, R$Sim_S_toOGF, fill=Taxa), 
             size = 5, stroke = 1, alpha=1, colour="black", fill="gray")+
  geom_point(data=B, aes(B$Age, B$Sim_S_toOGF, fill=Taxa), 
             size = 5, stroke = 1, alpha=1, colour="black", fill="gray")+
  geom_point(data=M, aes(M$Age, M$Sim_S_toOGF, fill=Taxa), 
             size = 5, stroke = 1, alpha=1, colour="black", fill="gray")
p
#2~some aesthetics
paes<-p+
  theme_bw()+
  ylab("Similarity to a reference forest \n (Sørensen index, SS)\n")+  
  xlab("\n Age (years after abandonment)")+
  scale_shape_manual(values=c(21,22,23,24,25,8))+
  scale_x_continuous(breaks = seq(0,250,by=50),
                     labels = c(seq(0,249,by=50), expression("250")),
                     limits = c(0,250), expand = c(0.02,0))+ 
  scale_y_continuous(breaks = seq(0,1.2, by=0.5),
                     labels = seq(0,1.2, by=0.5),
                     limits = c(0,1.5),expand = c(0.02,0))+
  scale_colour_manual(name="Taxa", 
                      values = c("Amphibians" = 'gray', "Reptiles"='gray', 
                                 "Birds"='gray', "Mammals"='gray'))+
  scale_fill_manual(name="Taxa",
                    values = c("Amphibians" = 'gray', "Reptiles"='gray', 
                               "Birds"='gray', "Mammals"='gray'))+
  theme(axis.text=element_text(size=22),axis.title=element_text(size=25),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.position="none")+
  facet_wrap(~Taxa, ncol=1, nrow=4, scales = "free_x")+ 
  theme(plot.margin = unit(c(1,1,1.5,1.2), "cm"))+
  theme(panel.margin = unit(1, "lines"))+
  theme(strip.background = element_blank())+
  theme(strip.text = element_text(size=20))
paes
#~3 fit the linear model
lpi<-paes+
  geom_hline(yintercept=0.9, linetype="dashed", colour="black", size=1.5)+
  geom_segment(data=A, aes(x=0, xend=110, y=(SSA$coefficients[1]+SSA$coefficients[2]),
                           yend=(SSA$coefficients[1]+(SSA$coefficients[2]*110))),
               colour="black", size=1,linetype="solid")+
  geom_segment(data=A, aes(x=0, xend=229, y=(lwPIA$coefficients[1]+lwPIA$coefficients[2]),
                           yend=(lwPIA$coefficients[1]+(lwPIA$coefficients[2]*229))),
               colour="gray", size=0.5,linetype="solid")+
  geom_segment(data=A, aes(x=0, xend=4, y=(upPIA$coefficients[1]+upPIA$coefficients[2]),
                           yend=(upPIA$coefficients[1]+(upPIA$coefficients[2]*4))),
               colour="gray", size=0.5,linetype="solid")+
  geom_segment(data=R, aes(x=0, xend=93, y=(SSR$coefficients[1]+SSR$coefficients[2]),
                           yend=(SSR$coefficients[1]+(SSR$coefficients[2]*93))),
               colour="black", size=1,linetype="solid")+
  geom_segment(data=R, aes(x=0, xend=249, y=(lwPIR$coefficients[1]+lwPIR$coefficients[2]),
                           yend=(lwPIR$coefficients[1]+(lwPIR$coefficients[2]*249))),
               colour="gray", size=0.5,linetype="solid")+
  geom_segment(data=R, aes(x=0, xend=0, y=(upPIR$coefficients[1]+upPIR$coefficients[2]),
                           yend=(upPIR$coefficients[1]+(upPIR$coefficients[2]*0))),
               colour="gray", size=0.5,linetype="solid")+
  geom_segment(data=B, aes(x=0, xend=36, y=(SSB$coefficients[1]+SSB$coefficients[2]),
                           yend=(SSB$coefficients[1]+(SSB$coefficients[2]*36))),
               colour="black", size=1,linetype="solid")+
  geom_segment(data=B, aes(x=0, xend=58, y=(lwPIB$coefficients[1]+lwPIB$coefficients[2]),
                           yend=(lwPIB$coefficients[1]+(lwPIB$coefficients[2]*58))),
               colour="gray", size=0.5,linetype="solid")+
  geom_segment(data=B, aes(x=0, xend=15, y=(upPIB$coefficients[1]+upPIB$coefficients[2]),
                           yend=(upPIB$coefficients[1]+(upPIB$coefficients[2]*15))),
               colour="gray", size=0.5,linetype="solid")+
  geom_segment(data=M, aes(x=0, xend=70, y=(SSM$coefficients[1]+SSM$coefficients[2]),
                           yend=(SSM$coefficients[1]+(SSM$coefficients[2]*70))),
               colour="black", size=1,linetype="solid")+
  geom_segment(data=M, aes(x=0, xend=150, y=(lwPIM$coefficients[1]+lwPIM$coefficients[2]),
                           yend=(lwPIM$coefficients[1]+(lwPIM$coefficients[2]*150))),
               colour="gray", size=0.5,linetype="solid")+
  geom_segment(data=M, aes(x=0, xend=5, y=(upPIM$coefficients[1]+upPIM$coefficients[2]),
                           yend=(upPIM$coefficients[1]+(upPIM$coefficients[2]*5))),
               colour="gray", size=0.5,linetype="solid")
lpi
#prediction interval or the time of recovery
rTime<-lpi+ 
  geom_segment(data=A, aes(x=4, xend=229, y=1.3, yend=1.3), 
               linetype="solid", colour="gray", size=1)+ 
  geom_segment(data=R, aes(x=0, xend=249, y=1.3, yend=1.3), 
               linetype="solid", colour="gray", size=1)+
  geom_segment(data=B, aes(x=15, xend=58, y=1.3, yend=1.3), 
               linetype="solid", colour="gray", size=1)+
  geom_segment(data=M, aes(x=5, xend=150, y=1.3, yend=1.3), 
               linetype="solid", colour="gray", size=1)+
  geom_point(data=A, aes(x=110, xend=110, y=1.3, yend=1.3), 
             colour="black", fill="black",shape=21, size=6)+
  geom_point(data=R, aes(x=93, xend=93, y=1.3, yend=1.3), 
             colour="black", fill="black",shape=21, size=6)+
  geom_point(data=B, aes(x=36, xend=36, y=1.3, yend=1.3), 
             colour="black", fill="black",shape=21, size=6)+
  geom_point(data=M, aes(x=70, xend=70, y=1.3, yend=1.3), 
             colour="black", fill="black",shape=21, size=6)
rTime#Figure of the recovery time
