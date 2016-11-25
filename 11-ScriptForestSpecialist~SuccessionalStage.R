##~~~~~~~~~~~~~~ Forest specialist ~~~~~~~~~~~~~~#
#Clear memory
rm(list=ls())
#working in my laptop, set the directory
setwd("C:/Users/Flaco/Dropbox/Tesis MSc/Capitulo 1 - Review/Articulos base/Results/GEB/Submitted/Resubmission")
#Call the file MainData.csv
mydata <- read.csv("C:/Users/Flaco/Dropbox/Tesis MSc/Capitulo 1 - Review/Articulos base/Results/GEB/Submitted/Resubmission/MainData.csv")
dataAll=mydata
#packages
library("ggplot2")

AmpFS<-subset(dataAll,Taxa=="Amphibians")
hist(AmpFS$p.ForSpec)
shapiro.test(AmpFS$p.ForSpec)
KAmpFS<-kruskal.test(AmpFS$p.ForSpec~AmpFS$Forest.cov)
KAmpFS
pairwise.t.test(AmpFS$p.ForSpec,AmpFS$Forest.cov, p.adj="bonf")

RepFS<-subset(dataAll,Taxa=="Reptiles")
RepFS<-subset(RepFS,Forest.cov!="OSF")
hist(RepFS$p.ForSpec)
shapiro.test(RepFS$p.ForSpec)
KRepFS<-kruskal.test(RepFS$p.ForSpec~RepFS$Forest.cov)
KRepFS
pairwise.t.test(RepFS$p.ForSpec,RepFS$Forest.cov, p.adj="bonf")

BirFS<-subset(dataAll,Taxa=="Birds")
BirFS<-subset(BirFS,Forest.cov!="OSF")
hist(BirFS$p.ForSpec)
shapiro.test(BirFS$p.ForSpec)
KBirFS<-kruskal.test(BirFS$p.ForSpec~BirFS$Forest.cov)
KBirFS
pairwise.t.test(BirFS$p.ForSpec,BirFS$Forest.cov, p.adj="bonf")#ES is significant different (p<0.05)

MamFS<-subset(dataAll,Taxa=="Mammals")
MamFS<-subset(MamFS,Forest.cov!="ES")
MamFS<-subset(MamFS,Forest.cov!="OSF")
hist(MamFS$p.ForSpec)
shapiro.test(MamFS$p.ForSpec)
KMamFS<-kruskal.test(MamFS$p.ForSpec~MamFS$Forest.cov)
KMamFS
pairwise.t.test(MamFS$p.ForSpec,MamFS$Forest.cov, p.adj="bonf")

#Boxplot
#To graph order
dataAll$Taxa<-factor(dataAll$Taxa, levels=c('Amphibians','Reptiles','Birds','Mammals'))
dataAll$Forest.cov<-factor(dataAll$Forest.cov, levels=c('ES','YSF','MSF','OSF','OGF'))
#Boxplot
FS<-ggplot(dataAll, aes(x=Forest.cov, y=p.ForSpec, alpha=Forest.cov))+
  geom_boxplot(fill="grey15", colour="black",
               outlier.colour = "grey15", outlier.shape = 1)+
  facet_wrap(~Taxa, ncol=1, nrow=4)+
  ylab("Proportion of forest specialist species \n")+  
  xlab(" \n Successional stages")+
  scale_alpha_discrete(range=c(0, 0.75), name="", labels=c("ES", 
                                                           "YSF",
                                                           "MSF",
                                                           "OSF",
                                                           "Reference"))+
  scale_x_discrete(labels = c("ES", 
                              "YSF", 
                              "MSF",
                              "OSF", 
                              "Ref."))+
  scale_y_continuous(breaks = seq(0,1.2, by=0.25),
                     labels = seq(0,1.2, by=0.25),
                     limits = c(0,1.2))+
  theme_bw()+
  theme(axis.text=element_text(size=22), axis.title=element_text(size=25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(colour=rgb(0,0,0), size = 18),
        legend.key.size=unit(.8, "cm"))+ 
  theme(legend.background = element_rect(fill=rgb(1,1,1), size=0.5, linetype=2))+
  theme(legend.position="none")+
  theme(panel.margin = unit(2, "lines"))+
  theme(strip.background = element_blank())+
  theme(strip.text = element_text(size=20))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

FS+ guides(alpha=guide_legend(nrow=4,byrow=TRUE))+
  annotate("text",x=1,y=1.15,label=c("a (4)","","",""), size=6)+
  annotate("text",x=2,y=1.15,label=c("ab (9)","","",""), size=6)+
  annotate("text",x=3,y=1.15,label=c("ab (9)","","",""), size=6)+
  annotate("text",x=4,y=1.15,label=c("-","","",""), size=6)+
  annotate("text",x=5,y=1.15,label=c("b (12)","","",""), size=6)+
  annotate("text",x=1,y=1.15,label=c("","a (3)","",""), size=6)+
  annotate("text",x=2,y=1.15,label=c("","a (8)","",""), size=6)+
  annotate("text",x=3,y=1.15,label=c("","ab (8)","",""), size=6)+
  annotate("text",x=4,y=1.15,label=c("","(1)","",""), size=6)+
  annotate("text",x=5,y=1.15,label=c("","b (10)","",""), size=6)+
  annotate("text",x=1,y=1.15,label=c("","","a (3)",""), size=6)+
  annotate("text",x=2,y=1.15,label=c("","","b (15)",""), size=6)+
  annotate("text",x=3,y=1.15,label=c("","","b (6)",""), size=6)+
  annotate("text",x=4,y=1.15,label=c("","","-",""), size=6)+
  annotate("text",x=5,y=1.15,label=c("","","b (18)",""), size=6)+
  annotate("text",x=1,y=1.15,label=c("","","","(1)"), size=6)+
  annotate("text",x=2,y=1.15,label=c("","","","(14)"), size=6)+
  annotate("text",x=3,y=1.15,label=c("","","","(2)"), size=6)+
  annotate("text",x=4,y=1.15,label=c("","","","-"), size=6)+
  annotate("text",x=5,y=1.15,label=c("","","","(14)"), size=6)
