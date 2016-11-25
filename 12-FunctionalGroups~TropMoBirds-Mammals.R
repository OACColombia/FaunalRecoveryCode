##~~~~~~~~~~~~~~ Tropical moist forest birds and mammals' functional groups ~~~~~~~~~~~~~~#
#Clear memory
rm(list=ls())
#packages used
library(ggplot2)
library(reshape2)
#working in my laptop, set the directory
setwd("C:/Users/Flaco/Dropbox/Tesis MSc/Capitulo 1 - Review/Articulos base/Results/GEB/Submitted/Resubmission")
#Call the file MainData.csv
mydata <- read.csv("C:/Users/Flaco/Dropbox/Tesis MSc/Capitulo 1 - Review/Articulos base/Results/GEB/Submitted/Resubmission/MainData.csv")
dataAll=mydata
dataMo<-subset(dataAll,Nbiome=="1.TropMo")
dataNoAmp<-subset(dataMo,Taxa!="Amphibians")
dataNoHerp<-subset(dataNoAmp,Taxa!="Reptiles")
data<-subset(dataNoHerp,Moist_Forest_specialist!="NA")
names(data)
newdata<-data[c(6,25,50,52,54,56,58,60,62)]#working only with successional stage, Taxa and the proportion of each functional group
head(newdata)
fungui<-melt(newdata, id=c("Forest.cov","Taxa"))
head(fungui)
colnames(fungui)<-c("Suc.Stage","Taxa","Guild","Prop")
head(fungui)
fungui<-subset(fungui, Suc.Stage!="OSF")
#Birds
fgb<-subset(fungui, Taxa =="Birds")
#Insectivorous
Ifgb<-subset(fgb, Guild=="p.In")
hist(Ifgb$Prop)
shapiro.test(Ifgb$Prop)
KIfgb<-kruskal.test(Ifgb$Prop~Ifgb$Suc.Stage)
KIfgb
pairwise.t.test(Ifgb$Prop,Ifgb$Suc.Stage, p.adj="bonf")#No differences
#Seed dispersers
SDfgb<-subset(fgb, Guild=="p.SD")
hist(SDfgb$Prop)
shapiro.test(SDfgb$Prop)
KSDfgb<-kruskal.test(SDfgb$Prop~SDfgb$Suc.Stage)
KSDfgb
pairwise.t.test(SDfgb$Prop,SDfgb$Suc.Stage, p.adj="bonf")#No differences
#Granivorous
Grfgb<-subset(fgb, Guild=="p.Gr")
hist(Grfgb$Prop)
shapiro.test(Grfgb$Prop)
KGrfgb<-kruskal.test(Grfgb$Prop~Grfgb$Suc.Stage)
KGrfgb
pairwise.t.test(Grfgb$Prop,Grfgb$Suc.Stage, p.adj="bonf")#No differences
#Predators
Prfgb<-subset(fgb, Guild=="p.Pr")
hist(Prfgb$Prop)
shapiro.test(Prfgb$Prop)
KPrfgb<-kruskal.test(Prfgb$Prop~Prfgb$Suc.Stage)
KPrfgb
pairwise.t.test(Prfgb$Prop,Prfgb$Suc.Stage, p.adj="bonf")#ES is really different!
#Pollinators
Pofgb<-subset(fgb, Guild=="p.Po")
hist(Pofgb$Prop)
shapiro.test(Pofgb$Prop)
KPofgb<-kruskal.test(Pofgb$Prop~Pofgb$Suc.Stage)
KPofgb
pairwise.t.test(Pofgb$Prop,Pofgb$Suc.Stage, p.adj="bonf")#No differences
#Scavengers
Scfgb<-subset(fgb, Guild=="p.Sc")
hist(Scfgb$Prop)
shapiro.test(Scfgb$Prop)
KScfgb<-kruskal.test(Scfgb$Prop~Scfgb$Suc.Stage)
KScfgb
pairwise.t.test(Scfgb$Prop,Scfgb$Suc.Stage, p.adj="bonf")#No differences
#Grazers
Gzfgb<-subset(fgb, Guild=="p.Gz")#no Grazers
#To order them in the figure
fgb$Suc.Stage <- factor(fgb$Suc.Stage,levels=c('ES','YSF','MSF','OGF'))
fgb$Guild <- factor(fgb$Guild, levels = c('p.In','p.SD','p.Gr','p.Pr','p.Po','p.Sc','p.Gz'))
#~~~~~Mammals
fgm<-subset(fungui, Taxa =="Mammals")
fgmF<-subset(fgm, Suc.Stage!="ES")
#~as ES does have just one study, We did not include in the analysis
summary(fgmF)
#Insectivorous
Ifgm<-subset(fgmF, Guild=="p.In")
hist(Ifgm$Prop)
shapiro.test(Ifgm$Prop)
KIfgm<-kruskal.test(Ifgm$Prop~Ifgm$Suc.Stage)
KIfgm
pairwise.t.test(Ifgm$Prop,Ifgm$Suc.Stage, p.adj="bonf")#No differences
#Seed dispersers
SDfgm<-subset(fgmF, Guild=="p.SD")
hist(SDfgm$Prop)
shapiro.test(SDfgm$Prop)
KSDfgm<-kruskal.test(SDfgm$Prop~SDfgm$Suc.Stage)
KSDfgm
pairwise.t.test(SDfgm$Prop,SDfgm$Suc.Stage, p.adj="bonf")
#Granivorous
Grfgm<-subset(fgmF, Guild=="p.Gr")
hist(Grfgm$Prop)
shapiro.test(Grfgm$Prop)
KGrfgm<-kruskal.test(Grfgm$Prop~Grfgm$Suc.Stage)
KGrfgm
pairwise.t.test(Grfgm$Prop,Grfgm$Suc.Stage, p.adj="bonf")
#Pollinators
Pofgm<-subset(fgmF, Guild=="p.Po")
hist(Pofgm$Prop)
shapiro.test(Pofgm$Prop)
KPofgm<-kruskal.test(Pofgm$Prop~Pofgm$Suc.Stage)
KPofgm
pairwise.t.test(Pofgm$Prop,Pofgm$Suc.Stage, p.adj="bonf")
#Predators
Prfgm<-subset(fgmF, Guild=="p.Pr")
hist(Prfgm$Prop)
shapiro.test(Prfgm$Prop)
KPrfgm<-kruskal.test(Prfgm$Prop~Prfgm$Suc.Stage)
KPrfgm
pairwise.t.test(Prfgm$Prop,Prfgm$Suc.Stage, p.adj="bonf")
#Scavengers
Scfgm<-subset(fgmF, Guild=="p.Sc")
hist(Scfgm$Prop)
shapiro.test(Scfgm$Prop)
KScfgm<-kruskal.test(Scfgm$Prop~Scfgm$Suc.Stage)
KScfgm
pairwise.t.test(Scfgm$Prop,Scfgm$Suc.Stage, p.adj="bonf")
#Grazers
Gzfgm<-subset(fgmF, Guild=="p.Gz")
hist(Gzfgm$Prop)
shapiro.test(Gzfgm$Prop)
KGzfgm<-kruskal.test(Gzfgm$Prop~Gzfgm$Suc.Stage)
KGzfgm
pairwise.t.test(Gzfgm$Prop,Gzfgm$Suc.Stage, p.adj="bonf")
#To order in graph
fgm$Suc.Stage <- factor(fgm$Suc.Stage,levels=c('ES','YSF','MSF','OGF'))
fgm$Guild <- factor(fgm$Guild, levels = c('p.In','p.SD','p.Gr','p.Pr','p.Po','p.Sc','p.Gz'))
#~~~Box plot-Figure of the Functional groups of Birds and Mammals 
ggplot(fungui, aes(x=Guild, y=Prop, alpha=Suc.Stage))+
  geom_boxplot(data=fgb, fill="grey15", colour="black",
               outlier.colour = "grey15", outlier.shape = 1)+
  geom_boxplot(data=fgm, fill="grey15", colour="black",
               outlier.colour = "grey15", outlier.shape = 1)+
  facet_wrap(~Taxa, ncol=1, nrow=2)+
  ylab("Proportion of species \n")+  
  xlab("Functional groups \n ")+
  scale_alpha_discrete(range=c(0, 0.75), name="", labels=c("ES", 
                                                           "YSF",
                                                           "MSF",
                                                           "Ref."))+
  scale_x_discrete(labels = c("\n Insectivores", 
                              "\n Seed dispersers", 
                              "\n Granivores",
                              "\n Predators", 
                              "\n Pollinators", 
                              "\n Scavengers", 
                              "\n Grazers"))+ 
  scale_y_continuous(breaks = seq(0,1.1, by=0.25),
                     labels = seq(0,1.1, by=0.25),
                     limits = c(0,1.1))+
  theme_bw()+
  theme(axis.text=element_text(size=22), axis.title=element_text(size=25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(legend.text = element_text(colour=rgb(0,0,0), size = 18),
        legend.key.size=unit(.8, "cm"))+
  theme(legend.background = element_rect( ))+
  theme(legend.background = element_rect(fill=rgb(1,1,1), size=0.5, linetype=2))+
  theme(legend.position=c(.94,.835))+
  theme(panel.margin = unit(2, "lines"))+
  theme(strip.background = element_blank())+
  theme(strip.text = element_text(size=20))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

