#Figure 2

#call the table saved
dat <- read.csv("//Users/orlando/Dropbox/Vertebrate Recovery - Nature/data/datbootsRRR.csv")
head(dat)

dat$label<-paste(dat$Succ.Stage,"(",dat$n,")")
head(dat)

#Lets graph all this work
colnames(dat)<- c("Bootstrap","n","Bmean","lower","upper", "FigurePart","Taxa", "Succ.Stage","Label")
dat$Label

dat$Succ.Stage <- factor(dat$Succ.Stage,levels=c('ES','YSF','MSF', 'OSF'))
dat$Taxa <- factor(dat$Taxa,levels=c('Mammals','Birds','Reptiles','Amphibians'))

ov<-subset(dat,FigurePart=="Overall")
a<-ggplot(ov, aes(y = Bmean, ymin = lower, ymax = upper,
                        x = Taxa, shape=Succ.Stage, fill=Taxa))+
  geom_vline(xintercept=c(1.5,2.5,3.5),color="darkgray", size=0.4)+
  geom_hline(yintercept = 0, linetype = "dashed",color="black", size=0.5)+
  geom_text(aes(label=ov$Label, y=-1.6, fill=Taxa), 
            position = position_dodge(width = 1), size=3.5)+
  #geom_pointrange(position = position_dodge(1.2), size=0.8,aes(fill=Taxa), colour="black", stroke=1)+
  geom_linerange(position = position_dodge(1), size=4, alpha=0.6, aes(color=Taxa))+
  geom_point(position = position_dodge(1), size=5, stroke = 1, alpha=0.8,  aes(fill=Taxa))+
  labs(title = "a",
  		subtitle = "Species richness")+
#  ggtitle('a')+
  ylab("")+  
  xlab("")+
  scale_shape_manual(name="Succ.Stage",
                     values = c("ES" = 21, "YSF"=22, 
                                "MSF"=24,"OSF"=23))+
  scale_y_continuous(breaks = seq(-1.5,0.9, by=0.5),
                     labels = seq(-1.5,0.9, by=0.5),
                     limits = c(-1.75,1), expand = c(0, 0))+
  scale_x_discrete(breaks=c('Amphibians','Reptiles','Birds','Mammals'),
                   labels=c('Amphibians','Reptiles','Birds','Mammals'), 
                   expand = c(0, 0))+
  scale_fill_manual(name="Taxa",
                    values = c("Amphibians" = "#397d34", "Reptiles"="#FFff02", 
                               "Birds"="#1f78b4","Mammals"="#FF7f00"))+
  scale_color_manual(name="Taxa",
                     values = c("Amphibians" = "#397d34", "Reptiles"="#FFff02", 
                               "Birds"="#1f78b4","Mammals"="#FF7f00"))+
  theme_bw()+
  theme(legend.position="none",axis.text.x=element_text(size=12, hjust = 0.5), axis.text.y=element_blank(), 
        plot.margin=unit(c(1,1,-2,1),"mm"),panel.margin.y = unit(0, "lines"), panel.grid.major = element_blank(), panel.grid.minor = 	element_blank())+ 
  theme(strip.background = element_blank(),plot.title=element_text(hjust=-0.025, size = 16, face = "bold"),
  		plot.subtitle=element_text(hjust=0.5,size=14),
        strip.text = element_blank())+
  coord_flip()+
  annotate("text",x=4,y=0.6,label = "Null", size=4,hjust=0.5)+
  annotate("text",x=3,y=0.6,label = "Age ↑", size=4,hjust=0.5)+
  annotate("text",x=2,y=0.6,label = "Geo. Condition - Islands ↓ \n Rainfall ↑ \n Age ↑", size=4,hjust=0.5)+
  annotate("text",x=1,y=0.6,label = "Prev. Land - Pasture ↓ \n Rainfall ↓", size=4,hjust=0.5)
a

#call the table saved for species composition
dat <- read.csv("//Users/orlando/Dropbox/Vertebrate Recovery - Nature/data/datbootsRRS.csv")
head(dat)

dat$label<-paste(dat$Succ.Stage,"(",dat$n,")")
head(dat)

#Lets graph all this work
colnames(dat)<- c("Bootstrap","n","Bmean","lower","upper", "FigurePart","Taxa", "Succ.Stage", "Label")
dat$Label

dat$Succ.Stage <- factor(dat$Succ.Stage,levels=c('ES','YSF','MSF', 'OSF'))
dat$Taxa <- factor(dat$Taxa,levels=c('Mammals','Birds','Reptiles','Amphibians'))

ov<-subset(dat,FigurePart=="Overall")
b<-ggplot(ov, aes(y = Bmean, ymin = lower, ymax = upper,
                        x = Taxa, shape=Succ.Stage, fill=Taxa))+
  geom_vline(xintercept=c(1.5,2.5,3.5),color="darkgray", size=0.4)+
  geom_hline(yintercept = 0, linetype = "dashed",color="black", size=0.5)+
  geom_text(aes(label=ov$Label, y=-1.6, fill=Taxa), 
            position = position_dodge(width = 1), size=3.5)+
  #geom_pointrange(position = position_dodge(1.2), size=0.8,aes(fill=Taxa), colour="black", stroke=1)+
  geom_linerange(position = position_dodge(1), size=4, alpha=0.6, aes(color=Taxa))+
  geom_point(position = position_dodge(1), size=5, stroke = 1, alpha=0.8,  aes(fill=Taxa))+
  labs(title = "b",
  		subtitle = "Species composition similarity")+
  ylab("")+  
  xlab("")+
  scale_shape_manual(name="Succ.Stage",
                     values = c("ES" = 21, "YSF"=22, 
                                "MSF"=24,"OSF"=23))+
  scale_y_continuous(breaks = seq(-1.5,0.1, by=0.5),
                     labels = seq(-1.5,0.1, by=0.5),
                     limits = c(-1.75,0.7), expand = c(0, 0))+
  scale_x_discrete(breaks=c('Amphibians','Reptiles','Birds','Mammals'),
                   labels=c('Amphibians','Reptiles','Birds','Mammals'), 
                   expand = c(0, 0))+
  scale_fill_manual(name="Taxa",
                    values = c("Amphibians" = "#397d34", "Reptiles"="#FFff02", 
                               "Birds"="#1f78b4","Mammals"="#FF7f00"))+
  scale_color_manual(name="Taxa",
                     values = c("Amphibians" = "#397d34", "Reptiles"="#FFff02", 
                               "Birds"="#1f78b4","Mammals"="#FF7f00"))+
  theme_bw()+
  theme(legend.position="none",axis.text.x=element_text(size=12, hjust = 0.5), axis.text.y=element_blank(), 
        plot.margin=unit(c(1,1,-2,1),"mm"),panel.margin.y = unit(0, "lines"), panel.grid.major = element_blank(), panel.grid.minor = 	element_blank())+ 
  theme(strip.background = element_blank(),plot.title=element_text(hjust=-0.025, size = 16, face = "bold"),
  		plot.subtitle=element_text(hjust=0.5,size=14),
        strip.text = element_blank())+
  coord_flip()+
  annotate("text",x=4,y=0.35,label = "Suc. Stage - MSF ↑ \n Landscape - PatchPtoA ↓   \n Geo. Condition - Islands ↑\n Landscape - nPatches ↑\n Rainfall ↓", size=4,hjust=0.5)+
  annotate("text",x=3,y=0.35,label = "Geo. Condition - Islands ↑ \n Landscape - nPatches ↑", size=4,hjust=0.5)+
  annotate("text",x=2,y=0.35,label = "Landscape - PatchPtoA ↑  \n Age ↑ \n Elevation ↑\n Landscape - PatchPtoA ↑", size=4,hjust=0.5)+
  annotate("text",x=1,y=0.35,label = "Rainfall ↓ \n Landscape - nPatches ↑ \n Landscape - p.forest ↓", size=4,hjust=0.5)
b

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Fit single model, age as single predictor for faunal recovery of species composition similarity
mydata <- read.csv("//Users/orlando/Dropbox/Vertebrate Recovery - Nature/data/MainData.csv")
#packages used
library(ggplot2)
library(reshape2)
library(gridExtra)
library(gtable)
library(grid)
library(ggrepel)
library(plyr)

#Exclude the forest with zero richness (the RRS can't be calculated)
data<-subset(mydata,RRS!="NA")

summary(data$Si)
attach(data)
head(data)
nrow(data)#186 observations
#SS
data <- data[!apply(data[,c("Si", 
                            "Mean.Age", 
                            "Taxa", 
                            "Biome")], 1, anyNA),]
nrow(data)#161, We lost rows of studies without list or NA value of similarity
attach(data)
#~~~Simple model Amphibians
A<-subset(data, Taxa=="Amphibians")
nrow(A)#42 observations
length(unique(A$Nstudy))#24 studies
SSA<-lm(A$Si~1+A$Mean.Age)
summary(SSA)
A <- data.frame(A, predict(SSA, A, interval="confidence"))
head(A)#fit = model; lwr=lower confidence; upr=upper confidence
A <- data.frame(A, predict(SSA, A, interval="predict"))
head(A)#fit.1 = model; lwr.1=lower prediction; upr.1=upper prediction
upPIA<-lm(A$upr.1~1+A$Mean.Age)#upper interval line
lwPIA<-lm(A$lwr.1~1+A$Mean.Age)#lower interval line
#~~~~~Reptiles
R<-subset(data, Taxa=="Reptiles")
nrow(R)#31 observations
length(unique(R$Nstudy))#16 studies
SSR<-lm(R$Si~1+R$Mean.Age)
summary(SSR)
R <- data.frame(R, predict(SSR, R, interval="confidence"))
head(R)#fit = model; lwr=lower confidence; upr=upper confidence
R <- data.frame(R, predict(SSR, R, interval="predict"))
head(R)#fit.1 = model; lwr.1=lower prediction; upr.1=upper prediction
upPIR<-lm(R$upr.1~1+R$Mean.Age)#upper interval line
lwPIR<-lm(R$lwr.1~1+R$Mean.Age)#lower interval line
#~~~~~Birds
B<-subset(data, Taxa=="Birds")
nrow(B)#49 observations
length(unique(B$Nstudy))#31 studies
SSB<-lm(B$Si~1+B$Mean.Age)
summary(SSB)
B <- data.frame(B, predict(SSB, B, interval="confidence"))
head(R)#fit = model; lwr=lower confidence; upr=upper confidence
B <- data.frame(B, predict(SSB, B, interval="predict"))
head(B)#fit.1 = model; lwr.1=lower prediction; upr.1=upper prediction
upPIB<-lm(B$upr.1~1+B$Mean.Age)#upper interval line
lwPIB<-lm(B$lwr.1~1+B$Mean.Age)#lower interval line


#~~~~~~Mammals
M<-subset(data, Taxa=="Mammals")
nrow(M)#39 observations
length(unique(M$Nstudy))#26 studies
SSM<-lm(M$Si~1+M$Mean.Age)
summary(SSM)
M <- data.frame(M, predict(SSM, M, interval="confidence"))
head(M)#fit = model; lwr=lower confidence; upr=upper confidence
M <- data.frame(M, predict(SSM, M, interval="predict"))
head(M)#fit.1 = model; lwr.1=lower prediction; upr.1=upper prediction
upPIM<-lm(M$upr.1~1+M$Mean.Age)#upper interval line
lwPIM<-lm(M$lwr.1~1+M$Mean.Age)#lower interval line

#Explaniation of the variation by age
SSA2<-anova(SSA)
SSA2ss<-SSA2$"Sum Sq"
print(cbind(SSA2,PctExp=SSA2ss/sum(SSA2ss)*100)) #PctExp: % explain variation

SSR2<-anova(SSR)
SSR2ss<-SSR2$"Sum Sq"
print(cbind(SSR2,PctExp=SSR2ss/sum(SSR2ss)*100)) #PctExp: % explain variation

SSB2<-anova(SSB)
SSB2ss<-SSB2$"Sum Sq"
print(cbind(SSB2,PctExp=SSB2ss/sum(SSB2ss)*100)) #PctExp: % explain variation

SSM2<-anova(SSM)
SSM2ss<-SSM2$"Sum Sq"
print(cbind(SSM2,PctExp=SSM2ss/sum(SSM2ss)*100)) #PctExp: % explain variation


#~~~Figure fit null model Sorensen Similarity and prediction time of recovery
library(ggplot2)
data$Taxa <- factor(data$Taxa,levels=c('Amphibians','Reptiles','Birds', 'Mammals'))
data$SuccStage <- factor(data$SuccStage,levels = c('ES','YSF','MSF','OSF'))
##changing the graph
#1~ points
p<-ggplot(data,aes(Mean.Age,Si, 
                   shape=SuccStage, fill=Taxa))+
  geom_point(size = 3, stroke = 1, alpha=0.8, colour="black")
p
#2~some aesthetics
paes<-p+
  theme_bw()+
  ylab("Similarity to a reference forest \n (Sørensen index, SS)")+  
  xlab("\n Age (years after abandonment)")+
  ggtitle('c')+
   scale_shape_manual(name="Succ.Stage",
                     values = c("ES" = 21, "YSF"=22, 
                                "MSF"=24,"OSF"=23))+
  scale_x_continuous(breaks = seq(0,150,by=30),
                     labels = seq(0,150,by=30),
                     limits = c(0,160), expand = c(0.02,0))+ 
  scale_y_continuous(breaks = seq(0,1.2, by=0.5),
                     labels = seq(0,1.2, by=0.5),
                     limits = c(0,1.4),expand = c(0.02,0))+
  scale_fill_manual(name="Taxa",
                    values = c("Amphibians" = "#397d34", "Reptiles"="#FFff02", 
                               "Birds"="#1f78b4","Mammals"="#FF7f00"))+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.position="none")+
  facet_wrap(~Taxa, ncol=2, nrow=2, scales = "free_x")+ 
  theme(plot.margin = unit(c(2,3,2,3),"mm"))+
  theme(panel.margin =  unit(1, "lines"),plot.title=element_text(hjust=-0.05, size = 16, face = "bold"))+
  theme(strip.background = element_blank())+
  theme(strip.text = element_text(size=14))
 paes

#~3 fit the linear model
lpi<-paes+
  geom_hline(yintercept=0.9, linetype="dashed", colour="black", size=0.5)+
  geom_segment(data=A, aes(x=0, xend=((0.9-SSA$coefficients[1])/SSA$coefficients[2]), 
                           y=(SSA$coefficients[1]),yend=0.9),
               colour="black", size=1.5,linetype="dotted")+
  geom_segment(data=A, aes(x=0, xend=((0.9-lwPIA$coefficients[1])/lwPIA$coefficients[2]), 
                           y=(lwPIA$coefficients[1]),yend=0.9),
               colour="black", size=1,linetype="dotted")+
  geom_segment(data=A, aes(x=0, xend=((0.9-upPIA$coefficients[1])/upPIA$coefficients[2]), 
                           y=(upPIA$coefficients[1]),yend=0.9),
               colour="black", size=1,linetype="dotted")+
  geom_segment(data=R, aes(x=0, xend=((0.9-SSR$coefficients[1])/SSR$coefficients[2]), 
                           y=(SSR$coefficients[1]),yend=0.9),
               colour="black", size=1.5,linetype="dotted")+
  geom_segment(data=R, aes(x=0, xend=((0.9-lwPIR$coefficients[1])/lwPIR$coefficients[2]), 
                           y=(lwPIR$coefficients[1]),yend=0.9),
               colour="black", size=1,linetype="dotted")+
  geom_segment(data=R, aes(x=0, xend=((0.9-upPIR$coefficients[1])/upPIR$coefficients[2]), 
                           y=(upPIR$coefficients[1]),yend=0.9),
               colour="black", size=1,linetype="dotted")+
  geom_segment(data=B, aes(x=0, xend=((0.9-SSB$coefficients[1])/SSB$coefficients[2]), 
                           y=(SSB$coefficients[1]),yend=0.9),
               colour="black", size=1.5,linetype="solid")+
  geom_segment(data=B, aes(x=0, xend=((0.9-lwPIB$coefficients[1])/lwPIB$coefficients[2]), 
                           y=(lwPIB$coefficients[1]),yend=0.9),
               colour="#1f78b4", size=0.5,linetype="solid")+
  geom_segment(data=B, aes(x=0, xend=((0.9-upPIB$coefficients[1])/upPIB$coefficients[2]), 
                           y=(upPIB$coefficients[1]),yend=0.9),
               colour="#1f78b4", size=0.5,linetype="solid")+
  geom_segment(data=M, aes(x=0, xend=((0.9-SSM$coefficients[1])/SSM$coefficients[2]), 
                           y=(SSM$coefficients[1]),yend=0.9),
               colour="black", size=1.5,linetype="dotted")+
  geom_segment(data=M, aes(x=0, xend=((0.9-lwPIM$coefficients[1])/lwPIM$coefficients[2]), 
                           y=(lwPIM$coefficients[1]),yend=0.9),
               colour="black", size=1,linetype="dotted")+
  geom_segment(data=M, aes(x=0, xend=((0.9-upPIM$coefficients[1])/upPIM$coefficients[2]), 
                           y=(upPIM$coefficients[1]),yend=0.9),
               colour="black", size=1,linetype="dotted")
lpi

#prediction interval or the time of recovery
c<-lpi+ 
  geom_segment(data=B, aes(x=((0.9-upPIB$coefficients[1])/upPIB$coefficients[2]),
                           xend=((0.9-lwPIB$coefficients[1])/lwPIB$coefficients[2]), y=1.2, yend=1.2), 
               linetype="solid", colour="#1f78b4", size=2)+
  geom_point(data=B, aes(x=((0.9-SSB$coefficients[1])/SSB$coefficients[2]), 
                         xend=((0.9-SSB$coefficients[1])/SSB$coefficients[2]), y=1.2, yend=1.2), 
             colour="black", fill="black",shape=21, size=6)
c#Figure of the recovery time
#Joint figures
ab=grid.arrange(arrangeGrob(arrangeGrob(a,b,ncol=2)),
                  left=textGrob("Successional stages by taxa", 
                                rot=90, vjust=1, 
                                  gp=gpar(fontsize=14)),
                  top=textGrob(""),
                  right=textGrob(""),
                  bottom=textGrob("Response ratio during secondary forest succession \n (bootstrapped effect size)", 
                                  gp=gpar(fontsize=14)))
#Figure2
grid.arrange(ab,c,nrow=2,heights = c(1.75,0.75))


##~~~~~~~~~~Forest specialist figure
#call the table saved
dat <- read.csv("//Users/orlando/Dropbox/Vertebrate Recovery - Nature/data/datbootsRRFSp.csv")
head(dat)

dat$label<-paste(dat$Succ.Stage,"(",dat$n,")")
head(dat)


#Lets graph all this work
colnames(dat)<- c("Bootstrap","n","Bmean","lower","upper", "FigurePart","Taxa", "Succ.Stage","Label")
dat$Label

dat$Succ.Stage <- factor(dat$Succ.Stage,levels=c('ES','YSF','MSF', 'OSF'))
dat$Taxa <- factor(dat$Taxa,levels=c('Mammals','Birds','Reptiles','Amphibians'))

Fig3<-ggplot(dat, aes(y = Bmean, ymin = lower, ymax = upper,
                        x = Taxa, shape=Succ.Stage, fill=Taxa))+
  geom_vline(xintercept=c(1.5,2.5,3.5),color="darkgray", size=0.4)+
  geom_hline(yintercept = 0, linetype = "dashed",color="black", size=0.5)+
  geom_text(aes(label=dat$Label, y=-3.25, fill=Taxa), 
            position = position_dodge(width = 1), size=3.5)+
  #geom_pointrange(position = position_dodge(1.2), size=0.8,aes(fill=Taxa), colour="black", stroke=1)+
  geom_linerange(position = position_dodge(1), size=4, alpha=0.6, aes(color=Taxa))+
  geom_point(position = position_dodge(1), size=5, stroke = 1, alpha=0.8,  aes(fill=Taxa))+
  labs(title = "",
  		subtitle = "Forest specialist species (tropical moist forest)")+
 ylab("\nResponse ratio during secondary forest succession \n (bootstrapped effect size)\n")+  
  xlab("\nSuccessional stages by taxa \n")+
  scale_shape_manual(name="Succ.Stage",
                     values = c("ES" = 21, "YSF"=22, 
                                "MSF"=24,"OSF"=23))+
  scale_y_continuous(breaks = seq(-3,0.9, by=0.5),
                     labels = seq(-3,0.9, by=0.5),
                     limits = c(-3.5,1), expand = c(0, 0))+
  scale_x_discrete(breaks=c('Amphibians','Reptiles','Birds','Mammals'),
                   labels=c('Amphibians','Reptiles','Birds','Mammals'), 
                   expand = c(0.025, 0))+
  scale_fill_manual(name="Taxa",
                    values = c("Amphibians" = "#397d34", "Reptiles"="#FFff02", 
                               "Birds"="#1f78b4","Mammals"="#FF7f00"))+
  scale_color_manual(name="Taxa",
                     values = c("Amphibians" = "#397d34", "Reptiles"="#FFff02", 
                               "Birds"="#1f78b4","Mammals"="#FF7f00"))+
  theme_bw()+
  theme(legend.position="none",axis.text.x=element_text(size=12, hjust = 0.5), axis.text.y=element_blank(), 
        plot.margin=unit(c(2,2,2,2),"mm"),panel.margin.y = unit(0, "lines"), panel.grid.major = element_blank(), panel.grid.minor = 	element_blank(), axis.title=element_text(size=14))+ 
  theme(strip.background = element_blank(),plot.title=element_text(hjust=-0.025, size = 16, face = "bold"),
  		plot.subtitle=element_text(hjust=0.5,size=14),
        strip.text = element_blank())+
  coord_flip()+
  annotate("text",x=4,y=0.6,label = "Null", size=4,hjust=0.5)+
  annotate("text",x=3,y=0.6,label = "Null", size=4,hjust=0.5)+
  annotate("text",x=2,y=0.6,label = "Age ↑", size=4,hjust=0.5)+
  annotate("text",x=1,y=0.6,label = "Rainfall ↓", size=4,hjust=0.5)
Fig3

##~~~~~~~~~~Functional Groups
#call the table saved
dat <- read.csv("//Users/orlando/Dropbox/Vertebrate Recovery - Nature/data/datbootsRRFunG.csv")
head(dat)

dat$label<-paste(dat$Succ.Stage,"(",dat$n,")")
head(dat)


#Lets graph all this work
colnames(dat)<- c("Bootstrap","n","Bmean","lower","upper", "FigurePart","Taxa", "Succ.Stage","Label")
dat$Label
birds<-subset(dat,Taxa=="Birds")

birds$FigurePart <- factor(birds$FigurePart,levels=c('Gz','Sc','Po','Pr','Gr','SD','In'))
birds$Succ.Stage <- factor(birds$Succ.Stage,levels=c('ES','YSF','MSF', 'OSF'))
Fb<-ggplot(birds,aes(y=Bmean,ymin=lower,ymax=upper,x=FigurePart,shape=Succ.Stage,fill=Taxa))+
 geom_text(aes(label=birds$Label, y=-1.7, fill=Taxa,shape=Succ.Stage,x=FigurePart), 
            position = position_dodge(width = 1), size=3.5)+
   geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5,6.5),color="darkgray", size=0.4)+
  labs(title = "a",
  		subtitle = "Birds (tropical moist forest)")+
  geom_hline(yintercept = 0, linetype = "dashed",color="black", size=0.5)+
  geom_linerange(position = position_dodge(1), size=4, alpha=0.6, aes(color=Taxa))+
  geom_point(position = position_dodge(1), size=5, stroke = 1, alpha=0.8,  aes(fill=Taxa))+
  scale_x_discrete(breaks=c('In','SD','Gr','Pr','Po','Sc','Gz'),
                   labels=c('Insectivores','Seed Dispersers','Granivores','Predators','Pollinators','Scavengers','Grazers'), 
                   expand = c(0, 0))+
  ylab("")+  
  xlab("")+
  scale_shape_manual(name="Succ.Stage",
                     values = c("ES" = 21, "YSF"=22, 
                                "MSF"=24,"OSF"=23))+
  scale_y_continuous(breaks = seq(-1.5,1.4, by=0.5),
                     labels = seq(-1.5,1.4, by=0.5),
                     limits = c(-1.95,1.7), expand = c(0, 0))+
    scale_fill_manual(name="Taxa",
                    values = c("Amphibians" = "#397d34", "Reptiles"="#FFdd02", 
                               "Birds"="#1f78b4","Mammals"="#FF7f00"))+
  scale_color_manual(name="Taxa",
                     values = c("Amphibians" = "#397d34", "Reptiles"="#FFdd02", 
                               "Birds"="#1f78b4","Mammals"="#FF7f00"))+
  theme_bw()+
  coord_flip()+
    theme(legend.position="none",axis.text.x=element_text(size=12, hjust = 0.5), axis.text.y=element_text(face="bold",size=12), 
        plot.margin=unit(c(1,1,1,1),"mm"),panel.margin.y = unit(0, "lines"), panel.grid.major = element_blank(), panel.grid.minor =	element_blank())+ 
  theme(strip.background = element_blank(),plot.title=element_text(hjust=-0.025, size = 16, face = "bold"),
  		plot.subtitle=element_text(hjust=0.5,size=14))+
  theme(strip.text = element_blank())+
  annotate("text",x=7,y=1.1,label = "Null", size=4,hjust=0.5)+
  annotate("text",x=6,y=1.1,label = "Landscape - Patch Area ↓ \n Age ↑ \n Elevation ↑", size=4,hjust=0.5)+
  annotate("text",x=5,y=1.1,label = "Landscape - p.Forest ↓ \n Landscape - PatchPtoA ↓", size=4,hjust=0.5)+
  annotate("text",x=4,y=1.1,label = "Null", size=4,hjust=0.5)+
  annotate("text",x=3,y=1.1,label = "Null", size=4,hjust=0.5)+
  annotate("text",x=2,y=1.1,label = "N.A.", size=4,hjust=0.5)+
  annotate("text",x=1,y=1.1,label = "N.A.", size=4,hjust=0.5)
Fb
#And mammals
mam<-subset(dat,Taxa=="Mammals")
mam$FigurePart <- factor(mam$FigurePart,levels=c('Gz','Sc','Po','Pr','Gr','SD','In'))
mam$Succ.Stage <- factor(mam$Succ.Stage,levels=c('ES','YSF','MSF', 'OSF'))
Fm<-ggplot(mam,aes(y=Bmean,ymin=lower,ymax=upper,x=FigurePart,shape=Succ.Stage,fill=Taxa))+
 geom_text(aes(label=mam$Label, y=-1.7, fill=Taxa,shape=Succ.Stage,x=FigurePart), 
            position = position_dodge(width = 1), size=3.5)+
   geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5,6.5),color="darkgray", size=0.4)+
  labs(title = "b",
  		subtitle = "Mammals (tropical moist forest)")+
  geom_hline(yintercept = 0, linetype = "dashed",color="black", size=0.5)+
  geom_linerange(position = position_dodge(1), size=4, alpha=0.6, aes(color=Taxa))+
  geom_point(position = position_dodge(1), size=5, stroke = 1, alpha=0.8,  aes(fill=Taxa))+
  scale_x_discrete(breaks=c('In','SD','Gr','Pr','Po','Sc','Gz'),
                   labels=c('','','','','','',''), 
                   expand = c(0, 0))+
  ylab("")+  
  xlab("")+
  scale_shape_manual(name="Succ.Stage",
                     values = c("ES" = 21, "YSF"=22, 
                                "MSF"=24,"OSF"=23))+
  scale_y_continuous(breaks = seq(-1.5,1.4, by=0.5),
                     labels = seq(-1.5,1.4, by=0.5),
                     limits = c(-1.95,1.7), expand = c(0, 0))+
    scale_fill_manual(name="Taxa",
                    values = c("Amphibians" = "#397d34", "Reptiles"="#FFdd02", 
                               "Birds"="#1f78b4","Mammals"="#FF7f00"))+
  scale_color_manual(name="Taxa",
                     values = c("Amphibians" = "#397d34", "Reptiles"="#FFdd02", 
                               "Birds"="#1f78b4","Mammals"="#FF7f00"))+
  theme_bw()+
  coord_flip()+
    theme(legend.position="none",axis.text.x=element_text(size=12, hjust = 0.5), #axis.text.y=element_blank(), 
        plot.margin=unit(c(1,1,1,1),"mm"),panel.margin.y = unit(0, "lines"), panel.grid.major = element_blank(), panel.grid.minor =	element_blank())+ 
  theme(strip.background = element_blank(),plot.title=element_text(hjust=-0.025, size = 16, face = "bold"),
  		plot.subtitle=element_text(hjust=0.5,size=14))+
  theme(strip.text = element_blank())+
  annotate("text",x=7,y=1,label = "Null", size=4,hjust=0.5)+
  annotate("text",x=6,y=1,label = "Rainfall ↓ \n Landscape - PatchPtoA ↓ \n Landscape - p.Forest ↓", size=4,hjust=0.5)+
  annotate("text",x=5,y=1,label = "Null", size=4,hjust=0.5)+
  annotate("text",x=4,y=1,label = "N.A.", size=4,hjust=0.5)+
  annotate("text",x=3,y=1,label = "N.A.", size=4,hjust=0.5)+
  annotate("text",x=2,y=1,label = "N.A.", size=4,hjust=0.5)+
  annotate("text",x=1,y=1,label = "Landscape - PatchPtoA ↓ \n Geo.Condition - Islands ↓", size=4,hjust=0.5)
Fm

#Joint figures
FSp=grid.arrange(arrangeGrob(Fb,Fm,ncol=2,widths = c(1,0.85)),
                  left=textGrob("Functional groups", 
                                rot=90, vjust=1,
                                  gp=gpar(fontsize=14)),
                  top=textGrob(""),
                  right=textGrob(""),
                  bottom=textGrob("Response ratio during secondary forest succession \n (bootstrapped effect size)\n ",
                  					hjust=0.25,
                                  gp=gpar(fontsize=14)))



