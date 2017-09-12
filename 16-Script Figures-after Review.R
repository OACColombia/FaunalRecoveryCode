#Figure 2
#Packages
#packages used
library(ggplot2)
library(reshape2)
library(gridExtra)
library(gtable)
library(grid)
library(ggrepel)
library(plyr)
library(png)
setwd("/Users/orlando/Dropbox/Vertebrate Recovery/data")

#call the table saved
dat <- read.csv("//Users/orlando/Dropbox/Vertebrate Recovery/data/datbootsRRR.csv")
head(dat)
ASil<-readPNG('ASil.png')
RSil<-readPNG('RSil.png')
BSil<-readPNG('BSil.png')
MSil<-readPNG('MSil.png')

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
  annotation_raster(ASil, ymin = -1.25,ymax= -1.4,xmin = 3.9,xmax = 4.1)+
  annotation_raster(RSil, ymin = -1.25,ymax= -1.4,xmin = 2.8,xmax = 3.2)+
  annotation_raster(BSil, ymin = -1.2,ymax= -1.4,xmin = 1.85,xmax = 2.2)+
  annotation_raster(MSil, ymin = -1.2,ymax= -1.4,xmin = 0.85,xmax = 1.2)+
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
dat <- read.csv("//Users/orlando/Dropbox/Vertebrate Recovery/data/datbootsRRS.csv")
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
  geom_hline(yintercept = -0.511, linetype = "dashed",color="gray", size=0.5)+
  geom_hline(yintercept = 0, linetype = "dashed",color="black", size=0.5)+
  geom_text(aes(label=ov$Label, y=-1.6, fill=Taxa), 
            position = position_dodge(width = 1), size=3.5)+
  #geom_pointrange(position = position_dodge(1.2), size=0.8,aes(fill=Taxa), colour="black", stroke=1)+
  geom_linerange(position = position_dodge(1), size=4, alpha=0.6, aes(color=Taxa))+
  geom_point(position = position_dodge(1), size=5, stroke = 1, alpha=0.8,  aes(fill=Taxa))+
  labs(title = "b",
  		subtitle = "Species compositional similarity")+
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
#Joint figures
ab=grid.arrange(arrangeGrob(arrangeGrob(a,b,ncol=2)),
                  left=textGrob("\nSuccessional stages by taxa\n", 
                                rot=90, vjust=1, 
                                  gp=gpar(fontsize=14)),
                  top=textGrob(""),
                  right=textGrob(""),
                  bottom=textGrob("Response ratio during secondary forest succession (bootstrapped effect size)\n", 
                                  gp=gpar(fontsize=14)))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Fit single model, age as single predictor for faunal recovery of species composition similarity
mydata <- read.csv("//Users/orlando/Dropbox/Vertebrate Recovery/data/MainData.csv")

##~~~~~~~~~~Forest specialist figure
#call the table saved
dat <- read.csv("//Users/orlando/Dropbox/Vertebrate Recovery/data/datbootsRRFSp.csv")
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
  annotation_raster(ASil, ymin = -2.2,ymax= -2.5,xmin = 3.8,xmax = 4.1)+
  annotation_raster(RSil, ymin = -2.25,ymax= -2.5,xmin = 2.7,xmax = 3.2)+
  annotation_raster(BSil, ymin = -2.2,ymax= -2.5,xmin = 1.85,xmax = 2.2)+
  annotation_raster(MSil, ymin = -2.2,ymax= -2.5,xmin = 0.85,xmax = 1.2)+
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
dat <- read.csv("//Users/orlando/Dropbox/Vertebrate Recovery/data/datbootsRRFunG.csv")
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
                  left=textGrob("\nFunctional groups", 
                                rot=90, vjust=1,
                                  gp=gpar(fontsize=14)),
                  top=textGrob(""),
                  right=textGrob(""),
                  bottom=textGrob("Response ratio during secondary forest succession \n (bootstrapped effect size)\n ",
                  					hjust=0.25,
                                  gp=gpar(fontsize=14)))



