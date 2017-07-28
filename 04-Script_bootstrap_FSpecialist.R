#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Bootstrapping of the Response ratio in Species richness
#Clear memory
rm(list=ls())
#Working in the office
setwd("//Users/orlando/Dropbox/Vertebrate Recovery - Nature/data")
#Call the file MainData.csv
mydata <- read.csv("//Users/orlando/Dropbox/Vertebrate Recovery - Nature/data/MainData.csv")

#Working in my house
setwd("c://Users/Flaco/Dropbox/Vertebrate Recovery - Nature/data")
#Call the file MainData.csv
mydata <- read.csv("c://Users/Flaco/Dropbox/Vertebrate Recovery - Nature/data/MainData.csv")
#packages used
library(ggplot2)
library(reshape2)
library(gridExtra)
library(gtable)
library(grid)
library(ggrepel)
library(plyr)


data<-subset(mydata,Biome=="Moist")

#Exclude the forest with zero richness (the RRFSp can't be calculated)
data<-subset(mydata,RRFSp!="NA")
nrow(data)#136
head(data)

#To select one comparison per study (avoid spacial pseudoreplication)
randomRows= function(df,n){
  return(df[sample(nrow(df),n),])
}

#B1: Amphibians in Early Succession (ES), moist forest specialist species
AmES<-subset(data,BGroupOverall=="AmphibiansES")
nrow(AmES)
median(AmES$RRFSp)
R<-10000
B1<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  AmES2<-ddply(AmES,.(Nstudy))
  boot.sample <- sample(AmES2$RRFSp, replace = TRUE)
  B1[kk] <- mean(boot.sample)
}
boxplot(B1)
quantile(B1,c(0.025,0.975))
hist(B1, breaks = 30)
dB1<-c(Bmean=mean(B1),quantile(B1,c(0.025,0.975)),n=nrow(AmES2))
dB1
#B2: Amphibians in Young Secondary Forest (YSF), overall
AmYSF<-subset(data,BGroupOverall=="AmphibiansYSF")
median(AmYSF$RRFSp)
nrow(AmYSF)
R<-10000
B2<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  AmYSF2<-ddply(AmYSF,.(Nstudy))
  boot.sample <- sample(AmYSF2$RRFSp, replace = TRUE)
  B2[kk] <- mean(boot.sample)
}
boxplot(B2)
quantile(B2,c(0.025,0.975))
hist(B2, breaks = 30)
dB2<-c(Bmean=mean(B2),quantile(B2,c(0.025,0.975)),n=nrow(AmYSF2))
dB2
#B3: Amphibians in Mid-successional Secondary Forest (MSF), overall
AmMSF<-subset(data,BGroupOverall=="AmphibiansMSF")
median(AmMSF$RRFSp)
nrow(AmMSF)
R<-10000
B3<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  AmMSF2<-ddply(AmMSF,.(Nstudy))
  boot.sample <- sample(AmMSF2$RRFSp, replace = TRUE)
  B3[kk] <- mean(boot.sample)
}
boxplot(B3)
quantile(B3,c(0.025,0.975))
hist(B3, breaks = 30)
dB3<-c(Bmean=mean(B3),quantile(B3,c(0.025,0.975)),n=nrow(AmMSF2))
dB3
#B4: Amphibians in Old Secondary Forest (OSF), overall
AmOSF<-subset(data,BGroupOverall=="AmphibiansOSF")
nrow(AmOSF)
B4<-rep(mean(AmOSF$RRFSp),10000)
dB4<-c(Bmean=mean(B4),"2.5%" = NaN,"97.5%" = NaN,n=nrow(AmOSF))
dB4
#B5: Reptiles in ES, overall
ReES<-subset(data,BGroupOverall=="ReptilesES")
nrow(ReES)
B5<-rep(mean(ReES$RRFSp),10000)
dB5<-c(Bmean=mean(B5),"2.5%" = NaN,"97.5%" = NaN,n=nrow(ReES))
dB5
#B6: Reptiles in YSF, overall
ReYSF<-subset(data,BGroupOverall=="ReptilesYSF")
nrow(ReYSF)
median(ReYSF$RRFSp)
R<-10000
B6<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  ReYSF2<-ddply(ReYSF,.(Nstudy))
  boot.sample <- sample(ReYSF2$RRFSp, replace = TRUE)
  B6[kk] <- mean(boot.sample)
}
boxplot(B6)
quantile(B6,c(0.025,0.975))
hist(B6, breaks = 30)
dB6<-c(Bmean=mean(B6),quantile(B6,c(0.025,0.975)),n=nrow(ReYSF2))
dB6
#B7: Reptiles in MSF, overall
ReMSF<-subset(data,BGroupOverall=="ReptilesMSF")
nrow(ReMSF)
median(ReMSF$RRFSp)
R<-10000
B7<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  ReMSF2<-ddply(ReMSF,.(Nstudy))
  boot.sample <- sample(ReMSF2$RRFSp, replace = TRUE)
  B7[kk] <- mean(boot.sample)
}
boxplot(B7)
quantile(B7,c(0.025,0.975))
hist(B7, breaks = 30)
dB7<-c(Bmean=mean(B7),quantile(B7,c(0.025,0.975)),n=nrow(ReMSF2))
dB7
#B8: Reptiles in OSF, overall
ReOSF<-subset(data,BGroupOverall=="ReptilesOSF")
nrow(ReOSF)
B8<-rep(mean(ReOSF$RRFSp),10000)
dB8<-c(Bmean=mean(B8),"2.5%" = NaN,"97.5%" = NaN,n=nrow(ReOSF))
dB8
#B9: Birds in ES, overall
BiES<-subset(data,BGroupOverall=="BirdsES")
nrow(BiES)
median(BiES$RRFSp)
R<-10000
B9<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  BiES2<-ddply(BiES,.(Nstudy))
  boot.sample <- sample(BiES2$RRFSp, replace = TRUE)
  B9[kk] <- mean(boot.sample)
}
boxplot(B9)
quantile(B9,c(0.025,0.975))
hist(B9, breaks = 30)
dB9<-c(Bmean=mean(B9),quantile(B9,c(0.025,0.975)),n=nrow(BiES2))
dB9
#B10: Birds in YSF, overall
BiYSF<-subset(data,BGroupOverall=="BirdsYSF")
median(BiYSF$RRFSp)
nrow(BiYSF)
R<-10000
B10<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  BiYSF2<-ddply(BiYSF,.(Nstudy))
  boot.sample <- sample(BiYSF2$RRFSp, replace = TRUE)
  B10[kk] <- mean(boot.sample)
}
boxplot(B10)
quantile(B10,c(0.025,0.975))
hist(B10, breaks = 30)
dB10<-c(Bmean=mean(B10),quantile(B10,c(0.025,0.975)),n=nrow(BiYSF2))
dB10
#B11: Birds in MSF, overall
BiMSF<-subset(data,BGroupOverall=="BirdsMSF")
median(BiMSF$RRFSp)
nrow(BiMSF)
R<-10000
B11<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  BiMSF2<-ddply(BiMSF,.(Nstudy))
  boot.sample <- sample(BiMSF2$RRFSp, replace = TRUE)
  B11[kk] <- mean(boot.sample)
}
boxplot(B11)
quantile(B11,c(0.025,0.975))
hist(B11, breaks = 30)
dB11<-c(Bmean=mean(B11),quantile(B11,c(0.025,0.975)),n=nrow(BiMSF2))
dB11
#B12: Birds in OSF, overall
BiOSF<-subset(data,BGroupOverall=="BirdsOSF")
nrow(BiOSF)
B12<-rep(mean(BiOSF$RRFSp),10000)
dB12<-c(Bmean=mean(B12),"2.5%" = NaN,"97.5%" = NaN,n=nrow(BiOSF))
dB12
#B13: Mammals in ES, overall
MaES<-subset(data,BGroupOverall=="MammalsES")
median(MaES$RRFSp)
nrow(MaES)
B13<-rep(mean(MaES$RRFSp),10000)
dB13<-c(Bmean=mean(B13),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MaES))
dB13
#B14: Mammals in YSF, overall
MaYSF<-subset(data,BGroupOverall=="MammalsYSF")
median(MaYSF$RRFSp)
nrow(MaYSF)
R<-10000
B14<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MaYSF2<-ddply(MaYSF,.(Nstudy))
  boot.sample <- sample(MaYSF2$RRFSp, replace = TRUE)
  B14[kk] <- mean(boot.sample)
}
boxplot(B14)
quantile(B14,c(0.025,0.975))
hist(B14, breaks = 30)
dB14<-c(Bmean=mean(B14),quantile(B14,c(0.025,0.975)),n=nrow(MaYSF2))
dB14
#B15: Mammals in MSF, overall
MaMSF<-subset(data,BGroupOverall=="MammalsMSF")
median(MaMSF$RRFSp)
nrow(MaMSF)
R<-10000
B15<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MaMSF2<-ddply(MaMSF,.(Nstudy))
  boot.sample <- sample(MaMSF2$RRFSp, replace = TRUE)
  B15[kk] <- mean(boot.sample)
}
boxplot(B15)
quantile(B15,c(0.025,0.975))
hist(B15, breaks = 30)
dB15<-c(Bmean=mean(B15),quantile(B15,c(0.025,0.975)),n=nrow(MaMSF2))
dB15
#B16: Mammals in OSF, overall
MaOSF<-subset(data,BGroupOverall=="MammalsOSF")
nrow(MaOSF)
B16<-rep(mean(MaOSF$RRFSp),10000)
dB16<-c(Bmean=mean(B16),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MaOSF))
dB16

#Join results
#Resume data of Bootstrap (Bootstrap mean, confidence limits[2.5%-97.5%],n)
dbR=data.frame(dB1,dB2,dB3,dB4,dB5,dB6,dB7,dB8,dB9,dB10, dB11,dB12,dB13,dB14,dB15,dB16)
t.dbR<-t(dbR)
datboots<-as.data.frame(t.dbR)
head(datboots)
datboots$FigurePart<-c(rep("FSpecialists",16))
taxdb<-c(rep("Amphibians",4),rep("Reptiles",4),rep("Birds",4),rep("Mammals",4))
Succ.stage<-c(rep("ES",1),rep("YSF",1),rep("MSF",1),rep("OSF",1))
datboots$Taxa<-c(rep(taxdb,1))
datboots$Succ.Stage<-c(rep(Succ.stage,4))
head(datboots)
tail(datboots)
write.table(datboots, "BootstrappFSpecialist.txt",quote=F, sep="\t")


#Lets graph all this work
colnames(datboots)<- c("Bmean","lower","upper", "n","FigurePart","Taxa", "Succ.Stage")
datboots$n

datboots$Succ.Stage <- factor(datboots$Succ.Stage,levels=c('ES','YSF','MSF', 'OSF'))
datboots$Taxa <- factor(datboots$Taxa,levels=c('Mammals','Birds','Reptiles','Amphibians'))

ggplot(datboots, aes(y = Bmean, ymin = lower, ymax = upper,
                 x = Taxa, shape=Succ.Stage, fill=Taxa))+
  geom_vline(xintercept=c(1.5,2.5,3.5),color="darkgray", size=0.4)+
  geom_hline(yintercept = 0, linetype = "dashed",color="black", size=0.5)+
  geom_text(aes(label=datboots$n, y=0.85, fill=Taxa), 
            position = position_dodge(width = 1), size=2.5)+
  #geom_pointrange(position = position_dodge(1.2), size=0.8,aes(fill=Taxa), colour="black", stroke=1)+
  geom_linerange(position = position_dodge(1), size=1.5, alpha=0.6, aes(color=Taxa))+
  geom_point(position = position_dodge(1), size=3.5, stroke = 1, alpha=0.8,  aes(fill=Taxa))+
  ylab("Response ratio of the vertebrate forest specialist species (tropical moist forest)\n  during secondary forest succession(bootstrapped effect size)")+  
  xlab("")+
  scale_shape_manual(name="Succ.Stage",
                     values = c("ES" = 21, "YSF"=22, 
                                "MSF"=24,"OSF"=23))+
  scale_y_continuous(breaks = seq(-3,0.9, by=0.5),
                     labels = seq(-3,0.9, by=0.5),
                     limits = c(-3.25,1.1), expand = c(0, 0))+
  scale_x_discrete(breaks=c('Amphibians','Reptiles','Birds','Mammals'),
                   labels=c('Amphibians','Reptiles','Birds','Mammals'), 
                   expand = c(0.025, 0))+
  scale_fill_manual(name="Taxa",
                    values = c("Amphibians" = "#397d34", "Reptiles"="#FFdd02", 
                               "Birds"="#1f78b4","Mammals"="#FF7f00"))+
  scale_color_manual(name="Taxa",
                     values = c("Amphibians" = "#397d34", "Reptiles"="#FFdd02", 
                               "Birds"="#1f78b4","Mammals"="#FF7f00"))+
  theme_bw()+
  theme(legend.position="none",axis.text.x=element_text(size=12, hjust = 0.5), axis.text.y=element_blank(), 
        plot.margin=unit(c(1,1,1,1),"mm"),panel.margin.y = unit(0, "lines"), panel.grid.major = element_blank(), panel.grid.minor = 	element_blank())+ 
  theme(strip.background = element_blank(),plot.title=element_text(hjust=-0.025, size = 16, face = "bold"))+
  theme(strip.text = element_blank())+
  coord_flip()+
  facet_wrap(~FigurePart, ncol = 3, nrow = 1)

