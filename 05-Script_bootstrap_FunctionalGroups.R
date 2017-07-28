#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Bootstrapping of the Response ratio of Functional groups (Birds and Mammals in Tropical Moist Forest)
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

#Insectivorous
#Exclude the forest with zero richness (the RR can't be calculated)
dataIn<-subset(mydata,RRIn!="NA")
nrow(dataIn)#70
head(dataIn)
#To select one comparison per study (avoid spacial pseudoreplication)
randomRows= function(df,n){
  return(df[sample(nrow(df),n),])
}

#B1: Birds Insectivores in Early Succession (ES)
BiInES<-subset(dataIn,BGroupOverall=="BirdsES")
nrow(BiInES)
median(BiInES$RRIn)
R<-10000
B1<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  BiInES2<-ddply(BiInES,.(Nstudy))
  boot.sample <- sample(BiInES2$RRIn, replace = TRUE)
  B1[kk] <- mean(boot.sample)
}
boxplot(B1)
quantile(B1,c(0.025,0.975))
hist(B1, breaks = 30)
dB1<-c(Bmean=mean(B1),quantile(B1,c(0.025,0.975)),n=nrow(BiInES2))
dB1

#B2: Birds Insectivores in Young Secondary Forest (YSF)
BiInYSF<-subset(dataIn,BGroupOverall=="BirdsYSF")
median(BiInYSF$RRIn)
nrow(BiInYSF)
R<-10000
B2<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  BiInYSF2<-ddply(BiInYSF,.(Nstudy))
  boot.sample <- sample(BiInYSF2$RRIn, replace = TRUE)
  B2[kk] <- mean(boot.sample)
}
boxplot(B2)
quantile(B2,c(0.025,0.975))
hist(B2, breaks = 30)
dB2<-c(Bmean=mean(B2),quantile(B2,c(0.025,0.975)),n=nrow(BiInYSF2))
dB2

#B3: Birds Insectivores in Mid-successional Secondary Forest (MSF)
BiInMSF<-subset(dataIn,BGroupOverall=="BirdsMSF")
median(BiInMSF$RRIn)
nrow(BiInMSF)
R<-10000
B3<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  BiInMSF2<-ddply(BiInMSF,.(Nstudy))
  boot.sample <- sample(BiInMSF2$RRIn, replace = TRUE)
  B3[kk] <- mean(boot.sample)
}
boxplot(B3)
quantile(B3,c(0.025,0.975))
hist(B3, breaks = 30)
dB3<-c(Bmean=mean(B3),quantile(B3,c(0.025,0.975)),n=nrow(BiInMSF2))
dB3

#B4: Birds Insectivores in Old Secondary Forest (OSF)
BiInOSF<-subset(dataIn,BGroupOverall=="AmphibiansOSF")
nrow(BiInOSF)
B4<-rep(mean(BiInOSF$RRIn),10000)
dB4<-c(Bmean=mean(B4),"2.5%" = NaN,"97.5%" = NaN,n=nrow(BiInOSF))
dB4

#B5: Mammals Insectivores in ES, overall
MaInES<-subset(dataIn,BGroupOverall=="MammalsES")
nrow(MaInES)
B5<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MaInES2<-ddply(MaInES,.(Nstudy))
  boot.sample <- sample(MaInES2$RRIn, replace = TRUE)
  B5[kk] <- mean(boot.sample)
}
boxplot(B5)
quantile(B5,c(0.025,0.975))
hist(B5, breaks = 30)
dB5<-c(Bmean=mean(B5),quantile(B5,c(0.025,0.975)),n=nrow(MaInES2))
dB5

#B6: Mammals Insectivores in YSF
MaInYSF<-subset(dataIn,BGroupOverall=="MammalsYSF")
nrow(MaInYSF)
median(MaInYSF$RRIn)
R<-10000
B6<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MaInYSF2<-ddply(MaInYSF,.(Nstudy))
  boot.sample <- sample(MaInYSF2$RRIn, replace = TRUE)
  B6[kk] <- mean(boot.sample)
}
boxplot(B6)
quantile(B6,c(0.025,0.975))
hist(B6, breaks = 30)
dB6<-c(Bmean=mean(B6),quantile(B6,c(0.025,0.975)),n=nrow(MaInYSF2))
dB6

#B7: Mammals Insectivores in MSF
MaInMSF<-subset(dataIn,BGroupOverall=="MammalsMSF")
nrow(MaInMSF)
median(MaInMSF$RRIn)
R<-10000
B7<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MaInMSF2<-ddply(MaInMSF,.(Nstudy))
  boot.sample <- sample(MaInMSF2$RRIn, replace = TRUE)
  B7[kk] <- mean(boot.sample)
}
boxplot(B7)
quantile(B7,c(0.025,0.975))
hist(B7, breaks = 30)
dB7<-c(Bmean=mean(B7),quantile(B7,c(0.025,0.975)),n=nrow(MaInMSF2))
dB7

#B8: Mammals Insectivores in OSF
MaInOSF<-subset(dataIn,BGroupOverall=="MammalsOSF")
nrow(MaInOSF)
B8<-rep(mean(MaInOSF$RRIn),10000)
dB8<-c(Bmean=mean(B8),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MaInOSF))
dB8

#Seed Dispersers
#Exclude the forest with zero richness (the RR can't be calculated)
dataSD<-subset(mydata,RRSD!="NA")
nrow(dataSD)#74
head(dataSD)

#B9: Birds Seed Dispersers in Early Succession (ES)
BiSDES<-subset(dataSD,BGroupOverall=="BirdsES")
nrow(BiSDES)
median(BiSDES$RRSD)
R<-10000
B9<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  BiSDES2<-ddply(BiSDES,.(Nstudy))
  boot.sample <- sample(BiSDES2$RRSD, replace = TRUE)
  B9[kk] <- mean(boot.sample)
}
boxplot(B9)
quantile(B9,c(0.025,0.975))
hist(B9, breaks = 30)
dB9<-c(Bmean=mean(B9),quantile(B9,c(0.025,0.975)),n=nrow(BiSDES2))
dB9

#B10: Birds Seed Dispersers in Young Secondary Forest (YSF)
BiSDYSF<-subset(dataSD,BGroupOverall=="BirdsYSF")
median(BiSDYSF$RRSD)
nrow(BiSDYSF)
R<-10000
B10<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  BiSDYSF2<-ddply(BiSDYSF,.(Nstudy))
  boot.sample <- sample(BiSDYSF2$RRSD, replace = TRUE)
  B10[kk] <- mean(boot.sample)
}
boxplot(B10)
quantile(B10,c(0.025,0.975))
hist(B10, breaks = 30)
dB10<-c(Bmean=mean(B10),quantile(B10,c(0.025,0.975)),n=nrow(BiSDYSF2))
dB10

#B11: Birds Seed Dispersers in Mid-successional Secondary Forest (MSF)
BiSDMSF<-subset(dataSD,BGroupOverall=="BirdsMSF")
median(BiSDMSF$RRSD)
nrow(BiSDMSF)
R<-10000
B11<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  BiSDMSF2<-ddply(BiSDMSF,.(Nstudy))
  boot.sample <- sample(BiSDMSF2$RRSD, replace = TRUE)
  B11[kk] <- mean(boot.sample)
}
boxplot(B11)
quantile(B11,c(0.025,0.975))
hist(B11, breaks = 30)
dB11<-c(Bmean=mean(B11),quantile(B11,c(0.025,0.975)),n=nrow(BiSDMSF2))
dB11

#B12: Birds Seed Dispersers in Old Secondary Forest (OSF)
BiSDOSF<-subset(dataSD,BGroupOverall=="AmphibiansOSF")
nrow(BiSDOSF)
B12<-rep(mean(BiSDOSF$RRSD),10000)
dB12<-c(Bmean=mean(B12),"2.5%" = NaN,"97.5%" = NaN,n=nrow(BiSDOSF))
dB12

#B13: Mammals Seed Dispersers in ES, overall
MaSDES<-subset(dataSD,BGroupOverall=="MammalsES")
nrow(MaSDES)
B13<-rep(mean(MaSDES $RRSD),10000)
dB13<-c(Bmean=mean(B13),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MaSDES))
dB13

#B14: Mammals Seed Dispersers in YSF
MaSDYSF<-subset(dataSD,BGroupOverall=="MammalsYSF")
nrow(MaSDYSF)
median(MaSDYSF$RRSD)
R<-10000
B14<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MaSDYSF2<-ddply(MaSDYSF,.(Nstudy))
  boot.sample <- sample(MaSDYSF2$RRSD, replace = TRUE)
  B14[kk] <- mean(boot.sample)
}
boxplot(B14)
quantile(B14,c(0.025,0.975))
hist(B14, breaks = 30)
dB14<-c(Bmean=mean(B14),quantile(B14,c(0.025,0.975)),n=nrow(MaSDYSF2))
dB14

#B15: Mammals Seed Dispersers in MSF
MaSDMSF<-subset(dataSD,BGroupOverall=="MammalsMSF")
nrow(MaSDMSF)
median(MaSDMSF$RRSD)
R<-10000
B15<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MaSDMSF2<-ddply(MaSDMSF,.(Nstudy))
  boot.sample <- sample(MaSDMSF2$RRSD, replace = TRUE)
  B15[kk] <- mean(boot.sample)
}
boxplot(B15)
quantile(B15,c(0.025,0.975))
hist(B15, breaks = 30)
dB15<-c(Bmean=mean(B15),quantile(B15,c(0.025,0.975)),n=nrow(MaSDMSF2))
dB15

#B16: Mammals Seed Dispersers in OSF
MaSDOSF<-subset(dataSD,BGroupOverall=="MammalsOSF")
nrow(MaSDOSF)
B16<-rep(mean(MaSDOSF$RRSD),10000)
dB16<-c(Bmean=mean(B16),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MaSDOSF))
dB16

#Granivores
#Exclude the forest with zero richness (the RR can't be calculated)
dataGr<-subset(mydata,RRGr!="NA")
nrow(dataGr)#57
head(dataGr)

#B17: Birds Granivores in Early Succession (ES)
BiGrES<-subset(dataGr,BGroupOverall=="BirdsES")
nrow(BiGrES)
median(BiGrES$RRGr)
R<-10000
B17<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  BiGrES2<-ddply(BiGrES,.(Nstudy))
  boot.sample <- sample(BiGrES2$RRGr, replace = TRUE)
  B17[kk] <- mean(boot.sample)
}
boxplot(B17)
quantile(B17,c(0.025,0.975))
hist(B17, breaks = 30)
dB17<-c(Bmean=mean(B17),quantile(B17,c(0.025,0.975)),n=nrow(BiGrES2))
dB17

#B18: Birds Granivores in Young Secondary Forest (YSF)
BiGrYSF<-subset(dataGr,BGroupOverall=="BirdsYSF")
median(BiGrYSF$RRGr)
nrow(BiGrYSF)
R<-10000
B18<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  BiGrYSF2<-ddply(BiGrYSF,.(Nstudy))
  boot.sample <- sample(BiGrYSF2$RRGr, replace = TRUE)
  B18[kk] <- mean(boot.sample)
}
boxplot(B18)
quantile(B18,c(0.025,0.975))
hist(B18, breaks = 30)
dB18<-c(Bmean=mean(B18),quantile(B18,c(0.025,0.975)),n=nrow(BiGrYSF2))
dB18

#B19: Birds Granivores in Mid-successional Secondary Forest (MSF)
BiGrMSF<-subset(dataGr,BGroupOverall=="BirdsMSF")
median(BiGrMSF$RRGr)
nrow(BiGrMSF)
R<-10000
B19<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  BiGrMSF2<-ddply(BiGrMSF,.(Nstudy))
  boot.sample <- sample(BiGrMSF2$RRGr, replace = TRUE)
  B19[kk] <- mean(boot.sample)
}
boxplot(B19)
quantile(B19,c(0.025,0.975))
hist(B19, breaks = 30)
dB19<-c(Bmean=mean(B19),quantile(B19,c(0.025,0.975)),n=nrow(BiGrMSF2))
dB19

#B20: Birds Granivores in Old Secondary Forest (OSF)
BiGrOSF<-subset(dataGr,BGroupOverall=="AmphibiansOSF")
nrow(BiGrOSF)
B20<-rep(mean(BiGrOSF$RRGr),10000)
dB20<-c(Bmean=mean(B20),"2.5%" = NaN,"97.5%" = NaN,n=nrow(BiGrOSF))
dB20

#B21: Mammals Granivores in ES, overall
MaGrES<-subset(dataGr,BGroupOverall=="MammalsES")
nrow(MaGrES)
B21<-rep(mean(MaGrES $RRGr),10000)
dB21<-c(Bmean=mean(B21),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MaGrES))
dB21

#B22: Mammals Granivores in YSF
MaGrYSF<-subset(dataGr,BGroupOverall=="MammalsYSF")
nrow(MaGrYSF)
median(MaGrYSF$RRGr)
R<-10000
B22<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MaGrYSF2<-ddply(MaGrYSF,.(Nstudy))
  boot.sample <- sample(MaGrYSF2$RRGr, replace = TRUE)
  B22[kk] <- mean(boot.sample)
}
boxplot(B22)
quantile(B22,c(0.025,0.975))
hist(B22, breaks = 30)
dB22<-c(Bmean=mean(B22),quantile(B22,c(0.025,0.975)),n=nrow(MaGrYSF2))
dB22

#B23: Mammals Granivores in MSF
MaGrMSF<-subset(dataGr,BGroupOverall=="MammalsMSF")
nrow(MaGrMSF)
B23<-rep(mean(MaGrMSF$RRGr),10000)
dB23<-c(Bmean=mean(B23),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MaGrMSF))
dB23

#B24: Mammals Granivores in OSF
MaGrOSF<-subset(dataGr,BGroupOverall=="MammalsOSF")
nrow(MaGrOSF)
B24<-rep(mean(MaGrOSF$RRGr),10000)
dB24<-c(Bmean=mean(B24),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MaGrOSF))
dB24


#Predators
#Exclude the forest with zero richness (the RR can't be calculated)
dataPr<-subset(mydata,RRPr!="NA")
nrow(dataPr)#42
head(dataPr)

#B25: Birds Predators in Early Succession (ES)
BiPrES<-subset(dataPr,BGroupOverall=="BirdsES")
nrow(BiPrES)
median(BiPrES$RRPr)
R<-10000
B25<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  BiPrES2<-ddply(BiPrES,.(Nstudy))
  boot.sample <- sample(BiPrES2$RRPr, replace = TRUE)
  B25[kk] <- mean(boot.sample)
}
boxplot(B25)
quantile(B25,c(0.025,0.975))
hist(B25, breaks = 30)
dB25<-c(Bmean=mean(B25),quantile(B25,c(0.025,0.975)),n=nrow(BiPrES2))
dB25

#B26: Birds Predators in Young Secondary Forest (YSF)
BiPrYSF<-subset(dataPr,BGroupOverall=="BirdsYSF")
median(BiPrYSF$RRPr)
nrow(BiPrYSF)
R<-10000
B26<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  BiPrYSF2<-ddply(BiPrYSF,.(Nstudy))
  boot.sample <- sample(BiPrYSF2$RRPr, replace = TRUE)
  B26[kk] <- mean(boot.sample)
}
boxplot(B26)
quantile(B26,c(0.025,0.975))
hist(B26, breaks = 30)
dB26<-c(Bmean=mean(B26),quantile(B26,c(0.025,0.975)),n=nrow(BiPrYSF2))
dB26

#B27: Birds Predators in Mid-successional Secondary Forest (MSF)
BiPrMSF<-subset(dataPr,BGroupOverall=="BirdsMSF")
median(BiPrMSF$RRPr)
nrow(BiPrMSF)
R<-10000
B27<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  BiPrMSF2<-ddply(BiPrMSF,.(Nstudy))
  boot.sample <- sample(BiPrMSF2$RRPr, replace = TRUE)
  B27[kk] <- mean(boot.sample)
}
boxplot(B27)
quantile(B27,c(0.025,0.975))
hist(B27, breaks = 30)
dB27<-c(Bmean=mean(B27),quantile(B27,c(0.025,0.975)),n=nrow(BiPrMSF2))
dB27

#B28: Birds Predators in Old Secondary Forest (OSF)
BiPrOSF<-subset(dataPr,BGroupOverall=="AmphibiansOSF")
nrow(BiPrOSF)
B28<-rep(mean(BiPrOSF$RRPr),10000)
dB28<-c(Bmean=mean(B28),"2.5%" = NaN,"97.5%" = NaN,n=nrow(BiPrOSF))
dB28

#B29: Mammals Predators in ES, overall
MaPrES<-subset(dataPr,BGroupOverall=="MammalsES")
nrow(MaPrES)
B29<-rep(mean(MaPrES $RRPr),10000)
dB29<-c(Bmean=mean(B29),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MaPrES))
dB29

#B30: Mammals Predators in YSF
MaPrYSF<-subset(dataPr,BGroupOverall=="MammalsYSF")
nrow(MaPrYSF)
median(MaPrYSF$RRPr)
R<-10000
B30<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MaPrYSF2<-ddply(MaPrYSF,.(Nstudy))
  boot.sample <- sample(MaPrYSF2$RRPr, replace = TRUE)
  B30[kk] <- mean(boot.sample)
}
boxplot(B30)
quantile(B30,c(0.025,0.975))
hist(B30, breaks = 30)
dB30<-c(Bmean=mean(B30),quantile(B30,c(0.025,0.975)),n=nrow(MaPrYSF2))
dB30

#B31: Mammals Predators in MSF
MaPrMSF<-subset(dataPr,BGroupOverall=="MammalsMSF")
nrow(MaPrMSF)
B31<-rep(mean(MaPrMSF$RRPr),10000)
dB31<-c(Bmean=mean(B31),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MaPrMSF))
dB31

#B32: Mammals Predators in OSF
MaPrOSF<-subset(dataPr,BGroupOverall=="MammalsOSF")
nrow(MaPrOSF)
B32<-rep(mean(MaPrOSF$RRPr),10000)
dB32<-c(Bmean=mean(B32),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MaPrOSF))
dB32

#Pollinators
#Exclude the forest with zero richness (the RR can't be calculated)
dataPo<-subset(mydata,RRPo!="NA")
nrow(dataPo)#32
head(dataPo)

#B33: Birds Pollinators in Early Succession (ES)
BiPoES<-subset(dataPo,BGroupOverall=="BirdsES")
nrow(BiPoES)
B33<-rep(mean(BiPoES$RRPo),10000)
dB33<-c(Bmean=mean(B33),"2.5%" = NaN,"97.5%" = NaN,n=nrow(BiPoES))
dB33

#B34: Birds Pollinators in Young Secondary Forest (YSF)
BiPoYSF<-subset(dataPo,BGroupOverall=="BirdsYSF")
median(BiPoYSF$RRPo)
nrow(BiPoYSF)
R<-10000
B34<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  BiPoYSF2<-ddply(BiPoYSF,.(Nstudy))
  boot.sample <- sample(BiPoYSF2$RRPo, replace = TRUE)
  B34[kk] <- mean(boot.sample)
}
boxplot(B34)
quantile(B34,c(0.025,0.975))
hist(B34, breaks = 30)
dB34<-c(Bmean=mean(B34),quantile(B34,c(0.025,0.975)),n=nrow(BiPoYSF))
dB34

#B35: Birds Pollinators in Mid-successional Secondary Forest (MSF)
BiPoMSF<-subset(dataPo,BGroupOverall=="BirdsMSF")
median(BiPoMSF$RRPo)
nrow(BiPoMSF)
R<-10000
B35<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  BiPoMSF2<-ddply(BiPoMSF,.(Nstudy))
  boot.sample <- sample(BiPoMSF2$RRPo, replace = TRUE)
  B35[kk] <- mean(boot.sample)
}
boxplot(B35)
quantile(B35,c(0.025,0.975))
hist(B35, breaks = 30)
dB35<-c(Bmean=mean(B35),quantile(B35,c(0.025,0.975)),n=nrow(BiPoMSF))
dB35

#B36: Birds Pollinators in Old Secondary Forest (OSF)
BiPoOSF<-subset(dataPo,BGroupOverall=="AmphibiansOSF")
nrow(BiPoOSF)
B36<-rep(mean(BiPoOSF$RRPo),10000)
dB36<-c(Bmean=mean(B36),"2.5%" = NaN,"97.5%" = NaN,n=nrow(BiPoOSF))
dB36

#B37: Mammals Pollinators in ES, overall
MaPoES<-subset(dataPo,BGroupOverall=="MammalsES")
nrow(MaPoES)
B37<-rep(mean(MaPoES $RRPo),10000)
dB37<-c(Bmean=mean(B37),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MaPoES))
dB37

#B38: Mammals Pollinators in YSF
MaPoYSF<-subset(dataPo,BGroupOverall=="MammalsYSF")
nrow(MaPoYSF)
median(MaPoYSF$RRPo)
R<-10000
B38<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MaPoYSF2<-ddply(MaPoYSF,.(Nstudy))
  boot.sample <- sample(MaPoYSF2$RRPo, replace = TRUE)
  B38[kk] <- mean(boot.sample)
}
boxplot(B38)
quantile(B38,c(0.025,0.975))
hist(B38, breaks = 30)
dB38<-c(Bmean=mean(B38),quantile(B38,c(0.025,0.975)),n=nrow(MaPoYSF2))
dB38

#B39: Mammals Pollinators in MSF
MaPoMSF<-subset(dataPo,BGroupOverall=="MammalsMSF")
nrow(MaPoMSF)
B39<-rep(mean(MaPoMSF$RRPo),10000)
dB39<-c(Bmean=mean(B39),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MaPoMSF))
dB39

#B40: Mammals Pollinators in OSF
MaPoOSF<-subset(dataPo,BGroupOverall=="MammalsOSF")
nrow(MaPoOSF)
B40<-rep(mean(MaPoOSF$RRPo),10000)
dB40<-c(Bmean=mean(B40),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MaPoOSF))
dB40

#Scavengers
#Exclude the forest with zero richness (the RR can't be calculated)
dataSc<-subset(mydata,RRSc!="NA")
nrow(dataSc)#14
head(dataSc)

#B41: Birds Scavengers in Early Succession (ES)
BiScES<-subset(dataSc,BGroupOverall=="BirdsES")
nrow(BiScES)
B41<-rep(mean(BiScES$RRSc),10000)
dB41<-c(Bmean=mean(B41),"2.5%" = NaN,"97.5%" = NaN,n=nrow(BiScES))
dB41

#B42: Birds Scavengers in Young Secondary Forest (YSF)
BiScYSF<-subset(dataSc,BGroupOverall=="BirdsYSF")
median(BiScYSF$RRSc)
nrow(BiScYSF)
R<-10000
B42<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  BiScYSF2<-ddply(BiScYSF,.(Nstudy))
  boot.sample <- sample(BiScYSF2$RRSc, replace = TRUE)
  B42[kk] <- mean(boot.sample)
}
boxplot(B42)
quantile(B42,c(0.025,0.975))
hist(B42, breaks = 30)
dB42<-c(Bmean=mean(B42),quantile(B42,c(0.025,0.975)),n=nrow(BiScYSF2))
dB42

#B43: Birds Scavengers in Mid-successional Secondary Forest (MSF)
BiScMSF<-subset(dataSc,BGroupOverall=="BirdsMSF")
median(BiScMSF$RRSc)
nrow(BiScMSF)
B43<-rep(mean(BiScMSF$RRSc),10000)
dB43<-c(Bmean=mean(B43),"2.5%" = NaN,"97.5%" = NaN,n=nrow(BiScMSF))
dB43

#B44: Birds Scavengers in Old Secondary Forest (OSF)
BiScOSF<-subset(dataSc,BGroupOverall=="AmphibiansOSF")
nrow(BiScOSF)
B44<-rep(mean(BiScOSF$RRSc),10000)
dB44<-c(Bmean=mean(B44),"2.5%" = NaN,"97.5%" = NaN,n=nrow(BiScOSF))
dB44

#B45: Mammals Scavengers in ES, overall
MaScES<-subset(dataSc,BGroupOverall=="MammalsES")
nrow(MaScES)
B45<-rep(mean(MaScES $RRSc),10000)
dB45<-c(Bmean=mean(B45),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MaScES))
dB45

#B46: Mammals Scavengers in YSF
MaScYSF<-subset(dataSc,BGroupOverall=="MammalsYSF")
nrow(MaScYSF)
B46<-rep(mean(MaScYSF$RRSc),10000)
dB46<-c(Bmean=mean(B46),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MaScYSF))
dB46

#B47: Mammals Scavengers in MSF
MaScMSF<-subset(dataSc,BGroupOverall=="MammalsMSF")
nrow(MaScMSF)
B47<-rep(mean(MaScMSF$RRSc),10000)
dB47<-c(Bmean=mean(B47),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MaScMSF))
dB47

#B48: Mammals Scavengers in OSF
MaScOSF<-subset(dataSc,BGroupOverall=="MammalsOSF")
nrow(MaScOSF)
B48<-rep(mean(MaScOSF$RRSc),10000)
dB48<-c(Bmean=mean(B48),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MaScOSF))
dB48
#Grazers
#Exclude the forest with zero richness (the RR can't be calculated)
dataGz<-subset(mydata,RRGz!="NA")
nrow(dataGz)#21
head(dataGz)

#B49: Birds Grazers in Early Succession (ES)
BiGzES<-subset(dataGz,BGroupOverall=="BirdsES")
nrow(BiGzES)
B49<-rep(mean(BiGzES$RRSc),10000)
dB49<-c(Bmean=mean(B49),"2.5%" = NaN,"97.5%" = NaN,n=nrow(BiGzES))
dB49

#B50: Birds Grazers in Young Secondary Forest (YSF)
BiGzYSF<-subset(dataGz,BGroupOverall=="BirdsYSF")
median(BiGzYSF$RRGz)
nrow(BiGzYSF)
B50<-rep(mean(BiGzYSF$RRSc),10000)
dB50<-c(Bmean=mean(B50),"2.5%" = NaN,"97.5%" = NaN,n=nrow(BiGzYSF))
dB50

#B51: Birds Grazers in Mid-successional Secondary Forest (MSF)
BiGzMSF<-subset(dataGz,BGroupOverall=="BirdsMSF")
median(BiGzMSF$RRGz)
nrow(BiGzMSF)
B51<-rep(mean(BiGzMSF$RRSc),10000)
dB51<-c(Bmean=mean(B51),"2.5%" = NaN,"97.5%" = NaN,n=nrow(BiGzMSF))
dB51

#B52: Birds Grazers in Old Secondary Forest (OSF)
BiGzOSF<-subset(dataGz,BGroupOverall=="AmphibiansOSF")
nrow(BiGzOSF)
B52<-rep(mean(BiGzOSF$RRGz),10000)
dB52<-c(Bmean=mean(B52),"2.5%" = NaN,"97.5%" = NaN,n=nrow(BiGzOSF))
dB52

#B53: Mammals Grazers in ES, overall
MaGzES<-subset(dataGz,BGroupOverall=="MammalsES")
nrow(MaGzES)
B53<-rep(mean(MaGzES$RRGz),10000)
dB53<-c(Bmean=mean(B53),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MaGzES))
dB53

#B54: Mammals Grazers in YSF
MaGzYSF<-subset(dataGz,BGroupOverall=="MammalsYSF")
nrow(MaGzYSF)
median(MaGzYSF$RRGz)
R<-10000
B54<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MaGzYSF2<-ddply(MaGzYSF,.(Nstudy))
  boot.sample <- sample(MaGzYSF2$RRGz, replace = TRUE)
  B54[kk] <- mean(boot.sample)
}
boxplot(B54)
quantile(B54,c(0.025,0.975))
hist(B54, breaks = 30)
dB54<-c(Bmean=mean(B54),quantile(B54,c(0.025,0.975)),n=nrow(MaGzYSF2))
dB54

#B55: Mammals Grazers in MSF
MaGzMSF<-subset(dataGz,BGroupOverall=="MammalsMSF")
nrow(MaGzMSF)
B55<-rep(mean(MaGzMSF$RRGz),10000)
dB55<-c(Bmean=mean(B55),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MaGzMSF))
dB55

#B56: Mammals Grazers in OSF
MaGzOSF<-subset(dataGz,BGroupOverall=="MammalsOSF")
nrow(MaGzOSF)
B56<-rep(mean(MaGzOSF$RRGz),10000)
dB56<-c(Bmean=mean(B56),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MaGzOSF))
dB56


#Join results
#Resume data of Bootstrap (Bootstrap mean, confidence limits[2.5%-97.5%],n)
dbR=data.frame(dB1,dB2,dB3,dB4, dB5,dB6,dB7,dB8,
				dB9,dB10,dB11,dB12, dB13,dB14,dB15,dB16,
				dB17,dB18,dB19,dB20, dB21,dB22,dB23,dB24,
				dB25,dB26,dB27,dB28, dB29,dB30,dB31,dB32,
				dB33,dB34,dB35,dB36, dB37,dB38,dB39,dB40,
				dB41,dB42,dB43,dB44, dB45,dB46,dB47,dB48,
				dB49,dB50,dB51,dB52, dB53,dB54,dB55,dB56)
t.dbR<-t(dbR)
datboots<-as.data.frame(t.dbR)
head(datboots)
datboots$FigurePart<-c(rep("In",8),rep("SD",8),rep("Gr",8),rep("Pr",8),rep("Po",8),rep("Sc",8),rep("Gz",8))
taxdb<-c(rep("Birds",4),rep("Mammals",4))
Succ.stage<-c(rep("ES",1),rep("YSF",1),rep("MSF",1),rep("OSF",1))
datboots$Taxa<-c(rep(taxdb,7))
datboots$Succ.Stage<-c(rep(Succ.stage,14))
head(datboots)
tail(datboots)
write.table(datboots, "BootstrappFunctionalGroups.txt",quote=F, sep="\t")


#Lets graph all this work
colnames(datboots)<- c("Bmean","lower","upper", "n","FigurePart","Taxa", "Succ.Stage")
datboots$n

birds<-subset(datboots,Taxa=="Birds")

birds$FigurePart <- factor(birds$FigurePart,levels=c('Gz','Sc','Po','Pr','Gr','SD','In'))
birds$Succ.Stage <- factor(birds$Succ.Stage,levels=c('ES','YSF','MSF', 'OSF'))
Fb<-ggplot(birds,aes(y=Bmean,ymin=lower,ymax=upper,x=FigurePart,shape=Succ.Stage,fill=Taxa))+
 geom_text(aes(label=birds$n, y=0.85, fill=Taxa,shape=Succ.Stage,x=FigurePart), 
            position = position_dodge(width = 1), size=2.5)+
   geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5,6.5),color="darkgray", size=0.4)+
  ggtitle('a')+
  geom_hline(yintercept = 0, linetype = "dashed",color="black", size=0.5)+
  geom_linerange(position = position_dodge(1), size=1.5, alpha=0.6, aes(color=Taxa))+
  geom_point(position = position_dodge(1), size=3.5, stroke = 1, alpha=0.8,  aes(fill=Taxa))+
  scale_x_discrete(breaks=c('In','SD','Gr','Pr','Po','Sc','Gz'),
                   labels=c('Insectivores','Seed Dispersers','Granivores','Predators','Pollinators','Scavengers','Grazers'), 
                   expand = c(0, 0))+
  ylab("")+  
  xlab("")+
  scale_shape_manual(name="Succ.Stage",
                     values = c("ES" = 21, "YSF"=22, 
                                "MSF"=24,"OSF"=23))+
  scale_y_continuous(breaks = seq(-1.5,0.9, by=0.5),
                     labels = seq(-1.5,0.9, by=0.5),
                     limits = c(-1.75,1.1), expand = c(0, 0))+
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
  theme(strip.background = element_blank(),plot.title=element_text(hjust=-0.025, size = 16, face = "bold"))+
  theme(strip.text = element_blank())
Fb
#And mammals
mam<-subset(datboots,Taxa=="Mammals")
mam$FigurePart <- factor(mam$FigurePart,levels=c('Gz','Sc','Po','Pr','Gr','SD','In'))
mam$Succ.Stage <- factor(mam$Succ.Stage,levels=c('ES','YSF','MSF', 'OSF'))
Fm<-ggplot(mam,aes(y=Bmean,ymin=lower,ymax=upper,x=FigurePart,shape=Succ.Stage,fill=Taxa))+
 geom_text(aes(label=mam$n, y=0.85, fill=Taxa,shape=Succ.Stage,x=FigurePart), 
            position = position_dodge(width = 1), size=2.5)+
   geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5,6.5),color="darkgray", size=0.4)+
  ggtitle('b')+
  geom_hline(yintercept = 0, linetype = "dashed",color="black", size=0.5)+
  geom_linerange(position = position_dodge(1), size=1.5, alpha=0.6, aes(color=Taxa))+
  geom_point(position = position_dodge(1), size=3.5, stroke = 1, alpha=0.8,  aes(fill=Taxa))+
  scale_x_discrete(breaks=c('In','SD','Gr','Pr','Po','Sc','Gz'),
                   labels=c('','','','','','',''), 
                   expand = c(0, 0))+
  ylab("")+  
  xlab("")+
  scale_shape_manual(name="Succ.Stage",
                     values = c("ES" = 21, "YSF"=22, 
                                "MSF"=24,"OSF"=23))+
  scale_y_continuous(breaks = seq(-1.5,0.9, by=0.5),
                     labels = seq(-1.5,0.9, by=0.5),
                     limits = c(-1.75,1.1), expand = c(0, 0))+
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
  theme(strip.background = element_blank(),plot.title=element_text(hjust=-0.025, size = 16, face = "bold"))+
  theme(strip.text = element_blank())
Fm

#Joint figures
FSp=grid.arrange(arrangeGrob(Fb,Fm,ncol=2,widths = c(1,0.85)),
                  left=textGrob("Functional groups (dietary preference)", 
                                rot=90, vjust=2),
                  top=textGrob(""),
                  right=textGrob(""),
                  bottom=textGrob("Response ratio during secondary forest succession \n (bootstrapped effect size)", 
                                  gp=gpar(fontsize=14)))



