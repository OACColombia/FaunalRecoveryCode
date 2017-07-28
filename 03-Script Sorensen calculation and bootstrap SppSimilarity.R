#Calculating of Sorensen similarity index for studies with species list (Table S2)
library(vegan)
#study<- read.csv("/Users/orlando/Dropbox/Vertebrate Recovery - Nature/data/Sorensen/6A.csv", row.names=1)
#study<- read.csv("C://Users/Flaco/Dropbox/Vertebrate Recovery - Nature/data/Sorensen/143M.csv", row.names=1)
Sor <- vegdist(study, binary=TRUE, diag = TRUE, upper=FALSE)#when study is each study with species list
Sor#dissimilarity of Sorensen
#all value was included in the MainData.csv matrix

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Bootstrapping of the Response ratio in Species composition similarity
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

#Exclude the forest with zero richness (the RRS can't be calculated)
data<-subset(mydata,RRS!="NA")
nrow(data)#186
head(data)
summary(data$BGroupOverall)

#To select one comparison per study (avoid spacial pseudoreplication)
randomRows= function(df,n){
  return(df[sample(nrow(df),n),])
}
#B1: Amphibians in Early Succession (ES), overall
AmES<-subset(data,BGroupOverall=="AmphibiansES")
nrow(AmES)
median(AmES$RRS)
R<-10000
B1<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  AmES2<-ddply(AmES,.(Nstudy))
    boot.sample <- sample(AmES$RRS, replace = TRUE)
  B1[kk] <- mean(boot.sample)
}
boxplot(B1)
quantile(B1,c(0.025,0.975))
hist(B1, breaks = 30)
dB1<-c(Bmean=mean(B1),quantile(B1,c(0.025,0.975)),n=nrow(AmES))
dB1
#B2: Amphibians in Young Secondary Forest (YSF), overall
AmYSF<-subset(data,BGroupOverall=="AmphibiansYSF")
median(AmYSF$RRS)
nrow(AmYSF)
R<-10000
B2<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  AmYSF2<-ddply(AmYSF,.(Nstudy))
  boot.sample <- sample(AmYSF2$RRS, replace = TRUE)
  B2[kk] <- mean(boot.sample)
}
boxplot(B2)
quantile(B2,c(0.025,0.975))
hist(B2, breaks = 30)
dB2<-c(Bmean=mean(B2),quantile(B2,c(0.025,0.975)),n=nrow(AmYSF2))
dB2
#B3: Amphibians in Mid-successional Secondary Forest (MSF), overall
AmMSF<-subset(data,BGroupOverall=="AmphibiansMSF")
median(AmMSF$RRS)
nrow(AmMSF)
R<-10000
B3<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  AmMSF2<-ddply(AmMSF,.(Nstudy))
  boot.sample <- sample(AmMSF2$RRS, replace = TRUE)
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
B4<-rep(mean(AmOSF$RRS),10000)
dB4<-c(Bmean=mean(B4),"2.5%" = NaN,"97.5%" = NaN,n=nrow(AmOSF))
dB4
#B5: Reptiles in ES, overall
ReES<-subset(data,BGroupOverall=="ReptilesES")
nrow(ReES)
median(ReES$RRS)
R<-10000
B5<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  ReES2<-ddply(ReES,.(Nstudy))
  boot.sample <- sample(ReES2$RRS, replace = TRUE)
  B5[kk] <- mean(boot.sample)
}
boxplot(B5)
quantile(B5,c(0.025,0.975))
hist(B5, breaks = 30)
dB5<-c(Bmean=mean(B5),quantile(B5,c(0.025,0.975)),n=nrow(ReES2))
dB5
#B6: Reptiles in YSF, overall
ReYSF<-subset(data,BGroupOverall=="ReptilesYSF")
nrow(ReYSF)
median(ReYSF$RRS)
R<-10000
B6<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  ReYSF2<-ddply(ReYSF,.(Nstudy))
  boot.sample <- sample(ReYSF2$RRS, replace = TRUE)
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
median(ReMSF$RRS)
R<-10000
B7<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  ReMSF2<-ddply(ReMSF,.(Nstudy))
  boot.sample <- sample(ReMSF2$RRS, replace = TRUE)
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
median(ReOSF$RRS)
nrow(ReOSF)
B8<-rep(mean(ReOSF$RRS),10000)
dB8<-c(Bmean=mean(B8),"2.5%" = NaN,"97.5%" = NaN,n=nrow(ReOSF))
dB8
#B9: Birds in ES, overall
BiES<-subset(data,BGroupOverall=="BirdsES")
nrow(BiES)
median(BiES$RRS)
R<-10000
B9<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  BiES2<-ddply(BiES,.(Nstudy))
  boot.sample <- sample(BiES2$RRS, replace = TRUE)
  B9[kk] <- mean(boot.sample)
}
boxplot(B9)
quantile(B9,c(0.025,0.975))
hist(B9, breaks = 30)
dB9<-c(Bmean=mean(B9),quantile(B9,c(0.025,0.975)),n=nrow(BiES2))
dB9
#B10: Birds in YSF, overall
BiYSF<-subset(data,BGroupOverall=="BirdsYSF")
median(BiYSF$RRS)
nrow(BiYSF)
R<-10000
B10<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  BiYSF2<-ddply(BiYSF,.(Nstudy))
  boot.sample <- sample(BiYSF2$RRS, replace = TRUE)
  B10[kk] <- mean(boot.sample)
}
boxplot(B10)
quantile(B10,c(0.025,0.975))
hist(B10, breaks = 30)
dB10<-c(Bmean=mean(B10),quantile(B10,c(0.025,0.975)),n=nrow(BiYSF2))
dB10
#B11: Birds in MSF, overall
BiMSF<-subset(data,BGroupOverall=="BirdsMSF")
median(BiMSF$RRS)
nrow(BiMSF)
R<-10000
B11<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  BiMSF2<-ddply(BiMSF,.(Nstudy))
  boot.sample <- sample(BiMSF2$RRS, replace = TRUE)
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
median(BiOSF$RRS)
nrow(BiOSF)
B12<-rep(mean(BiOSF$RRS),10000)
dB12<-c(Bmean=mean(B12),"2.5%" = NaN,"97.5%" = NaN,n=nrow(BiOSF))
dB12
#B13: Mammals in ES, overall
MaES<-subset(data,BGroupOverall=="MammalsES")
median(MaES$RRS)
nrow(MaES)
R<-10000
B13<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MaES2<-ddply(MaES,.(Nstudy))
  boot.sample <- sample(MaES2$RRS, replace = TRUE)
  B13[kk] <- mean(boot.sample)
}
boxplot(B13)
quantile(B13,c(0.025,0.975))
hist(B13, breaks = 30)
dB13<-c(Bmean=mean(B13),quantile(B13,c(0.025,0.975)),n=nrow(MaES))
dB13
#B14: Mammals in YSF, overall
MaYSF<-subset(data,BGroupOverall=="MammalsYSF")
median(MaYSF$RRS)
nrow(MaYSF)
R<-10000
B14<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MaYSF2<-ddply(MaYSF,.(Nstudy))
  boot.sample <- sample(MaYSF2$RRS, replace = TRUE)
  B14[kk] <- mean(boot.sample)
}
boxplot(B14)
quantile(B14,c(0.025,0.975))
hist(B14, breaks = 30)
dB14<-c(Bmean=mean(B14),quantile(B14,c(0.025,0.975)),n=nrow(MaYSF2))
dB14
#B15: Mammals in MSF, overall
MaMSF<-subset(data,BGroupOverall=="MammalsMSF")
median(MaMSF$RRS)
nrow(MaMSF)
R<-10000
B15<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MaMSF2<-ddply(MaMSF,.(Nstudy))
  boot.sample <- sample(MaMSF2$RRS, replace = TRUE)
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
B16<-rep(mean(MaOSF$RRS),10000)
dB16<-c(Bmean=mean(B16),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MaOSF))
dB16
#BIOMES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#B17: Amphibians in Early Succession (ES), Tropical Moist Forest
MoAmES<-subset(data,BGroupBiome=="MoistAmphibiansES")
median(MoAmES$RRS)
R<-10000
B17<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MoAmES2<-ddply(MoAmES,.(Nstudy))
  boot.sample <- sample(MoAmES2$RRS, replace = TRUE)
  B17[kk] <- mean(boot.sample)
}
boxplot(B17)
quantile(B17,c(0.025,0.975))
hist(B17, breaks = 30)
dB17<-c(Bmean=mean(B17),quantile(B17,c(0.025,0.975)),n=nrow(MoAmES2))
dB17
#B18: Amphibians in Young Secondary Forest (YSF), Tropical Moist Forest
MoAmYSF<-subset(data,BGroupBiome=="MoistAmphibiansYSF")
median(MoAmYSF$RRS)
nrow(MoAmYSF)
R<-10000
B18<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MoAmYSF2<-ddply(MoAmYSF,.(Nstudy))
  boot.sample <- sample(MoAmYSF2$RRS, replace = TRUE)
  B18[kk] <- mean(boot.sample)
}
boxplot(B18)
quantile(B18,c(0.025,0.975))
hist(B18, breaks = 30)
dB18<-c(Bmean=mean(B18),quantile(B18,c(0.025,0.975)),n=nrow(MoAmYSF2))
dB18
#B19: Amphibians in Mid-successional Secondary Forest (MSF), Tropical Moist Forest
MoAmMSF<-subset(data,BGroupBiome=="MoistAmphibiansMSF")
median(MoAmMSF$RRS)
R<-10000
B19<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MoAmMSF2<-ddply(MoAmMSF,.(Nstudy))
  boot.sample <- sample(MoAmMSF2$RRS, replace = TRUE)
  B19[kk] <- mean(boot.sample)
}
boxplot(B19)
quantile(B19,c(0.025,0.975))
hist(B19, breaks = 30)
dB19<-c(Bmean=mean(B19),quantile(B19,c(0.025,0.975)),n=nrow(MoAmMSF2))
dB19
#B20: Amphibians in Old Secondary Forest (OSF), overall
MoAmOSF<-subset(data,BGroupBiome=="MoistAmphibiansOSF")
nrow(MoAmOSF)
B20<-rep(mean(MoAmOSF$RRS),10000)
dB20<-c(Bmean=mean(B20),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MoAmOSF))
dB20
#B21: Reptiles in ES, Tropical Moist Forest
MoReES<-subset(data,BGroupBiome=="MoistReptilesES")
nrow(MoReES)
R<-10000
B21<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MoReES2<-ddply(MoReES,.(Nstudy))
  boot.sample <- sample(MoReES2$RRS, replace = TRUE)
  B21[kk] <- mean(boot.sample)
}
boxplot(B21)
quantile(B21,c(0.025,0.975))
hist(B21, breaks = 30)
dB21<-c(Bmean=mean(B21),quantile(B21,c(0.025,0.975)),n=nrow(MoReES2))
dB21
#B22: Reptiles in YSF, Tropical Moist Forest
MoReYSF<-subset(data,BGroupBiome=="MoistReptilesYSF")
median(MoReYSF$RRS)
R<-10000
B22<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MoReYSF2<-ddply(MoReYSF,.(Nstudy))
  boot.sample <- sample(MoReYSF2$RRS, replace = TRUE)
  B22[kk] <- mean(boot.sample)
}
boxplot(B22)
quantile(B22,c(0.025,0.975))
hist(B22, breaks = 30)
dB22<-c(Bmean=mean(B22),quantile(B22,c(0.025,0.975)),n=nrow(MoReYSF2))
dB22
#B23: Reptiles in MSF, Moist Forest
MoReMSF<-subset(data,BGroupBiome=="MoistReptilesMSF")
median(MoReMSF$RRS)
R<-10000
B23<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MoReMSF2<-ddply(MoReMSF,.(Nstudy))
  boot.sample <- sample(MoReMSF2$RRS, replace = TRUE)
  B23[kk] <- mean(boot.sample)
}
boxplot(B23)
quantile(B23,c(0.025,0.975))
hist(B23, breaks = 30)
dB23<-c(Bmean=mean(B23),quantile(B23,c(0.025,0.975)),n=nrow(MoReMSF2))
dB23
#B24: Reptiles in OSF, Moist
MoReOSF<-subset(data,BGroupBiome=="MoistReptilesOSF")
nrow(MoReOSF)
B24<-rep(mean(MoReOSF$RRS),10000)
dB24<-c(Bmean=mean(B24),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MoReOSF))
dB24
#B25: Birds in ES, Moist
MoBiES<-subset(data,BGroupBiome=="MoistBirdsES")
median(MoBiES$RRS)
nrow(MoBiES)
R<-10000
B25<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MoBiES2<-ddply(MoBiES,.(Nstudy))
  boot.sample <- sample(MoBiES2$RRS, replace = TRUE)
  B25[kk] <- mean(boot.sample)
}
boxplot(B25)
quantile(B25,c(0.025,0.975))
hist(B25, breaks = 30)
dB25<-c(Bmean=mean(B25),quantile(B25,c(0.025,0.975)),n=nrow(MoBiES2))
dB25
#B26: Birds in YSF, Moist
MoBiYSF<-subset(data,BGroupBiome=="MoistBirdsYSF")
median(MoBiYSF$RRS)
nrow(MoBiYSF)
R<-10000
B26<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MoBiYSF2<-ddply(MoBiYSF,.(Nstudy))
  boot.sample <- sample(MoBiYSF2$RRS, replace = TRUE)
  B26[kk] <- mean(boot.sample)
}
boxplot(B26)
quantile(B26,c(0.025,0.975))
hist(B26, breaks = 30)
dB26<-c(Bmean=mean(B26),quantile(B26,c(0.025,0.975)),n=nrow(MoBiYSF2))
dB26
#B27: Birds in MSF, Moist
MoBiMSF<-subset(data,BGroupBiome=="MoistBirdsMSF")
median(MoBiMSF$RRS)
nrow(MoBiMSF)
R<-10000
B27<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MoBiMSF2<-ddply(MoBiMSF,.(Nstudy))
  boot.sample <- sample(MoBiMSF2$RRS, replace = TRUE)
  B27[kk] <- mean(boot.sample)
}
boxplot(B27)
quantile(B27,c(0.025,0.975))
hist(B27, breaks = 30)
dB27<-c(Bmean=mean(B27),quantile(B27,c(0.025,0.975)),n=nrow(MoBiMSF2))
dB27
#B28: Birds in OSF, Moist
MoBiOSF<-subset(data,BGroupBiome=="MoistBirdsOSF")
nrow(MoBiOSF)
B28<-rep(mean(MoBiOSF$RRS),10000)
dB28<-c(Bmean=mean(B28),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MoBiOSF))
dB28
#B29: Mammals in ES, Moist
MoMaES<-subset(data,BGroupBiome=="MoistMammalsES")
median(MoMaES$RRS)
nrow(MoMaES)
R<-10000
B29<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MoMaES2<-ddply(MoMaES,.(Nstudy))
  boot.sample <- sample(MoMaES2$RRS, replace = TRUE)
  B29[kk] <- mean(boot.sample)
}
boxplot(B29)
quantile(B29,c(0.025,0.975))
hist(B29, breaks = 30)
dB29<-c(Bmean=mean(B29),quantile(B29,c(0.025,0.975)),n=nrow(MoMaES2))
dB29
#B30: Mammals in YSF, Moist
MoMaYSF<-subset(data,BGroupBiome=="MoistMammalsYSF")
median(MoMaYSF$RRS)
nrow(MoMaYSF)
R<-10000
B30<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MoMaYSF2<-ddply(MoMaYSF,.(Nstudy))
  boot.sample <- sample(MoMaYSF2$RRS, replace = TRUE)
  B30[kk] <- mean(boot.sample)
}
boxplot(B30)
quantile(B30,c(0.025,0.975))
hist(B30, breaks = 30)
dB30<-c(Bmean=mean(B30),quantile(B30,c(0.025,0.975)),n=nrow(MoMaYSF2))
dB30
#B31: Mammals in MSF, Moist
MoMaMSF<-subset(data,BGroupBiome=="MoistMammalsMSF")
median(MoMaMSF$RRS)
nrow(MoMaMSF)
R<-10000
B31<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  MoMaMSF2<-ddply(MoMaYSF,.(Nstudy))
  boot.sample <- sample(MoMaMSF2$RRS, replace = TRUE)
  B31[kk] <- mean(boot.sample)
}
boxplot(B31)
quantile(B31,c(0.025,0.975))
hist(B31, breaks = 30)
dB31<-c(Bmean=mean(B31),quantile(B31,c(0.025,0.975)),n=nrow(MoMaMSF2))
dB31
#B32: Mammals in OSF, Moist
MoMaOSF<-subset(data,BGroupBiome=="MoistMammalsOSF")
nrow(MoMaOSF)
B32<-rep(mean(MoMaOSF$RRS),10000)
dB32<-c(Bmean=mean(B32),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MoMaOSF))
dB32

##~~ Tropical Dry Forest
#B33: Amphibians in Early Succession (ES), Tropical Dry Forest
DrAmES<-subset(data,BGroupBiome=="DryAmphibiansES")
nrow(DrAmES)
B33<-rep(mean(DrAmES$RRS),10000)
dB33<-c(Bmean=mean(B33),"2.5%" = NaN,"97.5%" = NaN,n=nrow(DrAmES))
dB33
#B34: Amphibians in Young Secondary Forest (YSF), Dry
DrAmYSF<-subset(data,BGroupBiome=="DryAmphibiansYSF")
nrow(DrAmYSF)
B34<-rep(mean(DrAmYSF$RRS),10000)
dB34<-c(Bmean=mean(B34),"2.5%" = NaN,"97.5%" = NaN,n=nrow(DrAmYSF))
dB34
#B35: Amphibians in Mid-successional Secondary Forest (MSF), Dry
DrAmMSF<-subset(data,BGroupBiome=="DryAmphibiansMSF")
nrow(DrAmMSF)
B35<-rep(mean(DrAmMSF$RRS),10000)
dB35<-c(Bmean=mean(B35),"2.5%" = NaN,"97.5%" = NaN,n=nrow(DrAmMSF))
dB35
#B36: Amphibians in Old Secondary Forest (OSF), Dry (there is no data)
B36<-rep(NaN,10000)
dB36<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB36
#B37: Reptiles in ES, Dry
DrReES<-subset(data,BGroupBiome=="DryReptilesES")
nrow(DrReES)
B37<-rep(mean(DrReES$RRS),10000)
dB37<-c(Bmean=mean(B37),"2.5%" = NaN,"97.5%" = NaN,n=nrow(DrReES))
dB37
#B38: Reptiles in YSF, Dry
DrReYSF<-subset(data,BGroupBiome=="DryReptilesYSF")
nrow(DrReYSF)
B38<-rep(mean(DrReYSF$RRS),10000)
dB38<-c(Bmean=mean(B38),"2.5%" = NaN,"97.5%" = NaN,n=nrow(DrReYSF))
dB38
#B39: Reptiles in MSF, Dry
DrReMSF<-subset(data,BGroupBiome=="DryReptilesMSF")
nrow(DrReMSF)
B39<-rep(mean(DrReMSF$RRS),10000)
dB39<-c(Bmean=mean(B39),"2.5%" = NaN,"97.5%" = NaN,n=nrow(DrReMSF))
dB39
#B40: Reptiles in OSF, Dry
B40<-rep(NaN,10000)
dB40<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB40
#B41: Birds in ES, Dry
DrBiES<-subset(data,BGroupBiome=="DryBirdsES")
nrow(DrBiES)
B41<-rep(mean(DrBiES$RRS),10000)
dB41<-c(Bmean=mean(B41),"2.5%" = NaN,"97.5%" = NaN,n=nrow(DrBiES))
dB41
#B42: Birds in YSF, Dry
DrBiYSF<-subset(data,BGroupBiome=="DryBirdsYSF")
nrow(DrBiYSF)
B42<-rep(mean(DrBiYSF$RRS),10000)
dB42<-c(Bmean=mean(B42),"2.5%" = NaN,"97.5%" = NaN,n=nrow(DrBiYSF))
dB42
#B43: Birds in MSF, Dry
DrBiMSF<-subset(data,BGroupBiome=="DryBirdsMSF")
nrow(DrBiMSF)
B43<-rep(mean(DrBiMSF$RRS),10000)
dB43<-c(Bmean=mean(B43),"2.5%" = NaN,"97.5%" = NaN,n=nrow(DrBiMSF))
dB43
#B44: Birds in OSF, Dry (no data)
B44<-rep(NaN,10000)
dB44<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB44
#B45: Mammals in ES, Dry
DrMaES<-subset(data,BGroupBiome=="DryMammalsES")
nrow(DrMaES)
B45<-rep(mean(DrMaES$RRS),10000)
dB45<-c(Bmean=mean(B45),"2.5%" = NaN,"97.5%" = NaN,n=nrow(DrMaES))
dB45
#B46: Mammals in YSF, Dry
DrMaYSF<-subset(data,BGroupBiome=="DryMammalsYSF")
nrow(DrMaYSF)
B46<-rep(mean(DrMaYSF$RRS),10000)
dB46<-c(Bmean=mean(B46),"2.5%" = NaN,"97.5%" = NaN,n=nrow(DrMaYSF))
dB46
#B47: Mammals in MSF, Dry
DrMaMSF<-subset(data,BGroupBiome=="DryMammalsMSF")
nrow(DrMaMSF)
B47<-rep(mean(DrMaMSF$RRS),10000)
dB47<-c(Bmean=mean(B47),"2.5%" = NaN,"97.5%" = NaN,n=nrow(DrMaMSF))
dB47
#B48: Mammals in OSF, Dry
DrMaOSF<-subset(data,BGroupBiome=="DryMammalsOSF")
nrow(DrMaOSF)
B48<-rep(mean(DrMaOSF$RRS),10000)
dB48<-c(Bmean=mean(B48),"2.5%" = NaN,"97.5%" = NaN,n=nrow(DrMaOSF))
dB48

##~~ Tropical Grasslands and Savannahs (ej. Orinoquia)
#B49: Amphibians in Early Succession (ES), Savannahs
SaAmES<-subset(data,BGroupBiome=="SavannahsAmphibiansES")
nrow(SaAmES)
B49<-rep(mean(SaAmES$RRS),10000)
dB49<-c(Bmean=mean(B49),"2.5%" = NaN,"97.5%" = NaN,n=nrow(SaAmES))
dB49
#B50: Amphibians in Young Secondary Forest (YSF), Savannahs (there is no data)
B50<-rep(NaN,10000)
dB50<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB50
#B51: Amphibians in Mid-successional Secondary Forest (MSF), Savannahs
SaAmMSF<-subset(data,BGroupBiome=="SavannahsAmphibiansMSF")
nrow(SaAmMSF)
B51<-rep(mean(SaAmMSF$RRS),10000)
dB51<-c(Bmean=mean(B51),"2.5%" = NaN,"97.5%" = NaN,n=nrow(SaAmMSF))
dB51
#B52: Amphibians in Old Secondary Forest (OSF), Savannahs (there is no data)
B52<-rep(NaN,10000)
dB52<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB52
#B53: Reptiles in ES, Savannahs
SaReES<-subset(data,BGroupBiome=="SavannahsReptilesES")
nrow(SaReES)
B53<-rep(mean(SaReES$RRS),10000)
dB53<-c(Bmean=mean(B53),"2.5%" = NaN,"97.5%" = NaN,n=nrow(SaReES))
dB53
#B54: Reptiles in YSF, Savannahs
B54<-rep(NaN,10000)
dB54<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB54
#B55: Reptiles in MSF, Savannahs
SaReMSF<-subset(data,BGroupBiome=="SavannahsReptilesMSF")
nrow(SaReMSF)
B55<-rep(mean(SaReMSF$RRS),10000)
dB55<-c(Bmean=mean(B55),"2.5%" = NaN,"97.5%" = NaN,n=nrow(SaReMSF))
dB55
#B56: Reptiles in OSF, NT_Neotropical (there is no data)
B56<-rep(NaN,10000)
dB56<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB56
#B57: Birds in ES, Savannahs
SaBiES<-subset(data,BGroupBiome=="SavannahsBirdsES")
nrow(SaBiES)
B57<-rep(mean(SaBiES$RRS),10000)
dB57<-c(Bmean=mean(B57),"2.5%" = NaN,"97.5%" = NaN,n=nrow(SaBiES))
dB57
#B58: Birds in YSF, Savannahs
SaBiYSF<-subset(data,BGroupBiome=="SavannahsBirdsYSF")
nrow(SaBiYSF)
B58<-rep(mean(SaBiYSF$RRS),10000)
dB58<-c(Bmean=mean(B58),"2.5%" = NaN,"97.5%" = NaN,n=nrow(SaBiYSF))
dB58
#B59: Birds in MSF, Savannahs
SaBiMSF<-subset(data,BGroupBiome=="SavannahsBirdsMSF")
nrow(SaBiMSF)
B59<-rep(mean(SaBiMSF$RRS),10000)
dB59<-c(Bmean=mean(B59),"2.5%" = NaN,"97.5%" = NaN,n=nrow(SaBiMSF))
dB59
#B60: Birds in OSF, Savannahs (no data)
B60<-rep(NaN,10000)
dB60<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB60
#B61: Mammals in ES, Savannahs
SaMaES<-subset(data,BGroupBiome=="SavannahsMammalsES")
nrow(SaMaES)
B61<-rep(mean(SaMaES$RRS),10000)
dB61<-c(Bmean=mean(B61),"2.5%" = NaN,"97.5%" = NaN,n=nrow(SaMaES))
dB61
#B62: Mammals in YSF, Savannahs (no data)
B62<-rep(NaN,10000)
dB62<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB62
#B63: Mammals in MSF, Savannahs
SaMaMSF<-subset(data,BGroupBiome=="SavannahsMammalsMSF")
nrow(SaMaMSF)
B63<-rep(mean(SaMaMSF$RRS),10000)
dB63<-c(Bmean=mean(B63),"2.5%" = NaN,"97.5%" = NaN,n=nrow(SaMaMSF))
dB63
#B64: Mammals in OSF, Savannahs (no data)
B64<-rep(NaN,10000)
dB64<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB64


##~~ Geographical Condition~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#B65: Amphibians in Early Succession (ES), Continent
CoAmES<-subset(data,BGroupContIsla=="ContinentAmphibiansES")
median(CoAmES$RRS)
R<-10000
B65<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  CoAmES2<-ddply(CoAmES,.(Nstudy))
  boot.sample <- sample(CoAmES2$RRS, replace = TRUE)
  B65[kk] <- mean(boot.sample)
}
boxplot(B65)
quantile(B65,c(0.025,0.975))
hist(B65, breaks = 30)
dB65<-c(Bmean=mean(B65),quantile(B65,c(0.025,0.975)),n=nrow(CoAmES2))
dB65
#B66: Amphibians in Young Secondary Forest (YSF), Continent
CoAmYSF<-subset(data,BGroupContIsla=="ContinentAmphibiansYSF")
median(CoAmYSF$RRS)
R<-10000
B66<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  CoAmYSF2<-ddply(CoAmYSF,.(Nstudy))
  boot.sample <- sample(CoAmYSF2$RRS, replace = TRUE)
  B66[kk] <- mean(boot.sample)
}
boxplot(B66)
quantile(B66,c(0.025,0.975))
hist(B66, breaks = 30)
dB66<-c(Bmean=mean(B66),quantile(B66,c(0.025,0.975)),n=nrow(CoAmYSF2))
dB66
#B67: Amphibians in Mid-successional Secondary Forest (MSF), Continent
CoAmMSF<-subset(data,BGroupContIsla=="ContinentAmphibiansMSF")
median(CoAmMSF$RRS)
R<-10000
B67<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  CoAmMSF2<-ddply(CoAmMSF,.(Nstudy))
  boot.sample <- sample(CoAmMSF2$RRS, replace = TRUE)
  B67[kk] <- mean(boot.sample)
}
boxplot(B67)
quantile(B67,c(0.025,0.975))
hist(B67, breaks = 30)
dB67<-c(Bmean=mean(B67),quantile(B67,c(0.025,0.975)),n=nrow(CoAmMSF2))
dB67
#B68: Amphibians in Old Secondary Forest (OSF), Continent
CoAmOSF<-subset(data,BGroupContIsla=="ContinentAmphibiansOSF")
nrow(CoAmOSF)
B68<-rep(mean(CoAmOSF$RRS),10000)
dB68<-c(Bmean=mean(B68),"2.5%" = NaN,"97.5%" = NaN,n=nrow(CoAmOSF))
dB68
#B69: Reptiles in ES, Continent 
CoReES<-subset(data,BGroupContIsla=="ContinentReptilesES")
median(CoReES$RRS)
R<-10000
B69<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  CoReES2<-ddply(CoReES,.(Nstudy))
  boot.sample <- sample(CoReES2$RRS, replace = TRUE)
  B69[kk] <- mean(boot.sample)
}
boxplot(B69)
quantile(B69,c(0.025,0.975))
hist(B69, breaks = 30)
dB69<-c(Bmean=mean(B69),quantile(B69,c(0.025,0.975)),n=nrow(CoReES2))
dB69
#B70: Reptiles in YSF, Continent
CoReYSF<-subset(data,BGroupContIsla=="ContinentReptilesYSF")
median(CoReYSF$RRS)
R<-10000
B70<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  CoReYSF2<-ddply(CoReYSF,.(Nstudy))
  boot.sample <- sample(CoReYSF2$RRS, replace = TRUE)
  B70[kk] <- mean(boot.sample)
}
boxplot(B70)
quantile(B70,c(0.025,0.975))
hist(B70, breaks = 30)
dB70<-c(Bmean=mean(B70),quantile(B70,c(0.025,0.975)),n=nrow(CoReYSF2))
dB70
#B71: Reptiles in MSF, Continental
CoReMSF<-subset(data,BGroupContIsla=="ContinentReptilesMSF")
median(CoReMSF$RRS)
R<-10000
B71<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  CoReMSF2<-ddply(CoReMSF,.(Nstudy))
  boot.sample <- sample(CoReMSF2$RRS, replace = TRUE)
  B71[kk] <- mean(boot.sample)
}
boxplot(B71)
quantile(B71,c(0.025,0.975))
hist(B71, breaks = 30)
dB71<-c(Bmean=mean(B71),quantile(B71,c(0.025,0.975)),n=nrow(CoReMSF2))
dB71
#B72: Reptiles in OSF, Continent (there is no data)
B72<-rep(NaN,10000)
dB72<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB72
#B73: Birds in ES, Continent 
CoBiES<-subset(data,BGroupContIsla=="ContinentBirdsES")
median(CoBiES$RRS)
R<-10000
B73<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  CoBiES2<-ddply(CoBiES,.(Nstudy))
  boot.sample <- sample(CoBiES2$RRS, replace = TRUE)
  B73[kk] <- mean(boot.sample)
}
boxplot(B73)
quantile(B73,c(0.025,0.975))
hist(B73, breaks = 30)
dB73<-c(Bmean=mean(B73),quantile(B73,c(0.025,0.975)),n=nrow(CoBiES2))
dB73
#B74: Birds in YSF, Continent
CoBiYSF<-subset(data,BGroupContIsla=="ContinentBirdsYSF")
median(CoBiYSF$RRS)
nrow(CoBiYSF)
R<-10000
B74<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  CoBiYSF2<-ddply(CoBiYSF,.(Nstudy))
  boot.sample <- sample(CoBiYSF2$RRS, replace = TRUE)
  B74[kk] <- mean(boot.sample)
}
boxplot(B74)
quantile(B74,c(0.025,0.975))
hist(B74, breaks = 30)
dB74<-c(Bmean=mean(B74),quantile(B74,c(0.025,0.975)),n=nrow(CoBiYSF2))
dB74
#B75: Birds in MSF, Continent
CoBiMSF<-subset(data,BGroupContIsla=="ContinentBirdsMSF")
median(CoBiMSF$RRS)
nrow(CoBiMSF)
R<-10000
B75<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  CoBiMSF2<-ddply(CoBiMSF,.(Nstudy))
  boot.sample <- sample(CoBiMSF2$RRS, replace = TRUE)
  B75[kk] <- mean(boot.sample)
}
boxplot(B75)
quantile(B75,c(0.025,0.975))
hist(B75, breaks = 30)
dB75<-c(Bmean=mean(B75),quantile(B75,c(0.025,0.975)),n=nrow(CoBiMSF2))
dB75
#B76: Birds in OSF, Continent
CoBiOSF<-subset(data,BGroupContIsla=="ContinentBirdsOSF")
nrow(CoBiOSF)
B76<-rep(mean(CoBiOSF$RRS),10000)
dB76<-c(Bmean=mean(B76),"2.5%" = NaN,"97.5%" = NaN,n=nrow(CoBiOSF))
dB76
#B77: Mammals in ES,  Continent
CoMaES<-subset(data,BGroupContIsla=="ContinentMammalsES")
median(CoMaES$RRS)
nrow(CoMaES)
R<-10000
B77<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  CoMaES2<-ddply(CoMaES,.(Nstudy))
  boot.sample <- sample(CoMaES2$RRS, replace = TRUE)
  B77[kk] <- mean(boot.sample)
}
boxplot(B77)
quantile(B77,c(0.025,0.975))
hist(B77, breaks = 30)
dB77<-c(Bmean=mean(B77),quantile(B77,c(0.025,0.975)),n=nrow(CoMaES2))
dB77
#B78: Mammals in YSF, Continent
CoMaYSF<-subset(data,BGroupContIsla=="ContinentMammalsYSF")
median(CoMaYSF$RRS)
nrow(CoMaYSF)
R<-10000
B78<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  CoMaYSF2<-ddply(CoMaYSF,.(Nstudy))
  boot.sample <- sample(CoMaYSF2$RRS, replace = TRUE)
  B78[kk] <- mean(boot.sample)
}
boxplot(B78)
quantile(B78,c(0.025,0.975))
hist(B78, breaks = 30)
dB78<-c(Bmean=mean(B78),quantile(B78,c(0.025,0.975)),n=nrow(CoMaYSF2))
dB78
#B79: Mammals in MSF, Continent
CoMaMSF<-subset(data,BGroupContIsla=="ContinentMammalsMSF")
median(CoMaMSF$RRS)
nrow(CoMaMSF)
R<-10000
B79<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  CoMaMSF2<-ddply(CoMaMSF,.(Nstudy))
  boot.sample <- sample(CoMaMSF2$RRS, replace = TRUE)
  B79[kk] <- mean(boot.sample)
}
boxplot(B79)
quantile(B79,c(0.025,0.975))
hist(B79, breaks = 30)
dB79<-c(Bmean=mean(B79),quantile(B79,c(0.025,0.975)),n=nrow(CoMaMSF2))
dB79
#B80: Mammals in OSF, Continent
CoMaOSF<-subset(data,BGroupContIsla=="ContinentMammalsOSF")
nrow(CoMaOSF)
B80<-rep(mean(CoMaOSF$RRS),10000)
dB80<-c(Bmean=mean(B80),"2.5%" = NaN,"97.5%" = NaN,n=nrow(CoMaOSF))
dB80

##~~ Islands
#B81: Amphibians in Early Succession (ES), Islands
IsAmES<-subset(data,BGroupContIsla=="IslandAmphibiansES")
nrow(IsAmES)
B81<-rep(mean(IsAmES$RRS),10000)
dB81<-c(Bmean=mean(B81),"2.5%" = NaN,"97.5%" = NaN,n=nrow(IsAmES))
dB81
#B82: Amphibians in Young Secondary Forest (YSF), Islands
IsAmYSF<-subset(data,BGroupContIsla=="IslandAmphibiansYSF")
nrow(IsAmYSF)
B82<-rep(mean(IsAmYSF$RRS),10000)
dB82<-c(Bmean=mean(B82),"2.5%" = NaN,"97.5%" = NaN,n=nrow(IsAmYSF))
dB82
#B83: Amphibians in Mid-successional Secondary Forest (MSF), Island
IsAmMSF<-subset(data,BGroupContIsla=="IslandAmphibiansMSF")
nrow(IsAmMSF)
B83<-rep(mean(IsAmMSF$RRS),10000)
dB83<-c(Bmean=mean(B83),"2.5%" = NaN,"97.5%" = NaN,n=nrow(IsAmYSF))
dB83
#B84: Amphibians in Old Secondary Forest (OSF), Island
B84<-rep(NaN,10000)
dB84<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB84
#B85: Reptiles in ES, Island
IsReES<-subset(data,BGroupContIsla=="IslandReptilesES")
nrow(IsReES)
R<-10000
B85<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  IsReES2<-ddply(IsReES,.(Nstudy))
  boot.sample <- sample(IsReES2$RRS, replace = TRUE)
  B85[kk] <- mean(boot.sample)
}
boxplot(B85)
quantile(B85,c(0.025,0.975))
hist(B85, breaks = 30)
dB85<-c(Bmean=mean(B85),quantile(B85,c(0.025,0.975)),n=nrow(IsReES2))
dB85
#B86: Reptiles in YSF, Island
IsReYSF<-subset(data,BGroupContIsla=="IslandReptilesYSF")
nrow(IsReYSF)
B86<-rep(mean(IsReYSF$RRS),10000)
dB86<-c(Bmean=mean(B86),"2.5%" = NaN,"97.5%" = NaN,n=nrow(IsReYSF))
dB86
#B87: Reptiles in MSF, Island
IsReMSF<-subset(data,BGroupContIsla=="IslandReptilesMSF")
nrow(IsReMSF)
B87<-rep(mean(IsReMSF$RRS),10000)
dB87<-c(Bmean=mean(B87),"2.5%" = NaN,"97.5%" = NaN,n=nrow(IsReMSF))
dB87
#B88: Reptiles in OSF, Island
IsReOSF<-subset(data,BGroupContIsla=="IslandReptilesOSF")
nrow(IsReOSF)
B88<-rep(mean(IsReOSF$RRS),10000)
dB88<-c(Bmean=mean(B88),"2.5%" = NaN,"97.5%" = NaN,n=nrow(IsReOSF))
dB88
#B89: Birds in ES, Island
IsBiES<-subset(data,BGroupContIsla=="IslandBirdsES")
median(IsBiES$RRS)
nrow(IsBiES)
R<-10000
B89<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  IsBiES2<-ddply(IsBiES,.(Nstudy))
  boot.sample <- sample(IsBiES2$RRS, replace = TRUE)
  B89[kk] <- mean(boot.sample)
}
boxplot(B89)
quantile(B89,c(0.025,0.975))
hist(B89, breaks = 30)
dB89<-c(Bmean=mean(B89),quantile(B89,c(0.025,0.975)),n=nrow(IsBiES2))
dB89
#B90: Birds in YSF, Island
IsBiYSF<-subset(data,BGroupContIsla=="IslandBirdsYSF")
median(IsBiYSF$RRS)
nrow(IsBiYSF)
R<-10000
B90<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  IsBiYSF2<-ddply(IsBiYSF,.(Nstudy))
  boot.sample <- sample(IsBiYSF2$RRS, replace = TRUE)
  B90[kk] <- mean(boot.sample)
}
boxplot(B90)
quantile(B90,c(0.025,0.975))
hist(B90, breaks = 30)
dB90<-c(Bmean=mean(B90),quantile(B90,c(0.025,0.975)),n=nrow(IsBiYSF))
dB90
#B91: Birds in MSF, Island)
IsBiMSF<-subset(data,BGroupContIsla=="IslandBirdsMSF")
median(IsBiMSF$RRS)
nrow(IsBiMSF)
R<-10000
B91<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  IsBiMSF2<-ddply(IsBiMSF,.(Nstudy))
  boot.sample <- sample(IsBiMSF2$RRS, replace = TRUE)
  B91[kk] <- mean(boot.sample)
}
boxplot(B91)
quantile(B91,c(0.025,0.975))
hist(B91, breaks = 30)
dB91<-c(Bmean=mean(B91),quantile(B91,c(0.025,0.975)),n=nrow(IsBiMSF2))
dB91
#B92: Birds in OSF, Islands
B92<-rep(NaN,10000)
dB92<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB92
#B93: Mammals in ES,  Islands
IsMaES<-subset(data,BGroupContIsla=="IslandMammalsES")
nrow(IsMaES)
B93<-rep(mean(IsMaES$RRS),10000)
dB93<-c(Bmean=mean(B93),"2.5%" = NaN,"97.5%" = NaN,n=nrow(IsMaES))
dB93
#B94: Mammals in YSF, Island
IsMaYSF<-subset(data,BGroupContIsla=="IslandMammalsYSF")
nrow(IsMaYSF)
B94<-rep(mean(IsMaYSF$RRS),10000)
dB94<-c(Bmean=mean(B94),"2.5%" = NaN,"97.5%" = NaN,n=nrow(IsMaYSF))
dB94
#B95: Mammals in MSF, Island (no data)
B95<-rep(NaN,10000)
dB95<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB95
#B96: Mammals in OSF, Island (no data)
IsMaOSF<-subset(data,BGroupContIsla=="IslandMammalsOSF")
nrow(IsMaOSF)
B96<-rep(mean(IsMaOSF$RRS),10000)
dB96<-c(Bmean=mean(B96),"2.5%" = NaN,"97.5%" = NaN,n=nrow(IsMaOSF))
dB96

#~~~~Zoogeographical realm ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary(data$BGroupBRealm)
##~~ Panamanian Zoogeographic Realm
#B97: Amphibians in Early Succession (ES), Panamanian
PNAmES<-subset(data,BGroupBRealm=="PN_PanamanianAmphibiansES")
median(PNAmES$RRS)
R<-10000
B97<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  PNAmES2<-ddply(PNAmES,.(Nstudy))
  boot.sample <- sample(PNAmES2$RRS, replace = TRUE)
  B97[kk] <- mean(boot.sample)
}
boxplot(B97)
quantile(B97,c(0.025,0.975))
hist(B97, breaks = 30)
dB97<-c(Bmean=mean(B97),quantile(B97,c(0.025,0.975)),n=nrow(PNAmES2))
dB97
#B98: Amphibians in Young Secondary Forest (YSF), Panamanian
PNAmYSF<-subset(data,BGroupBRealm=="PN_PanamanianAmphibiansYSF")
median(PNAmYSF$RRS)
nrow(PNAmYSF)
R<-10000
B98<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  PNAmYSF2<-ddply(PNAmYSF,.(Nstudy))
  boot.sample <- sample(PNAmYSF2$RRS, replace = TRUE)
  B98[kk] <- mean(boot.sample)
}
boxplot(B98)
quantile(B98,c(0.025,0.975))
hist(B98, breaks = 30)
dB98<-c(Bmean=mean(B98),quantile(B98,c(0.025,0.975)),n=nrow(PNAmYSF2))
dB98
#B99: Amphibians in Mid-successional Secondary Forest (MSF), Panamanian
PNAmMSF<-subset(data,BGroupBRealm=="PN_PanamanianAmphibiansMSF")
median(PNAmMSF$RRS)
R<-10000
B99<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  PNAmMSF2<-ddply(PNAmMSF,.(Nstudy))
  boot.sample <- sample(PNAmMSF2$RRS, replace = TRUE)
  B99[kk] <- mean(boot.sample)
}
boxplot(B99)
quantile(B99,c(0.025,0.975))
hist(B99, breaks = 30)
dB99<-c(Bmean=mean(B99),quantile(B99,c(0.025,0.975)),n=nrow(PNAmMSF2))
dB99
#B100: Amphibians in Old Secondary Forest (OSF), PN_Panamanian
PNAmOSF<-subset(data,BGroupBRealm=="PN_PanamanianAmphibiansOSF")
nrow(PNAmOSF)
B100<-rep(mean(PNAmOSF$RRS),10000)
dB100<-c(Bmean=mean(B100),"2.5%" = NaN,"97.5%" = NaN,n=nrow(PNAmOSF))
dB100
#B101: Reptiles in ES, PN_Panamanian
PNReES<-subset(data,BGroupBRealm=="PN_PanamanianReptilesES")
median(PNReES$RRS)
R<-10000
B101<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  PNReES2<-ddply(PNReES,.(Nstudy))
  boot.sample <- sample(PNReES2$RRS, replace = TRUE)
  B101[kk] <- mean(boot.sample)
}
boxplot(B101)
quantile(B101,c(0.025,0.975))
hist(B101, breaks = 30)
dB101<-c(Bmean=mean(B101),quantile(B101,c(0.025,0.975)),n=nrow(PNReES2))
dB101
#B102: Reptiles in YSF, PN_Panamanian
PNReYSF<-subset(data,BGroupBRealm=="PN_PanamanianReptilesYSF")
median(PNReYSF$RRS)
R<-10000
B102<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  PNReYSF2<-ddply(PNReYSF,.(Nstudy))
  boot.sample <- sample(PNReYSF2$RRS, replace = TRUE)
  B102[kk] <- mean(boot.sample)
}
boxplot(B102)
quantile(B102,c(0.025,0.975))
hist(B102, breaks = 30)
dB102<-c(Bmean=mean(B102),quantile(B102,c(0.025,0.975)),n=nrow(PNReYSF2))
dB102
#B103: Reptiles in MSF, PN_Panamanian
PNReMSF<-subset(data,BGroupBRealm=="PN_PanamanianReptilesMSF")
median(PNReMSF$RRS)
R<-10000
B103<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  PNReMSF2<-ddply(PNReMSF,.(Nstudy))
  boot.sample <- sample(PNReMSF2$RRS, replace = TRUE)
  B103[kk] <- mean(boot.sample)
}
boxplot(B103)
quantile(B103,c(0.025,0.975))
hist(B103, breaks = 30)
dB103<-c(Bmean=mean(B103),quantile(B103,c(0.025,0.975)),n=nrow(PNReMSF2))
dB103
#B104: Reptiles in OSF, PN_Panamanian
PNReOSF<-subset(data,BGroupBRealm=="PN_PanamanianReptilesOSF")
nrow(PNReOSF)
B104<-rep(mean(PNReOSF$RRS),10000)
dB104<-c(Bmean=mean(B104),"2.5%" = NaN,"97.5%" = NaN,n=nrow(PNReOSF))
dB104
#B105: Birds in ES, PN_Panamanian
PNBiES<-subset(data,BGroupBRealm=="PN_PanamanianBirdsES")
nrow(PNBiES)
B105<-rep(mean(PNBiES$RRS),10000)
dB105<-c(Bmean=mean(B105),"2.5%" = NaN,"97.5%" = NaN,n=nrow(PNBiES))
dB105
#B106: Birds in YSF, PN_Panamanian
PNBiYSF<-subset(data,BGroupBRealm=="PN_PanamanianBirdsYSF")
median(PNBiYSF$RRS)
nrow(PNBiYSF)
R<-10000
B106<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  PNBiYSF2<-ddply(PNBiYSF,.(Nstudy))
  boot.sample <- sample(PNBiYSF2$RRS, replace = TRUE)
  B106[kk] <- mean(boot.sample)
}
boxplot(B106)
quantile(B106,c(0.025,0.975))
hist(B106, breaks = 30)
dB106<-c(Bmean=mean(B106),quantile(B106,c(0.025,0.975)),n=nrow(PNBiYSF2))
dB106
#B107: Birds in MSF, PN_Panamanian
PNBiMSF<-subset(data,BGroupBRealm=="PN_PanamanianBirdsMSF")
nrow(PNBiMSF)
B107<-rep(mean(PNBiMSF$RRS),10000)
dB107<-c(Bmean=mean(B107),"2.5%" = NaN,"97.5%" = NaN,n=nrow(PNBiMSF))
dB107
#B108: Birds in OSF, PN_Panamanian
PNBiOSF<-subset(data,BGroupBRealm=="PN_PanamanianBirdsOSF")
nrow(PNBiOSF)
B108<-rep(mean(PNBiOSF$RRS),10000)
dB108<-c(Bmean=mean(B108),"2.5%" = NaN,"97.5%" = NaN,n=nrow(PNBiOSF))
dB108
#B109: Mammals in ES, PN_Panamanian
PNMaES<-subset(data,BGroupBRealm=="PN_PanamanianMammalsES")
nrow(PNMaES)
B109<-rep(mean(PNMaES$RRS),10000)
dB109<-c(Bmean=mean(B109),"2.5%" = NaN,"97.5%" = NaN,n=nrow(PNMaES))
dB109
#B110: Mammals in YSF, PN_Panamanian
PNMaYSF<-subset(data,BGroupBRealm=="PN_PanamanianMammalsYSF")
median(PNMaYSF$RRS)
nrow(PNMaYSF)
R<-10000
B110<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  PNMaYSF2<-ddply(PNMaYSF,.(Nstudy))
  boot.sample <- sample(PNMaYSF2$RRS, replace = TRUE)
  B110[kk] <- mean(boot.sample)
}
boxplot(B110)
quantile(B110,c(0.025,0.975))
hist(B110, breaks = 30)
dB110<-c(Bmean=mean(B110),quantile(B110,c(0.025,0.975)),n=nrow(PNMaYSF2))
dB110
#B111: Mammals in MSF, PN_Panamanian
PNMaMSF<-subset(data,BGroupBRealm=="PN_PanamanianMammalsMSF")
nrow(PNMaMSF)
B111<-rep(mean(PNMaMSF$RRS),10000)
dB111<-c(Bmean=mean(B111),"2.5%" = NaN,"97.5%" = NaN,n=nrow(PNMaMSF))
dB111
#B112: Mammals in OSF, PN_Panamanian
PNMaOSF<-subset(data,BGroupBRealm=="PN_PanamanianMammalsOSF")
nrow(PNMaOSF)
B112<-rep(mean(PNMaOSF$RRS),10000)
dB112<-c(Bmean=mean(B112),"2.5%" = NaN,"97.5%" = NaN,n=nrow(PNMaOSF))
dB112

##~~ NT_Neotropical Biogeographic Realm
#B113: Amphibians in Early Succession (ES), NT_Neotropical
NTAmES<-subset(data,BGroupBRealm=="NT_NeotropicalAmphibiansES")
nrow(NTAmES)
R<-10000
B113<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  NTAmES2<-ddply(NTAmES,.(Nstudy))
  boot.sample <- sample(NTAmES2$RRS, replace = TRUE)
  B113[kk] <- mean(boot.sample)
}
boxplot(B113)
quantile(B113,c(0.025,0.975))
hist(B113, breaks = 30)
dB113<-c(Bmean=mean(B113),quantile(B113,c(0.025,0.975)),n=nrow(NTAmES2))
dB113
#B114: Amphibians in Young Secondary Forest (YSF), NT_Neotropical
NTAmYSF<-subset(data,BGroupBRealm=="NT_NeotropicalAmphibiansYSF")
nrow(NTAmYSF)
R<-10000
B114<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  NTAmYSF2<-ddply(NTAmYSF,.(Nstudy))
  boot.sample <- sample(NTAmYSF2$RRS, replace = TRUE)
  B114[kk] <- mean(boot.sample)
}
boxplot(B114)
quantile(B114,c(0.025,0.975))
hist(B114, breaks = 30)
dB114<-c(Bmean=mean(B114),quantile(B114,c(0.025,0.975)),n=nrow(NTAmYSF2))
dB114
#B115: Amphibians in Mid-successional Secondary Forest (MSF), NT_Neotropical
NTAmMSF<-subset(data,BGroupBRealm=="NT_NeotropicalAmphibiansMSF")
nrow(NTAmMSF)
B115<-rep(mean(NTAmMSF$RRS),10000)
dB115<-c(Bmean=mean(B115),"2.5%" = NaN,"97.5%" = NaN,n=nrow(NTAmMSF))
dB115
#B116: Amphibians in Old Secondary Forest (OSF), NT_Neotropical (there is no data)
B116<-rep(NaN,10000)
dB116<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB116
#B117: Reptiles in ES, NT_Neotropical
NTReES<-subset(data,BGroupBRealm=="NT_NeotropicalReptilesES")
nrow(NTReES)
B117<-rep(mean(NTReES$RRS),10000)
dB117<-c(Bmean=mean(B117),"2.5%" = NaN,"97.5%" = NaN,n=nrow(NTReES))
dB117
#B118: Reptiles in YSF, NT_Neotropical
NTReYSF<-subset(data,BGroupBRealm=="NT_NeotropicalReptilesYSF")
nrow(NTReYSF)
B118<-rep(mean(NTReYSF$RRS),10000)
dB118<-c(Bmean=mean(B118),"2.5%" = NaN,"97.5%" = NaN,n=nrow(NTReYSF))
dB118
#B119: Reptiles in MSF, NT_Neotropical (there is no data)
B119<-rep(NaN,10000)
dB119<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB119
#B120: Reptiles in OSF, NT_Neotropical (there is no data)
B120<-rep(NaN,10000)
dB120<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB120
#B121: Birds in ES, NT_Neotropical
NTBiES<-subset(data,BGroupBRealm=="NT_NeotropicalBirdsES")
nrow(NTBiES)
B121<-rep(mean(NTBiES$RRS),10000)
dB121<-c(Bmean=mean(B121),"2.5%" = NaN,"97.5%" = NaN,n=nrow(NTBiES))
dB121
#B122: Birds in YSF, NT_Neotropical
NTBiYSF<-subset(data,BGroupBRealm=="NT_NeotropicalBirdsYSF")
median(NTBiYSF$RRS)
nrow(NTBiYSF)
R<-10000
B122<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  NTBiYSF2<-ddply(NTBiYSF,.(Nstudy))
  boot.sample <- sample(NTBiYSF2$RRS, replace = TRUE)
  B122[kk] <- mean(boot.sample)
}
boxplot(B122)
quantile(B122,c(0.025,0.975))
hist(B122, breaks = 30)
dB122<-c(Bmean=mean(B122),quantile(B122,c(0.025,0.975)),n=nrow(NTBiYSF2))
dB122
#B123: Birds in MSF, NT_Neotropical
NTBiMSF<-subset(data,BGroupBRealm=="NT_NeotropicalBirdsMSF")
nrow(NTBiMSF)
B123<-rep(mean(NTBiMSF$RRS),10000)
dB123<-c(Bmean=mean(B123),"2.5%" = NaN,"97.5%" = NaN,n=nrow(NTBiMSF))
dB123
#B124: Birds in OSF, NT_Neotropical (no data)
NTBiOSF<-subset(data,BGroupBRealm=="NT_NeotropicalBirdsOSF")
nrow(NTBiOSF)
B124<-rep(mean(NTBiOSF$RRS),10000)
dB124<-c(Bmean=mean(B124),"2.5%" = NaN,"97.5%" = NaN,n=nrow(NTBiOSF))
dB124
#B125: Mammals in ES, NT_Neotropical
NTMaES<-subset(data,BGroupBRealm=="NT_NeotropicalMammalsES")
nrow(NTMaES)
B125<-rep(mean(NTMaES$RRS),10000)
dB125<-c(Bmean=mean(B125),"2.5%" = NaN,"97.5%" = NaN,n=nrow(NTMaES))
dB125
#B126: Mammals in YSF, NT_Neotropical
NTMaYSF<-subset(data,BGroupBRealm=="NT_NeotropicalMammalsYSF")
median(NTMaYSF$RRS)
nrow(NTMaYSF)
R<-10000
B126<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  NTMaYSF2<-ddply(NTMaYSF,.(Nstudy))
  boot.sample <- sample(NTMaYSF2$RRS, replace = TRUE)
  B126[kk] <- mean(boot.sample)
}
boxplot(B126)
quantile(B126,c(0.025,0.975))
hist(B126, breaks = 30)
dB126<-c(Bmean=mean(B126),quantile(B126,c(0.025,0.975)),n=nrow(NTMaYSF2))
dB126
#B127: Mammals in MSF, NT_Neotropical
NTMaMSF<-subset(data,BGroupBRealm=="NT_NeotropicalMammalsMSF")
nrow(NTMaMSF)
B127<-rep(mean(NTMaMSF$RRS),10000)
dB127<-c(Bmean=mean(B127),"2.5%" = NaN,"97.5%" = NaN,n=nrow(NTMaMSF))
dB127
#B128: Mammals in OSF, NT_Neotropical (no data)
B128<-rep(NaN,10000)
dB128<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB128

##~~ OR_Oriental Biogeographic Realm
#B129: Amphibians in Early Succession (ES), OR_Oriental
ORAmES<-subset(data,BGroupBRealm=="OR_OrientalAmphibiansES")
nrow(ORAmES)
B129<-rep(mean(ORAmES$RRS),10000)
dB129<-c(Bmean=mean(B129),"2.5%" = NaN,"97.5%" = NaN,n=nrow(ORAmES))
dB129
#B130: Amphibians in Young Secondary Forest (YSF), OR_Oriental
ORAmYSF<-subset(data,BGroupBRealm=="OR_OrientalAmphibiansYSF")
nrow(ORAmYSF)
B130<-rep(mean(ORAmYSF$RRS),10000)
dB130<-c(Bmean=mean(B130),"2.5%" = NaN,"97.5%" = NaN,n=nrow(ORAmYSF))
dB130
#B131: Amphibians in Mid-successional Secondary Forest (MSF), OR_Oriental
ORAmMSF<-subset(data,BGroupBRealm=="OR_OrientalAmphibiansMSF")
nrow(ORAmMSF)
B131<-rep(mean(ORAmMSF$RRS),10000)
dB131<-c(Bmean=mean(B131),"2.5%" = NaN,"97.5%" = NaN,n=nrow(ORAmMSF))
dB131
#B1322: Amphibians in Old Secondary Forest (OSF), OR_Oriental (no data)
B132<-rep(NaN,10000)
dB132<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB132
#B133: Reptiles in ES, OR_Oriental (no data)
B133<-rep(NaN,10000)
dB133<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB133
#B134: Reptiles in YSF, OR_Oriental
ORReYSF<-subset(data,BGroupBRealm=="OR_OrientalReptilesYSF")
nrow(ORReYSF)
B134<-rep(mean(ORReYSF$RRS),10000)
dB134<-c(Bmean=mean(B134),"2.5%" = NaN,"97.5%" = NaN,n=nrow(ORReYSF))
dB134
#B135: Reptiles in MSF, OR_Oriental
ORReMSF<-subset(data,BGroupBRealm=="OR_OrientalReptilesMSF")
nrow(ORReMSF)
B135<-rep(mean(ORReMSF$RRS),10000)
dB135<-c(Bmean=mean(B135),"2.5%" = NaN,"97.5%" = NaN,n=nrow(ORReMSF))
dB135
#B136: Reptiles in OSF, OR_Oriental (there is no data)
B136<-rep(NaN,10000)
dB136<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB136
#B137: Birds in ES, OR_Oriental
ORBiES<-subset(data,BGroupBRealm=="OR_OrientalBirdsES")
nrow(ORBiES)
B137<-rep(mean(ORBiES$RRS),10000)
dB137<-c(Bmean=mean(B137),"2.5%" = NaN,"97.5%" = NaN,n=nrow(ORBiES))
dB137
#B138: Birds in YSF, OR_Oriental
ORBiYSF<-subset(data,BGroupBRealm=="OR_OrientalBirdsYSF")
median(ORBiYSF$RRS)
nrow(ORBiYSF)
R<-10000
B138<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  ORBiYSF2<-ddply(ORBiYSF,.(Nstudy))
  boot.sample <- sample(ORBiYSF2$RRS, replace = TRUE)
  B138[kk] <- mean(boot.sample)
}
boxplot(B138)
quantile(B138,c(0.025,0.975))
hist(B138, breaks = 30)
dB138<-c(Bmean=mean(B138),quantile(B138,c(0.025,0.975)),n=nrow(ORBiYSF2))
dB138
#B139: Birds in MSF, OR_Oriental
ORBiMSF<-subset(data,BGroupBRealm=="OR_OrientalBirdsMSF")
median(ORBiMSF$RRS)
nrow(ORBiMSF)
R<-10000
B139<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
##Sample 1 data per study
  ORBiMSF2<-ddply(ORBiMSF,.(Nstudy))
  boot.sample <- sample(ORBiMSF2$RRS, replace = TRUE)
  B139[kk] <- mean(boot.sample)
}
boxplot(B139)
quantile(B139,c(0.025,0.975))
hist(B139, breaks = 30)
dB139<-c(Bmean=mean(B139),quantile(B139,c(0.025,0.975)),n=nrow(ORBiMSF2))
dB139
#B140: Birds in OSF, OR_Oriental (no data)
ORBiOSF<-subset(data,BGroupBRealm=="OR_OrientalBirdsOSF")
nrow(ORBiOSF)
B140<-rep(mean(ORBiOSF$RRS),10000)
dB140<-c(Bmean=mean(B140),"2.5%" = NaN,"97.5%" = NaN,n=nrow(ORBiOSF))
dB140
#B141: Mammals in ES,  OR_Oriental
ORMaES<-subset(data,BGroupBRealm=="OR_OrientalMammalsES")
nrow(ORMaES)
B141<-rep(mean(ORMaES$RRS),10000)
dB141<-c(Bmean=mean(B141),"2.5%" = NaN,"97.5%" = NaN,n=nrow(ORMaES))
dB141
#B142: Mammals in YSF, OR_Oriental
ORMaYSF<-subset(data,BGroupBRealm=="OR_OrientalMammalsYSF")
median(ORMaYSF$RRS)
nrow(ORMaYSF)
B142<-rep(mean(ORMaYSF$RRS),10000)
dB142<-c(Bmean=mean(B142),"2.5%" = NaN,"97.5%" = NaN,n=nrow(ORMaYSF))
dB142
#B143: Mammals in MSF, OR_Oriental
ORMaMSF<-subset(data,BGroupBRealm=="OR_OrientalMammalsMSF")
nrow(ORMaMSF)
B143<-rep(mean(ORMaMSF$RRS),10000)
dB143<-c(Bmean=mean(B143),"2.5%" = NaN,"97.5%" = NaN,n=nrow(ORMaMSF))
dB143
#B144: Mammals in OSF, OR_Oriental
ORMaOSF<-subset(data,BGroupBRealm=="OR_OrientalMammalsOSF")
nrow(ORMaOSF)
B144<-rep(mean(ORMaOSF$RRS),10000)
dB144<-c(Bmean=mean(B144),"2.5%" = NaN,"97.5%" = NaN,n=nrow(ORMaOSF))
dB144

##~~ AF_Afrotropical Biogeographic Realm
#B145: Amphibians in Early Succession (ES), AF_Afrotropical (no data)
B145<-rep(NaN,10000)
dB145<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB145
#B146: Amphibians in Young Secondary Forest (YSF), AF_Afrotropical
B146<-rep(NaN,10000)
dB146<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB146
#B147: Amphibians in Mid-successional Secondary Forest (MSF), AF_Afrotropical
AfAmMSF<-subset(data,BGroupBRealm=="AF_AfrotropicalAmphibiansMSF")
nrow(AfAmMSF)
B147<-rep(mean(AfAmMSF$RRS),10000)
dB147<-c(Bmean=mean(B147),"2.5%" = NaN,"97.5%" = NaN,n=nrow(AfAmMSF))
dB147
#B148: Amphibians in Old Secondary Forest (OSF), AF_Afrotropical (no data)
B148<-rep(NaN,10000)
dB148<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB148
#B149: Reptiles in ES, AF_Afrotropical (no data)
B149<-rep(NaN,10000)
dB149<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB149
#B150: Reptiles in YSF, AF_Afrotropical (no data)
B150<-rep(NaN,10000)
dB150<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB150
#B151: Reptiles in MSF, AF_Afrotropical (no data)
B151<-rep(NaN,10000)
dB151<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB151
#B152: Reptiles in OSF, AF_Afrotropical (there is no data)
B152<-rep(NaN,10000)
dB152<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB152
#B153: Birds in ES, AF_Afrotropical
AfBiES<-subset(data,BGroupBRealm=="AF_AfrotropicalBirdsES")
nrow(AfBiES)
B153<-rep(mean(AfBiES$RRS),10000)
dB153<-c(Bmean=mean(B153),"2.5%" = NaN,"97.5%" = NaN,n=nrow(AfBiES))
dB153
#B154: Birds in YSF, AF_Afrotropical
AfBiYSF<-subset(data,BGroupBRealm=="AF_AfrotropicalBirdsYSF")
nrow(AfBiYSF)
B154<-rep(mean(AfBiYSF$RRS),10000)
dB154<-c(Bmean=mean(B154),"2.5%" = NaN,"97.5%" = NaN,n=nrow(AfBiYSF))
dB154
#B155: Birds in MSF, AF_Afrotropical
AfBiMSF<-subset(data,BGroupBRealm=="AF_AfrotropicalBirdsMSF")
nrow(AfBiMSF)
B155<-rep(mean(AfBiMSF$RRS),10000)
dB155<-c(Bmean=mean(B155),"2.5%" = NaN,"97.5%" = NaN,n=nrow(AfBiMSF))
dB155
#B156: Birds in OSF, AF_Afrotropical
AfBiOSF<-subset(data,BGroupBRealm=="AF_AfrotropicalBirdsOSF")
nrow(AfBiOSF)
B156<-rep(mean(AfBiOSF$RRS),10000)
dB156<-c(Bmean=mean(B156),"2.5%" = NaN,"97.5%" = NaN,n=nrow(AfBiOSF))
dB156
#B157: Mammals in ES,  AF_Afrotropical
AfMaES<-subset(data,BGroupBRealm=="AF_AfrotropicalMammalsES")
nrow(AfMaES)
B157<-rep(mean(AfMaES$RRS),10000)
dB157<-c(Bmean=mean(B157),"2.5%" = NaN,"97.5%" = NaN,n=nrow(AfMaES))
dB157
#B158: Mammals in YSF, AF_Afrotropical
AfMaYSF<-subset(data,BGroupBRealm=="AF_AfrotropicalMammalsYSF")
nrow(AfMaYSF)
B158<-rep(mean(AfMaYSF$RRS),10000)
dB158<-c(Bmean=mean(B158),"2.5%" = NaN,"97.5%" = NaN,n=nrow(AfMaYSF))
dB158
#B159: Mammals in MSF, AF_Afrotropical
AfMaMSF<-subset(data,BGroupBRealm=="AF_AfrotropicalMammalsMSF")
nrow(AfMaMSF)
B159<-rep(mean(AfMaMSF$RRS),10000)
dB159<-c(Bmean=mean(B159),"2.5%" = NaN,"97.5%" = NaN,n=nrow(AfMaMSF))
dB159
#B160: Mammals in OSF, AF_Afrotropical (no data)
B160<-rep(NaN,10000)
dB160<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB160

##~~ NA_Neartic Biogeographic Realm
#B161: Amphibians in Early Succession (ES), NA_Neartic
NAAmES<-subset(data,BGroupBRealm=="NA_NearticAmphibiansES")
nrow(NAAmES)
B161<-rep(mean(NAAmES$RRS),10000)
dB161<-c(Bmean=mean(B161),"2.5%" = NaN,"97.5%" = NaN,n=nrow(NAAmES))
dB161
#B162: Amphibians in Young Secondary Forest (YSF), NA_Neartic
NAAmYSF<-subset(data,BGroupBRealm=="NA_NearticAmphibiansYSF")
nrow(NAAmYSF)
B162<-rep(NAAmYSF$RRS,10000)
dB162<-c(Bmean=mean(B162),"2.5%" = NaN,"97.5%" = NaN,n=nrow(NAAmYSF))
dB162
#B163: Amphibians in Mid-successional Secondary Forest (MSF), NA_Neartic
NAAmMSF<-subset(data,BGroupBRealm=="NA_NearticAmphibiansMSF")
nrow(NAAmMSF)
B163<-rep(NAAmMSF$RRS,10000)
dB163<-c(Bmean=mean(B163),"2.5%" = NaN,"97.5%" = NaN,n=nrow(NAAmMSF))
dB163
#B164: Amphibians in Old Secondary Forest (OSF), NA_Neartic (no data)
B164<-rep(NaN,10000)
dB164<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB164
#B165: Reptiles in ES, NA_Neartic
NAReES<-subset(data,BGroupBRealm=="NA_NearticReptilesES")
nrow(NAReES)
B165<-rep(NAReES$RRS,10000)
dB165<-c(Bmean=mean(B165),"2.5%" = NaN,"97.5%" = NaN,n=nrow(NAReES))
dB165
#B166: Reptiles in YSF, NA_Neartic
NAReYSF<-subset(data,BGroupBRealm=="NA_NearticReptilesYSF")
nrow(NAReYSF)
B166<-rep(NAReYSF$RRS,10000)
dB166<-c(Bmean=mean(B166),"2.5%" = NaN,"97.5%" = NaN,n=nrow(NAReYSF))
dB166
#B167: Reptiles in MSF, NA_Neartic
NAReMSF<-subset(data,BGroupBRealm=="NA_NearticReptilesMSF")
nrow(NAReMSF)
B167<-rep(NAReMSF$RRS,10000)
dB167<-c(Bmean=mean(B167),"2.5%" = NaN,"97.5%" = NaN,n=nrow(NAReMSF))
dB167
#B168: Reptiles in OSF, NA_Neartic (there is no data)
B168<-rep(NaN,10000)
dB168<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB168
#B169: Birds in ES, NA_Neartic (no data)
B169<-rep(NaN,10000)
dB169<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB169
#B170: Birds in YSF, NA_Neartic
NABiYSF<-subset(data,BGroupBRealm=="NA_NearticBirdsYSF")
median(NABiYSF$RRS)
nrow(NABiYSF)
B170<-rep(NABiYSF$RRS,10000)
dB170<-c(Bmean=mean(B170),"2.5%" = NaN,"97.5%" = NaN,n=nrow(NABiYSF))
dB170
#B1171: Birds in MSF, NA_Neartic
NABiMSF<-subset(data,BGroupBRealm=="NA_NearticBirdsMSF")
nrow(NABiMSF)
B171<-rep(NABiMSF$RRS,10000)
dB171<-c(Bmean=mean(B171),"2.5%" = NaN,"97.5%" = NaN,n=nrow(NABiMSF))
dB171
#B172: Birds in OSF, NA_Neartic (no data)
B172<-rep(NaN,10000)
dB172<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB172
#B173: Mammals in ES,  NA_Neartic
NAMaES<-subset(data,BGroupBRealm=="NA_NearticMammalsES")
nrow(NAMaES)
B173<-rep(NAMaES$RRS,10000)
dB173<-c(Bmean=mean(B173),"2.5%" = NaN,"97.5%" = NaN,n=nrow(NAMaES))
dB173
#B174: Mammals in YSF, NA_Neartic
NAMaYSF<-subset(data,BGroupBRealm=="NA_NearticMammalsYSF")
nrow(NAMaYSF)
B174<-rep(NAMaYSF$RRS,10000)
dB174<-c(Bmean=mean(B174),"2.5%" = NaN,"97.5%" = NaN,n=nrow(NAMaYSF))
dB174
#B175: Mammals in MSF, NA_Neartic
NAMaMSF<-subset(data,BGroupBRealm=="NA_NearticMammalsMSF")
nrow(NAMaMSF)
B175<-rep(NAMaMSF$RRS,10000)
dB175<-c(Bmean=mean(B175),"2.5%" = NaN,"97.5%" = NaN,n=nrow(NAMaMSF))
dB175
#B176: Mammals in OSF, NA_Neartic (no data)
B176<-rep(NaN,10000)
dB176<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB176

##~~ OC_Oceanian Biogeographic Realm
#B177: Amphibians in Early Succession (ES), OC_Oceanian
B177<-rep(NaN,10000)
dB177<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB177
#B178: Amphibians in Young Secondary Forest (YSF), OC_Oceanian
B178<-rep(NaN,10000)
dB178<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB178
#B179: Amphibians in Mid-successioMAl Secondary Forest (MSF), OC_Oceanian
B179<-rep(NaN,10000)
dB179<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB179
#B180: Amphibians in Old Secondary Forest (OSF), OC_Oceanian (no data)
B180<-rep(NaN,10000)
dB180<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB180

#B181: Reptiles in ES, OC_Oceanian (no data)
OCReES<-subset(data,BGroupBRealm=="OC_OceanianReptilesES")
nrow(OCReES)
B181<-rep(OCReES$RRS,10000)
dB181<-c(Bmean=mean(B181),"2.5%" = NaN,"97.5%" = NaN,n=nrow(OCReES))
dB181
#B182: Reptiles in YSF, OC_Oceanian
OCReYSF<-subset(data,BGroupBRealm=="OC_OceanianReptilesYSF")
nrow(OCReYSF)
B182<-rep(OCReYSF$RRS,10000)
dB182<-c(Bmean=mean(B182),"2.5%" = NaN,"97.5%" = NaN,n=nrow(OCReYSF))
dB182
#B183: Reptiles in MSF, OC_Oceanian
OCReMSF<-subset(data,BGroupBRealm=="OC_OceanianReptilesMSF")
nrow(OCReMSF)
B183<-rep(OCReMSF$RRS,10000)
dB183<-c(Bmean=mean(B183),"2.5%" = NaN,"97.5%" = NaN,n=nrow(OCReMSF))
dB183
#B184: Reptiles in OSF, OC_Oceanian (there is no data)
B184<-rep(NaN,10000)
dB184<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB184
#B185: Birds in ES, OC_Oceanian
OCBiES<-subset(data,BGroupBRealm=="OC_OceanianBirdsES")
nrow(OCBiES)
B185<-rep(OCBiES$RRS,10000)
dB185<-c(Bmean=mean(B185),"2.5%" = NaN,"97.5%" = NaN,n=nrow(OCBiES))
dB185
#B186: Birds in YSF, OC_Oceanian
OCBiYSF<-subset(data,BGroupBRealm=="OC_OceanianBirdsYSF")
nrow(OCBiYSF)
B186<-rep(OCBiYSF$RRS,10000)
dB186<-c(Bmean=mean(B186),"2.5%" = NaN,"97.5%" = NaN,n=nrow(OCBiYSF))
dB186
#B187: Birds in MSF, OC_Oceanian
OCBiMSF<-subset(data,BGroupBRealm=="OC_OceanianBirdsMSF")
nrow(OCBiMSF)
B187<-rep(OCBiMSF$RRS,10000)
dB187<-c(Bmean=mean(B187),"2.5%" = NaN,"97.5%" = NaN,n=nrow(OCBiMSF))
dB187
#B188: Birds in OSF, OC_Oceanian (no data)
B188<-rep(NaN,10000)
dB188<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB188
#B189: Mammals in ES,  OC_Oceanian
B189<-rep(NaN,10000)
dB189<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB189
#B190: Mammals in YSF, OC_Oceanian
B190<-rep(NaN,10000)
dB190<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB190
#B191: Mammals in MSF, OC_Oceanian
B191<-rep(NaN,10000)
dB191<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB191
#B192: Mammals in OSF, OC_Oceanian
B192<-rep(NaN,10000)
dB192<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB192

##~~ AU_Australian Biogeographic Realm
#B193: Amphibians in Early Succession (ES), AU_Australian
AUAmES<-subset(data,BGroupBRealm=="AU_AustralianAmphibiansES")
nrow(AUAmES)
B193<-rep(AUAmES$RRS,10000)
dB193<-c(Bmean=mean(B193),"2.5%" = NaN,"97.5%" = NaN,n=nrow(AUAmES))
dB193
#B194: Amphibians in Young Secondary Forest (YSF), AU_Australian
B194<-rep(NaN,10000)
dB194<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB194
#B195: Amphibians in Mid-successioMAl Secondary Forest (MSF), AU_Australian
AUAmMSF<-subset(data,BGroupBRealm=="AU_AustralianAmphibiansMSF")
nrow(AUAmMSF)
B195<-rep(AUAmMSF$RRS,10000)
dB195<-c(Bmean=mean(B195),"2.5%" = NaN,"97.5%" = NaN,n=nrow(AUAmMSF))
dB195
#B196: Amphibians in Old Secondary Forest (OSF), AU_Australian (no data)
B196<-rep(NaN,10000)
dB196<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB196
#B197: Reptiles in ES, AU_Australian (no data)
AUReES<-subset(data,BGroupBRealm=="AU_AustralianReptilesES")
nrow(AUReES)
B197<-rep(AUReES$RRS,10000)
dB197<-c(Bmean=mean(B197),"2.5%" = NaN,"97.5%" = NaN,n=nrow(AUReES))
dB197
#B198: Reptiles in YSF, AU_Australian (no data)
B198<-rep(NaN,10000)
dB198<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB198
#B199: Reptiles in MSF, AU_Australian
AUReMSF<-subset(data,BGroupBRealm=="AU_AustralianReptilesMSF")
nrow(AUReMSF)
B199<-rep(AUReMSF$RRS,10000)
dB199<-c(Bmean=mean(B199),"2.5%" = NaN,"97.5%" = NaN,n=nrow(AUReMSF))
dB199
#B200: Reptiles in OSF, AU_Australian (there is no data)
B200<-rep(NaN,10000)
dB200<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB200
#B201: Birds in ES, AU_Australian
AUBiES<-subset(data,BGroupBRealm=="AU_AustralianBirdsES")
nrow(AUBiES)
B201<-rep(AUBiES$RRS,10000)
dB201<-c(Bmean=mean(B201),"2.5%" = NaN,"97.5%" = NaN,n=nrow(AUBiES))
dB201
#B202: Birds in YSF, AU_Australian
B202<-rep(NaN,10000)
dB202<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB202
#B203: Birds in MSF, AU_Australian
AUBiMSF<-subset(data,BGroupBRealm=="AU_AustralianBirdsMSF")
nrow(AUBiMSF)
B203<-rep(AUBiMSF$RRS,10000)
dB203<-c(Bmean=mean(B203),"2.5%" = NaN,"97.5%" = NaN,n=nrow(AUBiMSF))
dB203
#B204: Birds in OSF, AU_Australian
B204<-rep(NaN,10000)
dB204<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB204
#B205: Mammals in ES,  AU_Australian
AUMaES<-subset(data,BGroupBRealm=="AU_AustralianMammalsES")
nrow(AUMaES)
B205<-rep(AUMaES$RRS,10000)
dB205<-c(Bmean=mean(B205),"2.5%" = NaN,"97.5%" = NaN,n=nrow(AUMaES))
dB205
#B206: Mammals in YSF, AU_Australian
B206<-rep(NaN,10000)
dB206<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB206
#B207: Mammals in MSF, AU_Australian
AUMaMSF<-subset(data,BGroupBRealm=="AU_AustralianMammalsMSF")
nrow(AUMaMSF)
B207<-rep(AUMaMSF$RRS,10000)
dB207<-c(Bmean=mean(B207),"2.5%" = NaN,"97.5%" = NaN,n=nrow(AUMaMSF))
dB207
#B208: Mammals in OSF, AU_Australian (no data)
B208<-rep(NaN,10000)
dB208<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB208

##~~ MA_Madagascan Biogeographic Realm
#B209: Amphibians in Early Succession (ES), MA_Madagascan
B209<-rep(NaN,10000)
dB209<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB209
#B210: Amphibians in Young Secondary Forest (YSF), MA_Madagascan
MAAmYSF<-subset(data,BGroupBRealm=="MA_MadagascanAmphibiansYSF")
nrow(MAAmYSF)
B210<-rep(MAAmYSF$RRS,10000)
dB210<-c(Bmean=mean(B210),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MAAmYSF))
dB210
#B211: Amphibians in Mid-successioMAl Secondary Forest (MSF), MA_Madagascan
MAAmMSF<-subset(data,BGroupBRealm=="MA_MadagascanAmphibiansMSF")
nrow(MAAmMSF)
B211<-rep(MAAmMSF$RRS,10000)
dB211<-c(Bmean=mean(B211),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MAAmMSF))
dB211
#B212: Amphibians in Old Secondary Forest (OSF), MA_Madagascan (no data)
B212<-rep(NaN,10000)
dB212<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB212
#B213: Reptiles in ES, MA_Madagascan
B213<-rep(NaN,10000)
dB213<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB213
#B214: Reptiles in YSF, MA_Madagascan
B214<-rep(NaN,10000)
dB214<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB214
#B215: Reptiles in MSF, MA_Madagascan (no data)
B215<-rep(NaN,10000)
dB215<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB215
#B216: Reptiles in OSF, MA_Madagascan (there is no data)
B216<-rep(NaN,10000)
dB216<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB216
#B217: Birds in ES, MA_Madagascan 
MABiES<-subset(data,BGroupBRealm=="MA_MadagascanBirdsES")
nrow(MABiES)
B217<-rep(MABiES$RRS,10000)
dB217<-c(Bmean=mean(B217),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MABiES))
dB217
#B218: Birds in YSF, MA_Madagascan
MABiYSF<-subset(data,BGroupBRealm=="MA_MadagascanBirdsYSF")
nrow(MABiYSF)
B218<-rep(MABiYSF$RRS,10000)
dB218<-c(Bmean=mean(B218),"2.5%" = NaN,"97.5%" = NaN,n=nrow(MABiYSF))
dB218
#B219: Birds in MSF, MA_Madagascan
B219<-rep(NaN,10000)
dB219<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB219
#B220: Birds in OSF, MA_Madagascan (no data)
B220<-rep(NaN,10000)
dB220<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB220
#B221: Mammals in ES,  MA_Madagascan
B221<-rep(NaN,10000)
dB221<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB221
#B222: Mammals in YSF, MA_Madagascan
B222<-rep(NaN,10000)
dB222<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB222
#B223: Mammals in MSF, MA_Madagascan
B223<-rep(NaN,10000)
dB223<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB223
#B224: Mammals in OSF, MA_Madagascan (no data)
B224<-rep(NaN,10000)
dB224<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB224


#Join results
#Resume data of Bootstrap (Bootstrap mean, confidence limits[2.5%-97.5%],n)
dbR=data.frame(dB1,dB2,dB3,dB4,dB5,dB6,dB7,dB8,dB9,dB10, dB11,dB12,dB13,dB14,dB15,dB16,dB17,dB18,dB19,dB20, 
               dB21,dB22,dB23,dB24,dB25,dB26,dB27,dB28,dB29,dB30, dB31,dB32,dB33,dB34,dB35,dB36,dB37,dB38,dB39,dB40,
               dB41,dB42,dB43,dB44,dB45,dB46,dB47,dB48,dB49,dB50, dB51,dB52,dB53,dB54,dB55,dB56,dB57,dB58,dB59,dB60,
               dB61,dB62,dB63,dB64,dB65,dB66,dB67,dB68,dB69,dB70, dB71,dB72,dB73,dB74,dB75,dB76,dB77,dB78,dB79,dB80,
               dB81,dB82,dB83,dB84,dB85,dB86,dB87,dB88,dB89,dB90, dB91,dB92,dB93,dB94,dB95,dB96,dB97,dB98,dB99,dB100,
               dB101,dB102,dB103,dB104,dB105,dB106,dB107,dB108,dB109,dB110, dB111,dB112,dB113,dB114,dB115,dB116,dB117,dB118,dB119,dB120, 
               dB121,dB122,dB123,dB124,dB125,dB126,dB127,dB128,dB129,dB130, dB131,dB132,dB133,dB134,dB135,dB136,dB137,dB138,dB139,dB140,
               dB141,dB142,dB143,dB144,dB145,dB146,dB147,dB148,dB149,dB150, dB151,dB152,dB153,dB154,dB155,dB156,dB157,dB158,dB159,dB160,
               dB161,dB162,dB163,dB164,dB165,dB166,dB167,dB168,dB169,dB170, dB171,dB172,dB173,dB174,dB175,dB176,dB177,dB178,dB179,dB180,
               dB181,dB182,dB183,dB184,dB185,dB186,dB187,dB188,dB189,dB190, dB191,dB192,dB193,dB194,dB195,dB196,dB197,dB198,dB199,dB200,
               dB201,dB202,dB203,dB204,dB205,dB206,dB207,dB208,dB209,dB210, dB211,dB212,dB213,dB214,dB215,dB216,dB217,dB218,dB219,dB220, 
               dB221,dB222,dB223,dB224)
t.dbR<-t(dbR)
datboots<-as.data.frame(t.dbR)
head(datboots)
datboots$FigurePart<-c(rep("Overall",16),rep("Moist",16),rep("Dry",16),rep("Savannahs",16),rep("Continent",16),rep("Island",16),
                       rep("PN_Panamanian",16),rep("NT_Neotropical",16),rep("OR_Oriental",16),rep("AF_Afrotropical",16),
                       rep("NA_Neartic",16),rep("OC_Oceanian",16),rep("AU_Australian",16),rep("MA_Madagascan",16))
taxdb<-c(rep("Amphibians",4),rep("Reptiles",4),rep("Birds",4),rep("Mammals",4))
Succ.stage<-c(rep("ES",1),rep("YSF",1),rep("MSF",1),rep("OSF",1))
datboots$Taxa<-c(rep(taxdb,14))
datboots$Succ.Stage<-c(rep(Succ.stage,4))
head(datboots)
tail(datboots)
write.table(datboots, "BootstrappSppCompositionSimilarity.txt",quote=F, sep="\t")


#Lets graph all this work
colnames(datboots)<- c("Bmean","lower","upper", "n","FigurePart","Taxa", "Succ.Stage")
datboots$n

datboots$FigurePart <- factor(datboots$FigurePart, levels=c('Overall',
                                                            'Moist', 'Dry','Savannahs',
                                                            'Continent','Island',
                                                            'PN_Panamanian','NT_Neotropical',
                                                            'OR_Oriental','AF_Afrotropical',
                                                            'NA_Neartic','OC_Oceanian',
                                                            'AU_Australian','MA_Madagascan'))
datboots$Succ.Stage <- factor(datboots$Succ.Stage,levels=c('ES','YSF','MSF', 'OSF'))
datboots$Taxa <- factor(datboots$Taxa,levels=c('Mammals','Birds','Reptiles','Amphibians'))

#Subsets
#Overall effect
ov<-subset(datboots,FigurePart=="Overall")
a<-ggplot(ov, aes(y = Bmean, ymin = lower, ymax = upper,
                  x = Taxa, shape=Succ.Stage, fill=Taxa))+
  geom_vline(xintercept=c(1.5,2.5,3.5),color="darkgray", size=0.4)+
  geom_hline(yintercept = 0, linetype = "dashed",color="black", size=0.5)+
  geom_text(aes(label=ov$n, y=0.2, fill=Taxa), 
            position = position_dodge(width = 1), size=2.5)+
  #geom_pointrange(position = position_dodge(1.2), size=0.8,aes(fill=Taxa), colour="black", stroke=1)+
  geom_linerange(position = position_dodge(1), size=1.5, alpha=0.6, aes(color=Taxa))+
  geom_point(position = position_dodge(1), size=3.5, stroke = 1, alpha=0.8, aes(fill=Taxa))+
  ggtitle('a')+
  ylab("")+  
  xlab("")+
  scale_shape_manual(name="Succ.Stage",
                     values = c("ES" = 21, "YSF"=22, 
                                "MSF"=24,"OSF"=23))+
  scale_y_continuous(breaks = seq(-1.5,0.1, by=0.5),
                     labels = seq(-1.5,0.1, by=0.5),
                     limits = c(-1.75,0.3), expand = c(0, 0))+
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
  theme(legend.position="none",axis.text.x=element_text(size=12, hjust = 1), axis.text.y=element_blank(), 
        plot.margin=unit(c(1,1,-2,1),"mm"),panel.margin.y = unit(0, "lines"), panel.grid.major = element_blank(), panel.grid.minor = 	element_blank())+ 
  theme(strip.background = element_blank(),plot.title=element_text(hjust=-0.025, size = 16, face = "bold"),
        strip.text = element_blank())+
  coord_flip()+
  facet_wrap(~FigurePart, ncol = 2)+annotate("text",x=4.25,y=-1.6,label = "Overall", size=5,hjust=0)
a

#Continent vs Island
cont<-datboots[65:96,]
b<-ggplot(cont, aes(y = Bmean, ymin = lower, ymax = upper,
                    x = Taxa, shape=Succ.Stage, fill=Taxa))+
  geom_vline(xintercept=c(1.5,2.5,3.5),color="darkgray", size=0.4)+
  geom_hline(yintercept = 0, linetype = "dashed",color="black", size=0.5)+
  geom_text(aes(label=cont$n, y=0.2, fill=Taxa), 
            position = position_dodge(width = 1), size=2.5)+
  #geom_pointrange(position = position_dodge(1.2), size=0.8,aes(fill=Taxa), colour="black", stroke=1)+
  geom_linerange(position = position_dodge(1), size=1.5, alpha=0.6, aes(color=Taxa))+
  geom_point(position = position_dodge(1), size=3.5, stroke = 1, alpha=0.8, aes(fill=Taxa))+
  ggtitle('b')+
  ylab("")+  
  xlab("")+
  scale_shape_manual(name="Succ.Stage",
                     values = c("ES" = 21, "YSF"=22, 
                                "MSF"=24,"OSF"=23))+
  scale_y_continuous(breaks = seq(-1.5,0.1, by=0.5),
                     labels = seq(-1.5,0.1, by=0.5),
                     limits = c(-1.75,0.3), expand = c(0, 0))+
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
  theme(legend.position="none",axis.text.x=element_text(size=12, hjust = 1), axis.text.y=element_blank(), 
        plot.margin=unit(c(1,1,-2,1),"mm"),panel.margin = unit(0, "lines"), panel.grid.major = element_blank(), panel.grid.minor = 	element_blank())+ 
  theme(strip.background = element_blank(),plot.title=element_text(hjust=-0.025, size = 16, face = "bold"))+
  theme(strip.text = element_blank())+
  coord_flip()+
  facet_wrap(~FigurePart, ncol = 2, nrow = 1)+
  annotate("text",x=4.25,y=-1.6,label = c("Continent","Island"), size=5,hjust=0)
b

#Biomes
biom<-datboots[17:64,]
c<-ggplot(biom, aes(y = Bmean, ymin = lower, ymax = upper,
                    x = Taxa, shape=Succ.Stage, fill=Taxa))+
  geom_vline(xintercept=c(1.5,2.5,3.5),color="darkgray", size=0.4)+
  geom_hline(yintercept = 0, linetype = "dashed",color="black", size=0.5)+
  geom_text(aes(label=biom$n, y=0.2, fill=Taxa), 
            position = position_dodge(width = 1), size=2.5)+
  #geom_pointrange(position = position_dodge(1.2), size=0.8,aes(fill=Taxa), colour="black", stroke=1)+
  geom_linerange(position = position_dodge(1), size=1.5, alpha=0.6, aes(color=Taxa))+
  geom_point(position = position_dodge(1), size=3.5, stroke = 1, alpha=0.8, aes(fill=Taxa))+
  ggtitle('c')+
  ylab("")+  
  xlab("")+
  scale_shape_manual(name="Succ.Stage",
                     values = c("ES" = 21, "YSF"=22, 
                                "MSF"=24,"OSF"=23))+
  scale_y_continuous(breaks = seq(-1.5,0.1, by=0.5),
                     labels = seq(-1.5,0.1, by=0.5),
                     limits = c(-1.75,0.3), expand = c(0, 0))+
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
  theme(legend.position="none",axis.text.x=element_text(size=12, hjust = 1), axis.text.y=element_blank(), 
        plot.margin=unit(c(1,1,-2,1),"mm"),panel.margin.y = unit(0, "lines"), panel.grid.major = element_blank(), panel.grid.minor = 	element_blank())+ 
  theme(strip.background = element_blank(),plot.title=element_text(hjust=-0.025, size = 16, face = "bold"))+
  theme(strip.text = element_blank())+
  coord_flip()+
  facet_wrap(~FigurePart, ncol = 3, nrow = 1)+
  annotate("text",x=4.25,y=-1.6,label = c("Moist","Dry","Savannahs"), size=5,hjust=0)
c

#Zoogeographical regions
zoog<-datboots[97:223,]
d<-ggplot(zoog, aes(y = Bmean, ymin = lower, ymax = upper,
                    x = Taxa, shape=Succ.Stage, fill=Taxa))+
  geom_vline(xintercept=c(1.5,2.5,3.5),color="darkgray", size=0.4)+
  geom_hline(yintercept = 0, linetype = "dashed",color="black", size=0.5)+
  geom_text(aes(label=zoog$n, y=0.2, fill=Taxa), 
            position = position_dodge(width = 1), size=2.5)+
  #geom_pointrange(position = position_dodge(1.2), size=0.8,aes(fill=Taxa), colour="black", stroke=1)+
  geom_linerange(position = position_dodge(1), size=1.5, alpha=0.6, aes(color=Taxa))+
  geom_point(position = position_dodge(1), size=3.5, stroke = 1, alpha=0.8, aes(fill=Taxa))+
  ggtitle('d')+
  ylab("")+  
  xlab("")+
  scale_shape_manual(name="Succ.Stage",
                     values = c("ES" = 21, "YSF"=22, 
                                "MSF"=24,"OSF"=23))+
  scale_y_continuous(breaks = seq(-1.5,0.1, by=0.5),
                     labels = seq(-1.5,0.1, by=0.5),
                     limits = c(-1.75,0.3), expand = c(0, 0))+
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
  theme(legend.position="none",axis.text.x=element_text(size=12, hjust = 1), axis.text.y=element_blank(), 
        plot.margin=unit(c(1,1,-2,1),"mm"),panel.margin.y = unit(0, "lines"), panel.grid.major = element_blank(), panel.grid.minor = 	element_blank())+ 
  theme(strip.background = element_blank(),plot.title=element_text(hjust=-0.025, size = 16, face = "bold"))+
  theme(strip.text = element_blank())+
  coord_flip()+
  facet_wrap(~FigurePart, ncol = 4, nrow = 2)+
  annotate("text",x=4.25,y=-1.6,label = c('PN','NT','OR','AF',
                                          'NA','OC',
                                          'AU','MA'), size=5, hjust=0)
d


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Fit single model, age as single predictor for faunal recovery of species composition similarity
#Clear memory
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
  ggtitle('e')+
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
  theme(plot.margin = unit(c(1,1,1,1),"mm"))+
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
e<-lpi+ 
  geom_segment(data=B, aes(x=((0.9-upPIB$coefficients[1])/upPIB$coefficients[2]),
                           xend=((0.9-lwPIB$coefficients[1])/lwPIB$coefficients[2]), y=1.2, yend=1.2), 
               linetype="solid", colour="#1f78b4", size=2)+
  geom_point(data=B, aes(x=((0.9-SSB$coefficients[1])/SSB$coefficients[2]), 
                         xend=((0.9-SSB$coefficients[1])/SSB$coefficients[2]), y=1.2, yend=1.2), 
             colour="black", fill="black",shape=21, size=6)
e#Figure of the recovery time


#Joint figures
abc=grid.arrange(arrangeGrob(arrangeGrob(a, b, ncol=2,widths = c(1,1)),(c)))
abcd=grid.arrange(arrangeGrob(arrangeGrob(abc),(d)),
                  left=textGrob("", 
                                rot=90, vjust=1),
                  top=textGrob(""),
                  right=textGrob(""),
                  bottom=textGrob("Response ratio of the vertebrate species composition similarity during secondary forest succession \n (bootstrapped effect size)", 
                                  gp=gpar(fontsize=14)))

grid.arrange(abcd,e,nrow=2,heights = c(1.75,0.65))

