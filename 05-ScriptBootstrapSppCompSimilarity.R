#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Bootstrapping of the Response ratio in Species composition similarity
#Clear memory
rm(list=ls())
setwd("c:/Users/Flaco/Dropbox/Tesis MSc/Capitulo 1 - Review/Articulos base/Results/GEB/Submitted/Resubmission")
#Call the file MainData.csv
mydata <- read.csv("C:/Users/Flaco/Dropbox/Tesis MSc/Capitulo 1 - Review/Articulos base/Results/GEB/Submitted/Resubmission/MainData.csv")
#packages used
library(ggplot2)
library(reshape2)
library(gridExtra)
library(gtable)
library(grid)

#Exclude the forest with no list of species (the RRS can't be calculated), and the OGF data
data1<-subset(mydata,RRS!="NA")
data<-subset(data1,Forest.cov!="OGF")
nrow(data)#170
head(data)
summary(data$BGroupOverall)
#B1: Amphibians in Early Succession (ES), overall
AmES<-subset(data,BGroupOverall=="AmphibiansES")
median(AmES$RRS)
R<-10000
B1<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
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
  boot.sample <- sample(AmYSF$RRS, replace = TRUE)
  B2[kk] <- mean(boot.sample)
}
boxplot(B2)
quantile(B2,c(0.025,0.975))
hist(B2, breaks = 30)
dB2<-c(Bmean=mean(B2),quantile(B2,c(0.025,0.975)),n=nrow(AmYSF))
dB2
#B3: Amphibians in Mid-successional Secondary Forest (MSF), overall
AmMSF<-subset(data,BGroupOverall=="AmphibiansMSF")
median(AmMSF$RRS)
R<-10000
B3<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(AmMSF$RRS, replace = TRUE)
  B3[kk] <- mean(boot.sample)
}
boxplot(B3)
quantile(B3,c(0.025,0.975))
hist(B3, breaks = 30)
dB3<-c(Bmean=mean(B3),quantile(B3,c(0.025,0.975)),n=nrow(AmMSF))
dB3
#B4: Amphibians in Old Secondary Forest (OSF), overall (there is no data)
B4<-rep(NaN,10000)
dB4<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB4
#B5: Reptiles in ES, overall
ReES<-subset(data,BGroupOverall=="ReptilesES")
median(ReES$RRS)
R<-10000
B5<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(ReES$RRS, replace = TRUE)
  B5[kk] <- mean(boot.sample)
}
boxplot(B5)
quantile(B5,c(0.025,0.975))
hist(B5, breaks = 30)
dB5<-c(Bmean=mean(B5),quantile(B5,c(0.025,0.975)),n=nrow(ReES))
dB5
#B6: Reptiles in YSF, overall
ReYSF<-subset(data,BGroupOverall=="ReptilesYSF")
median(ReYSF$RRS)
R<-10000
B6<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(ReYSF$RRS, replace = TRUE)
  B6[kk] <- mean(boot.sample)
}
boxplot(B6)
quantile(B6,c(0.025,0.975))
hist(B6, breaks = 30)
dB6<-c(Bmean=mean(B6),quantile(B6,c(0.025,0.975)),n=nrow(ReYSF))
dB6
#B7: Reptiles in MSF, overall
ReMSF<-subset(data,BGroupOverall=="ReptilesMSF")
median(ReMSF$RRS)
R<-10000
B7<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(ReMSF$RRS, replace = TRUE)
  B7[kk] <- mean(boot.sample)
}
boxplot(B7)
quantile(B7,c(0.025,0.975))
hist(B7, breaks = 30)
dB7<-c(Bmean=mean(B7),quantile(B7,c(0.025,0.975)),n=nrow(ReMSF))
dB7
#B8: Reptiles in OSF, overall
ReOSF<-subset(data,BGroupOverall=="ReptilesOSF")
median(ReOSF$RRS)
nrow(ReOSF)
R<-10000
B8<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(ReOSF$RRS, replace = TRUE)
  B8[kk] <- mean(boot.sample)
}
boxplot(B8)
quantile(B8,c(0.025,0.975))
hist(B8, breaks = 30)
dB8<-c(Bmean=mean(B8),quantile(B8,c(0.025,0.975)),n=nrow(ReOSF))
dB8
#B9: Birds in ES, overall
BiES<-subset(data,BGroupOverall=="BirdsES")
median(BiES$RRS)
nrow(BiES)
R<-10000
B9<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(BiES$RRS, replace = TRUE)
  B9[kk] <- mean(boot.sample)
}
boxplot(B9)
quantile(B9,c(0.025,0.975))
hist(B9, breaks = 30)
dB9<-c(Bmean=mean(B9),quantile(B9,c(0.025,0.975)),n=nrow(BiES))
dB9
#B10: Birds in YSF, overall
BiYSF<-subset(data,BGroupOverall=="BirdsYSF")
median(BiYSF$RRS)
nrow(BiYSF)
R<-10000
B10<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(BiYSF$RRS, replace = TRUE)
  B10[kk] <- mean(boot.sample)
}
boxplot(B10)
quantile(B10,c(0.025,0.975))
hist(B10, breaks = 30)
dB10<-c(Bmean=mean(B10),quantile(B10,c(0.025,0.975)),n=nrow(BiYSF))
dB10
#B11: Birds in MSF, overall
BiMSF<-subset(data,BGroupOverall=="BirdsMSF")
median(BiMSF$RRS)
nrow(BiMSF)
R<-10000
B11<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(BiMSF$RRS, replace = TRUE)
  B11[kk] <- mean(boot.sample)
}
boxplot(B11)
quantile(B11,c(0.025,0.975))
hist(B11, breaks = 30)
dB11<-c(Bmean=mean(B11),quantile(B11,c(0.025,0.975)),n=nrow(BiMSF))
dB11
#B12: Birds in OSF, overall (no data)
B12<-rep(NaN,10000)
dB12<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB12
#B13: Mammals in ES, overall
MaES<-subset(data,BGroupOverall=="MammalsES")
median(MaES$RRS)
nrow(MaES)
R<-10000
B13<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(MaES$RRS, replace = TRUE)
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
  boot.sample <- sample(MaYSF$RRS, replace = TRUE)
  B14[kk] <- mean(boot.sample)
}
boxplot(B14)
quantile(B14,c(0.025,0.975))
hist(B14, breaks = 30)
dB14<-c(Bmean=mean(B14),quantile(B14,c(0.025,0.975)),n=nrow(MaYSF))
dB14
#B15: Mammals in MSF, overall
MaMSF<-subset(data,BGroupOverall=="MammalsMSF")
median(MaMSF$RRS)
nrow(MaMSF)
R<-10000
B15<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(MaMSF$RRS, replace = TRUE)
  B15[kk] <- mean(boot.sample)
}
boxplot(B15)
quantile(B15,c(0.025,0.975))
hist(B15, breaks = 30)
dB15<-c(Bmean=mean(B15),quantile(B15,c(0.025,0.975)),n=nrow(MaMSF))
dB15
#B16: Mammals in OSF, overall
MaOSF<-subset(data,BGroupOverall=="MammalsOSF")
median(MaOSF$RRS)
nrow(MaOSF)
R<-10000
B16<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(MaOSF$RRS, replace = TRUE)
  B16[kk] <- mean(boot.sample)
}
boxplot(B16)
quantile(B16,c(0.025,0.975))
hist(B16, breaks = 30)
dB16<-c(Bmean=mean(B16),quantile(B16,c(0.025,0.975)),n=nrow(MaOSF))
dB16

##~~ Neotropical Biogeographic Realm
#B17: Amphibians in Early Succession (ES), Neotropic
NAmES<-subset(data,BGroupBRealm=="NeotropicAmphibiansES")
median(NAmES$RRS)
R<-10000
B17<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(NAmES$RRS, replace = TRUE)
  B17[kk] <- mean(boot.sample)
}
boxplot(B17)
quantile(B17,c(0.025,0.975))
hist(B17, breaks = 30)
dB17<-c(Bmean=mean(B17),quantile(B17,c(0.025,0.975)),n=nrow(NAmES))
dB17
#B18: Amphibians in Young Secondary Forest (YSF), Neotropic
NAmYSF<-subset(data,BGroupBRealm=="NeotropicAmphibiansYSF")
median(NAmYSF$RRS)
nrow(NAmYSF)
R<-10000
B18<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(NAmYSF$RRS, replace = TRUE)
  B18[kk] <- mean(boot.sample)
}
boxplot(B18)
quantile(B18,c(0.025,0.975))
hist(B18, breaks = 30)
dB18<-c(Bmean=mean(B18),quantile(B18,c(0.025,0.975)),n=nrow(NAmYSF))
dB18
#B19: Amphibians in Mid-successional Secondary Forest (MSF), Neotropic
NAmMSF<-subset(data,BGroupBRealm=="NeotropicAmphibiansMSF")
median(NAmMSF$RRS)
R<-10000
B19<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(NAmMSF$RRS, replace = TRUE)
  B19[kk] <- mean(boot.sample)
}
boxplot(B19)
quantile(B19,c(0.025,0.975))
hist(B19, breaks = 30)
dB19<-c(Bmean=mean(B19),quantile(B19,c(0.025,0.975)),n=nrow(NAmMSF))
dB19
#B20: Amphibians in Old Secondary Forest (OSF), Neotropic (there is no data)
B20<-rep(NaN,10000)
dB20<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB20
#B21: Reptiles in ES, Neotropic
NReES<-subset(data,BGroupBRealm=="NeotropicReptilesES")
median(NReES$RRS)
R<-10000
B21<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(NReES$RRS, replace = TRUE)
  B21[kk] <- mean(boot.sample)
}
boxplot(B21)
quantile(B21,c(0.025,0.975))
hist(B21, breaks = 30)
dB21<-c(Bmean=mean(B21),quantile(B21,c(0.025,0.975)),n=nrow(NReES))
dB21
#B22: Reptiles in YSF, Neotropic
NReYSF<-subset(data,BGroupBRealm=="NeotropicReptilesYSF")
median(NReYSF$RRS)
R<-10000
B22<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(NReYSF$RRS, replace = TRUE)
  B22[kk] <- mean(boot.sample)
}
boxplot(B22)
quantile(B22,c(0.025,0.975))
hist(B22, breaks = 30)
dB22<-c(Bmean=mean(B22),quantile(B22,c(0.025,0.975)),n=nrow(NReYSF))
dB22
#B23: Reptiles in MSF, Neotropic
NReMSF<-subset(data,BGroupBRealm=="NeotropicReptilesMSF")
median(NReMSF$RRS)
R<-10000
B23<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(NReMSF$RRS, replace = TRUE)
  B23[kk] <- mean(boot.sample)
}
boxplot(B23)
quantile(B23,c(0.025,0.975))
hist(B23, breaks = 30)
dB23<-c(Bmean=mean(B23),quantile(B23,c(0.025,0.975)),n=nrow(NReMSF))
dB23
#B24: Reptiles in OSF, Neotropic
NReOSF<-subset(data,BGroupBRealm=="NeotropicReptilesOSF")
median(NReOSF$RRS)
nrow(NReOSF)
R<-10000
B24<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(NReOSF$RRS, replace = TRUE)
  B24[kk] <- mean(boot.sample)
}
boxplot(B24)
quantile(B24,c(0.025,0.975))
hist(B24, breaks = 30)
dB24<-c(Bmean=mean(B24),quantile(B24,c(0.025,0.975)),n=nrow(NReOSF))
dB24
#B25: Birds in ES, Neotropic
NBiES<-subset(data,BGroupBRealm=="NeotropicBirdsES")
median(NBiES$RRS)
nrow(NBiES)
R<-10000
B25<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(NBiES$RRS, replace = TRUE)
  B25[kk] <- mean(boot.sample)
}
boxplot(B25)
quantile(B25,c(0.025,0.975))
hist(B25, breaks = 30)
dB25<-c(Bmean=mean(B25),quantile(B25,c(0.025,0.975)),n=nrow(NBiES))
dB25
#B26: Birds in YSF, Neotropic
NBiYSF<-subset(data,BGroupBRealm=="NeotropicBirdsYSF")
median(NBiYSF$RRS)
nrow(NBiYSF)
R<-10000
B26<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(NBiYSF$RRS, replace = TRUE)
  B26[kk] <- mean(boot.sample)
}
boxplot(B26)
quantile(B26,c(0.025,0.975))
hist(B26, breaks = 30)
dB26<-c(Bmean=mean(B26),quantile(B26,c(0.025,0.975)),n=nrow(NBiYSF))
dB26
#B27: Birds in MSF, Neotropic
NBiMSF<-subset(data,BGroupBRealm=="NeotropicBirdsMSF")
median(NBiMSF$RRS)
nrow(NBiMSF)
R<-10000
B27<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(NBiMSF$RRS, replace = TRUE)
  B27[kk] <- mean(boot.sample)
}
boxplot(B27)
quantile(B27,c(0.025,0.975))
hist(B27, breaks = 30)
dB27<-c(Bmean=mean(B27),quantile(B27,c(0.025,0.975)),n=nrow(NBiMSF))
dB27
#B28: Birds in OSF, Neotropic (no data)
B28<-rep(NaN,10000)
dB28<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB28
#B29: Mammals in ES, Neotropic
NMaES<-subset(data,BGroupBRealm=="NeotropicMammalsES")
median(NMaES$RRS)
nrow(NMaES)
R<-10000
B29<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(NMaES$RRS, replace = TRUE)
  B29[kk] <- mean(boot.sample)
}
boxplot(B29)
quantile(B29,c(0.025,0.975))
hist(B29, breaks = 30)
dB29<-c(Bmean=mean(B29),quantile(B29,c(0.025,0.975)),n=nrow(NMaES))
dB29
#B30: Mammals in YSF, Neotropic
NMaYSF<-subset(data,BGroupBRealm=="NeotropicMammalsYSF")
median(NMaYSF$RRS)
nrow(NMaYSF)
R<-10000
B30<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(NMaYSF$RRS, replace = TRUE)
  B30[kk] <- mean(boot.sample)
}
boxplot(B30)
quantile(B30,c(0.025,0.975))
hist(B30, breaks = 30)
dB30<-c(Bmean=mean(B30),quantile(B30,c(0.025,0.975)),n=nrow(NMaYSF))
dB30
#B31: Mammals in MSF, Neotropic
NMaMSF<-subset(data,BGroupBRealm=="NeotropicMammalsMSF")
median(NMaMSF$RRS)
nrow(NMaMSF)
R<-10000
B31<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(NMaMSF$RRS, replace = TRUE)
  B31[kk] <- mean(boot.sample)
}
boxplot(B31)
quantile(B31,c(0.025,0.975))
hist(B31, breaks = 30)
dB31<-c(Bmean=mean(B31),quantile(B31,c(0.025,0.975)),n=nrow(NMaMSF))
dB31
#B32: Mammals in OSF, Neotropic
NMaOSF<-subset(data,BGroupBRealm=="NeotropicMammalsOSF")
median(NMaOSF$RRS)
nrow(NMaOSF)
R<-10000
B32<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(NMaOSF$RRS, replace = TRUE)
  B32[kk] <- mean(boot.sample)
}
boxplot(B32)
quantile(B32,c(0.025,0.975))
hist(B32, breaks = 30)
dB32<-c(Bmean=mean(B32),quantile(B32,c(0.025,0.975)),n=nrow(NMaOSF))
dB32

##~~ Australasia Biogeographic Realm
#B33: Amphibians in Early Succession (ES), Australasia
AuAmES<-subset(data,BGroupBRealm=="AustralasiaAmphibiansES")
median(AuAmES$RRS)
R<-10000
B33<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(AuAmES$RRS, replace = TRUE)
  B33[kk] <- mean(boot.sample)
}
boxplot(B33)
quantile(B33,c(0.025,0.975))
hist(B33, breaks = 30)
dB33<-c(Bmean=mean(B33),quantile(B33,c(0.025,0.975)),n=nrow(AuAmES))
dB33
#B34: Amphibians in Young Secondary Forest (YSF), Australasia (there is no data)
B34<-rep(NaN,10000)
dB34<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB34
#B35: Amphibians in Mid-successional Secondary Forest (MSF), Australasia
AuAmMSF<-subset(data,BGroupBRealm=="AustralasiaAmphibiansMSF")
median(AuAmMSF$RRS)
R<-10000
B35<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(AuAmMSF$RRS, replace = TRUE)
  B35[kk] <- mean(boot.sample)
}
boxplot(B35)
quantile(B35,c(0.025,0.975))
hist(B35, breaks = 30)
dB35<-c(Bmean=mean(B35),quantile(B35,c(0.025,0.975)),n=nrow(AuAmMSF))
dB35
#B36: Amphibians in Old Secondary Forest (OSF), Australasia (there is no data)
B36<-rep(NaN,10000)
dB36<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB36
#B37: Reptiles in ES, Australasia
AuReES<-subset(data,BGroupBRealm=="AustralasiaReptilesES")
median(AuReES$RRS)
R<-10000
B37<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(AuReES$RRS, replace = TRUE)
  B37[kk] <- mean(boot.sample)
}
boxplot(B37)
quantile(B37,c(0.025,0.975))
hist(B37, breaks = 30)
dB37<-c(Bmean=mean(B37),quantile(B37,c(0.025,0.975)),n=nrow(AuReES))
dB37
#B38: Reptiles in YSF, Australasia
AuReYSF<-subset(data,BGroupBRealm=="AustralasiaReptilesYSF")
median(AuReYSF$RRS)
R<-10000
B38<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(AuReYSF$RRS, replace = TRUE)
  B38[kk] <- mean(boot.sample)
}
boxplot(B38)
quantile(B38,c(0.025,0.975))
hist(B38, breaks = 30)
dB38<-c(Bmean=mean(B38),quantile(B38,c(0.025,0.975)),n=nrow(AuReYSF))
dB38
#B39: Reptiles in MSF, Australasia
AuReMSF<-subset(data,BGroupBRealm=="AustralasiaReptilesMSF")
median(AuReMSF$RRS)
R<-10000
B39<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(AuReMSF$RRS, replace = TRUE)
  B39[kk] <- mean(boot.sample)
}
boxplot(B39)
quantile(B39,c(0.025,0.975))
hist(B39, breaks = 30)
dB39<-c(Bmean=mean(B39),quantile(B39,c(0.025,0.975)),n=nrow(AuReMSF))
dB39
#B40: Reptiles in OSF, Australasia (there is no data)
B40<-rep(NaN,10000)
dB40<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB40
#B41: Birds in ES, Australasia
AuBiES<-subset(data,BGroupBRealm=="AustralasiaBirdsES")
median(AuBiES$RRS)
nrow(AuBiES)
R<-10000
B41<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(AuBiES$RRS, replace = TRUE)
  B41[kk] <- mean(boot.sample)
}
boxplot(B41)
quantile(B41,c(0.025,0.975))
hist(B41, breaks = 30)
dB41<-c(Bmean=mean(B41),quantile(B41,c(0.025,0.975)),n=nrow(AuBiES))
dB41
#B42: Birds in YSF, Australasia
AuBiYSF<-subset(data,BGroupBRealm=="AustralasiaBirdsYSF")
median(AuBiYSF$RRS)
nrow(AuBiYSF)
R<-10000
B42<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(AuBiYSF$RRS, replace = TRUE)
  B42[kk] <- mean(boot.sample)
}
boxplot(B42)
quantile(B42,c(0.025,0.975))
hist(B42, breaks = 30)
dB42<-c(Bmean=mean(B42),quantile(B42,c(0.025,0.975)),n=nrow(AuBiYSF))
dB42
#B43: Birds in MSF, Australasia
AuBiMSF<-subset(data,BGroupBRealm=="AustralasiaBirdsMSF")
median(AuBiMSF$RRS)
nrow(AuBiMSF)
R<-10000
B43<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(AuBiMSF$RRS, replace = TRUE)
  B43[kk] <- mean(boot.sample)
}
boxplot(B43)
quantile(B43,c(0.025,0.975))
hist(B43, breaks = 30)
dB43<-c(Bmean=mean(B43),quantile(B43,c(0.025,0.975)),n=nrow(AuBiMSF))
dB43
#B44: Birds in OSF, Australasia (no data)
B44<-rep(NaN,10000)
dB44<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB44
#B45: Mammals in ES, Australasia
AuMaES<-subset(data,BGroupBRealm=="AustralasiaMammalsES")
median(AuMaES$RRS)
nrow(AuMaES)
R<-10000
B45<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(AuMaES$RRS, replace = TRUE)
  B45[kk] <- mean(boot.sample)
}
boxplot(B45)
quantile(B45,c(0.025,0.975))
hist(B45, breaks = 30)
dB45<-c(Bmean=mean(B45),quantile(B45,c(0.025,0.975)),n=nrow(AuMaES))
dB45
#B46: Mammals in YSF, Australasia (no data)
B46<-rep(NaN,10000)
dB46<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB46
#B47: Mammals in MSF, Australasia
AuMaMSF<-subset(data,BGroupBRealm=="AustralasiaMammalsMSF")
median(AuMaMSF$RRS)
nrow(AuMaMSF)
R<-10000
B47<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(AuMaMSF$RRS, replace = TRUE)
  B47[kk] <- mean(boot.sample)
}
boxplot(B47)
quantile(B47,c(0.025,0.975))
hist(B47, breaks = 30)
dB47<-c(Bmean=mean(B47),quantile(B47,c(0.025,0.975)),n=nrow(AuMaMSF))
dB47
#B48: Mammals in OSF, Australasia (no data)
B48<-rep(NaN,10000)
dB48<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB48


##~~ Indo Malay Biogeographic Realm
#B49: Amphibians in Early Succession (ES), Indo Malay (no data)
B49<-rep(NaN,10000)
dB49<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB49
#B50: Amphibians in Young Secondary Forest (YSF), Indo Malay
IMAmYSF<-subset(data,BGroupBRealm=="Indo_MalayAmphibiansYSF")
median(IMAmYSF$RRS)
R<-10000
B50<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(IMAmYSF$RRS, replace = TRUE)
  B50[kk] <- mean(boot.sample)
}
boxplot(B50)
quantile(B50,c(0.025,0.975))
hist(B50, breaks = 30)
dB50<-c(Bmean=mean(B50),quantile(B50,c(0.025,0.975)),n=nrow(IMAmYSF))
dB50
#B51: Amphibians in Mid-successional Secondary Forest (MSF), Indo Malay
IMAmMSF<-subset(data,BGroupBRealm=="Indo_MalayAmphibiansMSF")
median(IMAmMSF$RRS)
R<-10000
B51<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(IMAmMSF$RRS, replace = TRUE)
  B51[kk] <- mean(boot.sample)
}
boxplot(B51)
quantile(B51,c(0.025,0.975))
hist(B51, breaks = 30)
dB51<-c(Bmean=mean(B51),quantile(B51,c(0.025,0.975)),n=nrow(IMAmMSF))
dB51
#B52: Amphibians in Old Secondary Forest (OSF), Indo Malay (no data)
B52<-rep(NaN,10000)
dB52<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB52
#B53: Reptiles in ES, Indo Malay (no data)
B53<-rep(NaN,10000)
dB53<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB53
#B54: Reptiles in YSF, Indo Malay
IMReYSF<-subset(data,BGroupBRealm=="Indo_MalayReptilesYSF")
median(IMReYSF$RRS)
R<-10000
B54<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(IMReYSF$RRS, replace = TRUE)
  B54[kk] <- mean(boot.sample)
}
boxplot(B54)
quantile(B54,c(0.025,0.975))
hist(B54, breaks = 30)
dB54<-c(Bmean=mean(B54),quantile(B54,c(0.025,0.975)),n=nrow(IMReYSF))
dB54
#B55: Reptiles in MSF, Indo Malay
IMReMSF<-subset(data,BGroupBRealm=="Indo_MalayReptilesMSF")
median(IMReMSF$RRS)
R<-10000
B55<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(IMReMSF$RRS, replace = TRUE)
  B55[kk] <- mean(boot.sample)
}
boxplot(B55)
quantile(B55,c(0.025,0.975))
hist(B55, breaks = 30)
dB55<-c(Bmean=mean(B55),quantile(B55,c(0.025,0.975)),n=nrow(IMReMSF))
dB55
#B56: Reptiles in OSF, Indo Malay (there is no data)
B56<-rep(NaN,10000)
dB56<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB56
#B57: Birds in ES, Indo Malay (there is no data)
B57<-rep(NaN,10000)
dB57<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB57
#B58: Birds in YSF, Indo Malay
IMBiYSF<-subset(data,BGroupBRealm=="Indo_MalayBirdsYSF")
median(IMBiYSF$RRS)
nrow(IMBiYSF)
R<-10000
B58<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(IMBiYSF$RRS, replace = TRUE)
  B58[kk] <- mean(boot.sample)
}
boxplot(B58)
quantile(B58,c(0.025,0.975))
hist(B58, breaks = 30)
dB58<-c(Bmean=mean(B58),quantile(B58,c(0.025,0.975)),n=nrow(IMBiYSF))
dB58
#B59: Birds in MSF, Indo Malay
IMBiMSF<-subset(data,BGroupBRealm=="Indo_MalayBirdsMSF")
median(IMBiMSF$RRS)
nrow(IMBiMSF)
R<-10000
B59<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(IMBiMSF$RRS, replace = TRUE)
  B59[kk] <- mean(boot.sample)
}
boxplot(B59)
quantile(B59,c(0.025,0.975))
hist(B59, breaks = 30)
dB59<-c(Bmean=mean(B59),quantile(B59,c(0.025,0.975)),n=nrow(IMBiMSF))
dB59
#B60: Birds in OSF, Indo Malay (no data)
B60<-rep(NaN,10000)
dB60<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB60
#B61: Mammals in ES,  Indo Malay (no data)
B61<-rep(NaN,10000)
dB61<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB61
#B62: Mammals in YSF, Indo Malay
IMMaYSF<-subset(data,BGroupBRealm=="Indo_MalayMammalsYSF")
median(IMMaYSF$RRS)
nrow(IMMaYSF)
R<-10000
B62<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(IMMaYSF$RRS, replace = TRUE)
  B62[kk] <- mean(boot.sample)
}
boxplot(B62)
quantile(B62,c(0.025,0.975))
hist(B62, breaks = 30)
dB62<-c(Bmean=mean(B62),quantile(B62,c(0.025,0.975)),n=nrow(IMMaYSF))
dB62
#B63: Mammals in MSF, Indo Malay
IMMaMSF<-subset(data,BGroupBRealm=="Indo_MalayMammalsMSF")
median(IMMaMSF$RRS)
nrow(IMMaMSF)
R<-10000
B63<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(IMMaMSF$RRS, replace = TRUE)
  B63[kk] <- mean(boot.sample)
}
boxplot(B63)
quantile(B63,c(0.025,0.975))
hist(B63, breaks = 30)
dB63<-c(Bmean=mean(B63),quantile(B63,c(0.025,0.975)),n=nrow(IMMaMSF))
dB63
#B64: Mammals in OSF, Australasia (no data)
B64<-rep(NaN,10000)
dB64<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB64

##~~ Afrotropic Biogeographic Realm
#B65: Amphibians in Early Succession (ES), Afrotropic (no data)
B65<-rep(NaN,10000)
dB65<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB65
#B66: Amphibians in Young Secondary Forest (YSF), Afrotropic
AfAmYSF<-subset(data,BGroupBRealm=="AfrotropicAmphibiansYSF")
median(AfAmYSF$RRS)
R<-10000
B66<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(AfAmYSF$RRS, replace = TRUE)
  B66[kk] <- mean(boot.sample)
}
boxplot(B66)
quantile(B66,c(0.025,0.975))
hist(B66, breaks = 30)
dB66<-c(Bmean=mean(B66),quantile(B66,c(0.025,0.975)),n=nrow(AfAmYSF))
dB66
#B67: Amphibians in Mid-successional Secondary Forest (MSF), Afrotropic
AfAmMSF<-subset(data,BGroupBRealm=="AfrotropicAmphibiansMSF")
median(AfAmMSF$RRS)
R<-10000
B67<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(AfAmMSF$RRS, replace = TRUE)
  B67[kk] <- mean(boot.sample)
}
boxplot(B67)
quantile(B67,c(0.025,0.975))
hist(B67, breaks = 30)
dB67<-c(Bmean=mean(B67),quantile(B67,c(0.025,0.975)),n=nrow(AfAmMSF))
dB67
#B68: Amphibians in Old Secondary Forest (OSF), Afrotropic (no data)
B68<-rep(NaN,10000)
dB68<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB68
#B69: Reptiles in ES, Afrotropic (no data)
B69<-rep(NaN,10000)
dB69<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB69
#B70: Reptiles in YSF, Afrotropic (no data)
B70<-rep(NaN,10000)
dB70<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB70
#B71: Reptiles in MSF, Afrotropic (no data)
B71<-rep(NaN,10000)
dB71<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB71
#B72: Reptiles in OSF, Afrotropic (there is no data)
B72<-rep(NaN,10000)
dB72<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB72
#B73: Birds in ES, Afrotropic
AfBiES<-subset(data,BGroupBRealm=="AfrotropicBirdsES")
median(AfBiES$RRS)
nrow(AfBiES)
R<-10000
B73<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(AfBiES$RRS, replace = TRUE)
  B73[kk] <- mean(boot.sample)
}
boxplot(B73)
quantile(B73,c(0.025,0.975))
hist(B73, breaks = 30)
dB73<-c(Bmean=mean(B73),quantile(B73,c(0.025,0.975)),n=nrow(AfBiES))
dB73
#B74: Birds in YSF, Afrotropic
AfBiYSF<-subset(data,BGroupBRealm=="AfrotropicBirdsYSF")
median(AfBiYSF$RRS)
nrow(AfBiYSF)
R<-10000
B74<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(AfBiYSF$RRS, replace = TRUE)
  B74[kk] <- mean(boot.sample)
}
boxplot(B74)
quantile(B74,c(0.025,0.975))
hist(B74, breaks = 30)
dB74<-c(Bmean=mean(B74),quantile(B74,c(0.025,0.975)),n=nrow(AfBiYSF))
dB74
#B75: Birds in MSF, Afrotropic (no data)
B75<-rep(NaN,10000)
dB75<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB75
#B76: Birds in OSF, Afrotropic
B76<-rep(NaN,10000)
dB76<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB76
#B77: Mammals in ES,  Afrotropic (no data)
B77<-rep(NaN,10000)
dB77<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB77
#B78: Mammals in YSF, Afrotropic
AfMaYSF<-subset(data,BGroupBRealm=="AfrotropicMammalsYSF")
median(AfMaYSF$RRS)
nrow(AfMaYSF)
R<-10000
B78<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(AfMaYSF$RRS, replace = TRUE)
  B78[kk] <- mean(boot.sample)
}
boxplot(B78)
quantile(B78,c(0.025,0.975))
hist(B78, breaks = 30)
dB78<-c(Bmean=mean(B78),quantile(B78,c(0.025,0.975)),n=nrow(AfMaYSF))
dB78
#B79: Mammals in MSF, Afrotropic (no data)
B79<-rep(NaN,10000)
dB79<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB79
#B80: Mammals in OSF, Afrotropic (no data)
B80<-rep(NaN,10000)
dB80<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB80

#B81: Amphibians in Early Succession (ES), Tropical Moist Forest
MoAmES<-subset(data,BGroupBiome=="1.TropMoAmphibiansES")
median(MoAmES$RRS)
R<-10000
B81<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(MoAmES$RRS, replace = TRUE)
  B81[kk] <- mean(boot.sample)
}
boxplot(B81)
quantile(B81,c(0.025,0.975))
hist(B81, breaks = 30)
dB81<-c(Bmean=mean(B81),quantile(B81,c(0.025,0.975)),n=nrow(MoAmES))
dB81
#B82: Amphibians in Young Secondary Forest (YSF), Tropical Moist Forest
MoAmYSF<-subset(data,BGroupBiome=="1.TropMoAmphibiansYSF")
median(MoAmYSF$RRS)
nrow(MoAmYSF)
R<-10000
B82<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(MoAmYSF$RRS, replace = TRUE)
  B82[kk] <- mean(boot.sample)
}
boxplot(B82)
quantile(B82,c(0.025,0.975))
hist(B82, breaks = 30)
dB82<-c(Bmean=mean(B82),quantile(B82,c(0.025,0.975)),n=nrow(MoAmYSF))
dB82
#B83: Amphibians in Mid-successional Secondary Forest (MSF), Tropical Moist Forest
MoAmMSF<-subset(data,BGroupBiome=="1.TropMoAmphibiansMSF")
median(MoAmMSF$RRS)
R<-10000
B83<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(MoAmMSF$RRS, replace = TRUE)
  B83[kk] <- mean(boot.sample)
}
boxplot(B83)
quantile(B83,c(0.025,0.975))
hist(B83, breaks = 30)
dB83<-c(Bmean=mean(B83),quantile(B83,c(0.025,0.975)),n=nrow(MoAmMSF))
dB83
#B84: Amphibians in Old Secondary Forest (OSF), overall (there is no data)
B84<-rep(NaN,10000)
dB84<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB84
#B85: Reptiles in ES, Tropical Moist Forest
MoReES<-subset(data,BGroupBiome=="1.TropMoReptilesES")
median(MoReES$RRS)
R<-10000
B85<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(MoReES$RRS, replace = TRUE)
  B85[kk] <- mean(boot.sample)
}
boxplot(B85)
quantile(B85,c(0.025,0.975))
hist(B85, breaks = 30)
dB85<-c(Bmean=mean(B85),quantile(B85,c(0.025,0.975)),n=nrow(MoReES))
dB85
#B66: Reptiles in YSF, Tropical Moist Forest
MoReYSF<-subset(data,BGroupBiome=="1.TropMoReptilesYSF")
median(MoReYSF$RRS)
R<-10000
B86<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(MoReYSF$RRS, replace = TRUE)
  B86[kk] <- mean(boot.sample)
}
boxplot(B86)
quantile(B86,c(0.025,0.975))
hist(B86, breaks = 30)
dB86<-c(Bmean=mean(B86),quantile(B86,c(0.025,0.975)),n=nrow(MoReYSF))
dB86
#B87: Reptiles in MSF, Moist Forest
MoReMSF<-subset(data,BGroupBiome=="1.TropMoReptilesMSF")
median(MoReMSF$RRS)
R<-10000
B87<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(MoReMSF$RRS, replace = TRUE)
  B87[kk] <- mean(boot.sample)
}
boxplot(B87)
quantile(B87,c(0.025,0.975))
hist(B87, breaks = 30)
dB87<-c(Bmean=mean(B87),quantile(B87,c(0.025,0.975)),n=nrow(MoReMSF))
dB87
#B88: Reptiles in OSF, Moist
MoReOSF<-subset(data,BGroupBiome=="1.TropMoReptilesOSF")
median(MoReOSF$RRS)
nrow(MoReOSF)
R<-10000
B88<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(MoReOSF$RRS, replace = TRUE)
  B88[kk] <- mean(boot.sample)
}
boxplot(B88)
quantile(B88,c(0.025,0.975))
hist(B88, breaks = 30)
dB88<-c(Bmean=mean(B88),quantile(B88,c(0.025,0.975)),n=nrow(MoReOSF))
dB88
#B89: Birds in ES, Moist
MoBiES<-subset(data,BGroupBiome=="1.TropMoBirdsES")
median(MoBiES$RRS)
nrow(MoBiES)
R<-10000
B89<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(MoBiES$RRS, replace = TRUE)
  B89[kk] <- mean(boot.sample)
}
boxplot(B89)
quantile(B89,c(0.025,0.975))
hist(B89, breaks = 30)
dB89<-c(Bmean=mean(B89),quantile(B89,c(0.025,0.975)),n=nrow(BiES))
dB89
#B90: Birds in YSF, Moist
MoBiYSF<-subset(data,BGroupBiome=="1.TropMoBirdsYSF")
median(MoBiYSF$RRS)
nrow(MoBiYSF)
R<-10000
B90<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(MoBiYSF$RRS, replace = TRUE)
  B90[kk] <- mean(boot.sample)
}
boxplot(B90)
quantile(B90,c(0.025,0.975))
hist(B90, breaks = 30)
dB90<-c(Bmean=mean(B90),quantile(B90,c(0.025,0.975)),n=nrow(MoBiYSF))
dB90
#B91: Birds in MSF, Moist
MoBiMSF<-subset(data,BGroupBiome=="1.TropMoBirdsMSF")
median(MoBiMSF$RRS)
nrow(MoBiMSF)
R<-10000
B91<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(MoBiMSF$RRS, replace = TRUE)
  B91[kk] <- mean(boot.sample)
}
boxplot(B91)
quantile(B91,c(0.025,0.975))
hist(B91, breaks = 30)
dB91<-c(Bmean=mean(B91),quantile(B91,c(0.025,0.975)),n=nrow(MoBiMSF))
dB91
#B92: Birds in OSF, Moist
B92<-rep(NaN,10000)
dB92<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB92
#B93: Mammals in ES, Moist
MoMaES<-subset(data,BGroupBiome=="1.TropMoMammalsES")
median(MoMaES$RRS)
nrow(MoMaES)
R<-10000
B93<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(MoMaES$RRS, replace = TRUE)
  B93[kk] <- mean(boot.sample)
}
boxplot(B93)
quantile(B93,c(0.025,0.975))
hist(B93, breaks = 30)
dB93<-c(Bmean=mean(B93),quantile(B93,c(0.025,0.975)),n=nrow(MoMaES))
dB93
#B94: Mammals in YSF, Moist
MoMaYSF<-subset(data,BGroupBiome=="1.TropMoMammalsYSF")
median(MoMaYSF$RRS)
nrow(MoMaYSF)
R<-10000
B94<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(MoMaYSF$RRS, replace = TRUE)
  B94[kk] <- mean(boot.sample)
}
boxplot(B94)
quantile(B94,c(0.025,0.975))
hist(B94, breaks = 30)
dB94<-c(Bmean=mean(B94),quantile(B94,c(0.025,0.975)),n=nrow(MoMaYSF))
dB94
#B95: Mammals in MSF, Moist
MoMaMSF<-subset(data,BGroupBiome=="1.TropMoMammalsMSF")
median(MoMaMSF$RRS)
nrow(MoMaMSF)
R<-10000
B95<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(MoMaMSF$RRS, replace = TRUE)
  B95[kk] <- mean(boot.sample)
}
boxplot(B95)
quantile(B95,c(0.025,0.975))
hist(B95, breaks = 30)
dB95<-c(Bmean=mean(B95),quantile(B95,c(0.025,0.975)),n=nrow(MoMaMSF))
dB95
#B96: Mammals in OSF, Moist
B96<-rep(NaN,10000)
dB96<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB96

##~~ Tropical Dry Forest
#B97: Amphibians in Early Succession (ES), Tropical Dry Forest
DrAmES<-subset(data,BGroupBiome=="2.TropDryAmphibiansES")
median(DrAmES$RRS)
R<-10000
B97<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(DrAmES$RRS, replace = TRUE)
  B97[kk] <- mean(boot.sample)
}
boxplot(B97)
quantile(B97,c(0.025,0.975))
hist(B97, breaks = 30)
dB97<-c(Bmean=mean(B97),quantile(B97,c(0.025,0.975)),n=nrow(DrAmES))
dB97
#B98: Amphibians in Young Secondary Forest (YSF), Dry
DrAmYSF<-subset(data,BGroupBiome=="2.TropDryAmphibiansYSF")
median(DrAmYSF$RRS)
nrow(DrAmYSF)
R<-10000
B98<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(DrAmYSF$RRS, replace = TRUE)
  B98[kk] <- mean(boot.sample)
}
boxplot(B98)
quantile(B98,c(0.025,0.975))
hist(B98, breaks = 30)
dB98<-c(Bmean=mean(B98),quantile(B98,c(0.025,0.975)),n=nrow(DrAmYSF))
dB98
#B99: Amphibians in Mid-successional Secondary Forest (MSF), Dry
DrAmMSF<-subset(data,BGroupBiome=="2.TropDryAmphibiansMSF")
median(DrAmMSF$RRS)
R<-10000
B99<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(DrAmMSF$RRS, replace = TRUE)
  B99[kk] <- mean(boot.sample)
}
boxplot(B99)
quantile(B99,c(0.025,0.975))
hist(B99, breaks = 30)
dB99<-c(Bmean=mean(B99),quantile(B99,c(0.025,0.975)),n=nrow(DrAmMSF))
dB99
#B100: Amphibians in Old Secondary Forest (OSF), Dry (there is no data)
B100<-rep(NaN,10000)
dB100<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB100
#B101: Reptiles in ES, Dry
DrReES<-subset(data,BGroupBiome=="2.TropDryReptilesES")
median(DrReES$RRS)
R<-10000
B101<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(DrReES$RRS, replace = TRUE)
  B101[kk] <- mean(boot.sample)
}
boxplot(B101)
quantile(B101,c(0.025,0.975))
hist(B101, breaks = 30)
dB101<-c(Bmean=mean(B101),quantile(B101,c(0.025,0.975)),n=nrow(DrReES))
dB101
#B102: Reptiles in YSF, Dry
DrReYSF<-subset(data,BGroupBiome=="2.TropDryReptilesYSF")
median(DrReYSF$RRS)
R<-10000
B102<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(DrReYSF$RRS, replace = TRUE)
  B102[kk] <- mean(boot.sample)
}
boxplot(B102)
quantile(B102,c(0.025,0.975))
hist(B102, breaks = 30)
dB102<-c(Bmean=mean(B102),quantile(B102,c(0.025,0.975)),n=nrow(DrReYSF))
dB102
#B103: Reptiles in MSF, Dry
DrReMSF<-subset(data,BGroupBiome=="2.TropDryReptilesMSF")
median(DrReMSF$RRS)
R<-10000
B103<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(DrReMSF$RRS, replace = TRUE)
  B103[kk] <- mean(boot.sample)
}
boxplot(B103)
quantile(B103,c(0.025,0.975))
hist(B103, breaks = 30)
dB103<-c(Bmean=mean(B103),quantile(B103,c(0.025,0.975)),n=nrow(DrReMSF))
dB103
#B104: Reptiles in OSF, Dry
B104<-rep(NaN,10000)
dB104<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB104
#B105: Birds in ES, Dry
B105<-rep(NaN,10000)
dB105<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB105
#B26: Birds in YSF, Dry
DrBiYSF<-subset(data,BGroupBiome=="2.TropDryBirdsYSF")
median(DrBiYSF$RRS)
nrow(DrBiYSF)
R<-10000
B106<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(DrBiYSF$RRS, replace = TRUE)
  B106[kk] <- mean(boot.sample)
}
boxplot(B106)
quantile(B106,c(0.025,0.975))
hist(B106, breaks = 30)
dB106<-c(Bmean=mean(B106),quantile(B106,c(0.025,0.975)),n=nrow(DrBiYSF))
dB106
#B107: Birds in MSF, Dry
DrBiMSF<-subset(data,BGroupBiome=="2.TropDryBirdsMSF")
median(DrBiMSF$RRS)
nrow(DrBiMSF)
R<-10000
B107<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(DrBiMSF$RRS, replace = TRUE)
  B107[kk] <- mean(boot.sample)
}
boxplot(B107)
quantile(B107,c(0.025,0.975))
hist(B107, breaks = 30)
dB107<-c(Bmean=mean(B107),quantile(B107,c(0.025,0.975)),n=nrow(DrBiMSF))
dB107
#B108: Birds in OSF, Dry (no data)
B108<-rep(NaN,10000)
dB108<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB108
#B109: Mammals in ES, Dry
DrMaES<-subset(data,BGroupBiome=="2.TropDryMammalsES")
median(DrMaES$RRS)
nrow(DrMaES)
R<-10000
B109<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(DrMaES$RRS, replace = TRUE)
  B109[kk] <- mean(boot.sample)
}
boxplot(B109)
quantile(B109,c(0.025,0.975))
hist(B109, breaks = 30)
dB109<-c(Bmean=mean(B109),quantile(B109,c(0.025,0.975)),n=nrow(DrMaES))
dB109
#B110: Mammals in YSF, Dry
DrMaYSF<-subset(data,BGroupBiome=="2.TropDryMammalsYSF")
median(DrMaYSF$RRS)
nrow(DrMaYSF)
R<-10000
B110<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(DrMaYSF$RRS, replace = TRUE)
  B110[kk] <- mean(boot.sample)
}
boxplot(B110)
quantile(B110,c(0.025,0.975))
hist(B110, breaks = 30)
dB110<-c(Bmean=mean(B110),quantile(B110,c(0.025,0.975)),n=nrow(DrMaYSF))
dB110
#B111: Mammals in MSF, Dry
DrMaMSF<-subset(data,BGroupBiome=="2.TropDryMammalsMSF")
median(DrMaMSF$RRS)
nrow(DrMaMSF)
R<-10000
B111<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(DrMaMSF$RRS, replace = TRUE)
  B111[kk] <- mean(boot.sample)
}
boxplot(B111)
quantile(B111,c(0.025,0.975))
hist(B111, breaks = 30)
dB111<-c(Bmean=mean(B111),quantile(B111,c(0.025,0.975)),n=nrow(DrMaMSF))
dB111
#B112: Mammals in OSF, Dry
DrMaOSF<-subset(data,BGroupBiome=="2.TropDryMammalsOSF")
median(DrMaOSF$RRS)
nrow(DrMaOSF)
R<-10000
B112<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(DrMaOSF$RRS, replace = TRUE)
  B112[kk] <- mean(boot.sample)
}
boxplot(B112)
quantile(B112,c(0.025,0.975))
hist(B112, breaks = 30)
dB112<-c(Bmean=mean(B112),quantile(B112,c(0.025,0.975)),n=nrow(DrMaOSF))
dB112

##~~ Tropical Grasslands and Savannahs (ej. Orinoquia)
#B113: Amphibians in Early Succession (ES), Savannahs
SaAmES<-subset(data,BGroupBiome=="3.TropGrSSAmphibiansES")
median(SaAmES$RRS)
R<-10000
B113<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(SaAmES$RRS, replace = TRUE)
  B113[kk] <- mean(boot.sample)
}
boxplot(B113)
quantile(B113,c(0.025,0.975))
hist(B113, breaks = 30)
dB113<-c(Bmean=mean(B113),quantile(B113,c(0.025,0.975)),n=nrow(SaAmES))
dB113
#B114: Amphibians in Young Secondary Forest (YSF), Savannahs (there is no data)
B114<-rep(NaN,10000)
dB114<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB114
#B115: Amphibians in Mid-successional Secondary Forest (MSF), Savannahs
SaAmMSF<-subset(data,BGroupBiome=="3.TropGrSSAmphibiansMSF")
median(SaAmMSF$RRS)
R<-10000
B115<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(SaAmMSF$RRS, replace = TRUE)
  B115[kk] <- mean(boot.sample)
}
boxplot(B115)
quantile(B115,c(0.025,0.975))
hist(B115, breaks = 30)
dB115<-c(Bmean=mean(B115),quantile(B115,c(0.025,0.975)),n=nrow(SaAmMSF))
dB115
#B116: Amphibians in Old Secondary Forest (OSF), Savannahs (there is no data)
B116<-rep(NaN,10000)
dB116<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB116
#B117: Reptiles in ES, Savannahs
SaReES<-subset(data,BGroupBiome=="3.TropGrSSReptilesES")
median(SaReES$RRS)
R<-10000
B117<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(SaReES$RRS, replace = TRUE)
  B117[kk] <- mean(boot.sample)
}
boxplot(B117)
quantile(B117,c(0.025,0.975))
hist(B117, breaks = 30)
dB117<-c(Bmean=mean(B117),quantile(B117,c(0.025,0.975)),n=nrow(SaReES))
dB117
#B118: Reptiles in YSF, Savannahs
B118<-rep(NaN,10000)
dB118<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB118
#B119: Reptiles in MSF, Savannahs
SaReMSF<-subset(data,BGroupBiome=="3.TropGrSSReptilesMSF")
median(SaReMSF$RRS)
R<-10000
B119<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(SaReMSF$RRS, replace = TRUE)
  B119[kk] <- mean(boot.sample)
}
boxplot(B119)
quantile(B119,c(0.025,0.975))
hist(B119, breaks = 30)
dB119<-c(Bmean=mean(B119),quantile(B119,c(0.025,0.975)),n=nrow(SaReMSF))
dB119
#B120: Reptiles in OSF, Australasia (there is no data)
B120<-rep(NaN,10000)
dB120<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB120
#B121: Birds in ES, Savannahs
SaBiES<-subset(data,BGroupBiome=="3.TropGrSSBirdsES")
median(SaBiES$RRS)
nrow(SaBiES)
R<-10000
B121<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(SaBiES$RRS, replace = TRUE)
  B121[kk] <- mean(boot.sample)
}
boxplot(B121)
quantile(B121,c(0.025,0.975))
hist(B121, breaks = 30)
dB121<-c(Bmean=mean(B121),quantile(B121,c(0.025,0.975)),n=nrow(SaBiES))
dB121
#B122: Birds in YSF, Savannahs
B122<-rep(NaN,10000)
dB122<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB122
#B123: Birds in MSF, Savannahs
SaBiMSF<-subset(data,BGroupBiome=="3.TropGrSSBirdsMSF")
median(SaBiMSF$RRS)
nrow(SaBiMSF)
R<-10000
B123<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(SaBiMSF$RRS, replace = TRUE)
  B123[kk] <- mean(boot.sample)
}
boxplot(B123)
quantile(B123,c(0.025,0.975))
hist(B123, breaks = 30)
dB123<-c(Bmean=mean(B123),quantile(B123,c(0.025,0.975)),n=nrow(SaBiMSF))
dB123
#B124: Birds in OSF, Savannahs (no data)
B124<-rep(NaN,10000)
dB124<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB124
#B125: Mammals in ES, Savannahs
SaMaES<-subset(data,BGroupBiome=="3.TropGrSSMammalsES")
median(SaMaES$RRS)
nrow(SaMaES)
R<-10000
B125<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(SaMaES$RRS, replace = TRUE)
  B125[kk] <- mean(boot.sample)
}
boxplot(B125)
quantile(B125,c(0.025,0.975))
hist(B125, breaks = 30)
dB125<-c(Bmean=mean(B125),quantile(B125,c(0.025,0.975)),n=nrow(SaMaES))
dB125
#B126: Mammals in YSF, Savannahs (no data)
B126<-rep(NaN,10000)
dB126<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB126
#B127: Mammals in MSF, Savannahs
SaMaMSF<-subset(data,BGroupBiome=="3.TropGrSSMammalsMSF")
median(SaMaMSF$RRS)
nrow(SaMaMSF)
R<-10000
B127<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(SaMaMSF$RRS, replace = TRUE)
  B127[kk] <- mean(boot.sample)
}
boxplot(B127)
quantile(B127,c(0.025,0.975))
hist(B127, breaks = 30)
dB127<-c(Bmean=mean(B127),quantile(B127,c(0.025,0.975)),n=nrow(SaMaMSF))
dB127
#B128: Mammals in OSF, Savannahs (no data)
B128<-rep(NaN,10000)
dB128<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB128


##~~ Geographical Condition
#B129: Amphibians in Early Succession (ES), Continent
CoAmES<-subset(data,BGroupContIsla=="ContinentAmphibiansES")
median(CoAmES$RRS)
R<-10000
B129<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(CoAmES$RRS, replace = TRUE)
  B129[kk] <- mean(boot.sample)
}
boxplot(B129)
quantile(B129,c(0.025,0.975))
hist(B129, breaks = 30)
dB129<-c(Bmean=mean(B129),quantile(B129,c(0.025,0.975)),n=nrow(CoAmES))
dB129
#B130: Amphibians in Young Secondary Forest (YSF), Continent
CoAmYSF<-subset(data,BGroupContIsla=="ContinentAmphibiansYSF")
median(CoAmYSF$RRS)
R<-10000
B130<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(CoAmYSF$RRS, replace = TRUE)
  B130[kk] <- mean(boot.sample)
}
boxplot(B130)
quantile(B130,c(0.025,0.975))
hist(B130, breaks = 30)
dB130<-c(Bmean=mean(B130),quantile(B130,c(0.025,0.975)),n=nrow(CoAmYSF))
dB130
#B131: Amphibians in Mid-successional Secondary Forest (MSF), Continent
CoAmMSF<-subset(data,BGroupContIsla=="ContinentAmphibiansMSF")
median(CoAmMSF$RRS)
R<-10000
B131<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(CoAmMSF$RRS, replace = TRUE)
  B131[kk] <- mean(boot.sample)
}
boxplot(B131)
quantile(B131,c(0.025,0.975))
hist(B131, breaks = 30)
dB131<-c(Bmean=mean(B131),quantile(B131,c(0.025,0.975)),n=nrow(CoAmMSF))
dB131
#B132: Amphibians in Old Secondary Forest (OSF), Continent (no data)
B132<-rep(NaN,10000)
dB132<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB132
#B133: Reptiles in ES, Continent 
CoReES<-subset(data,BGroupContIsla=="ContinentReptilesES")
median(CoReES$RRS)
R<-10000
B133<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(CoReES$RRS, replace = TRUE)
  B133[kk] <- mean(boot.sample)
}
boxplot(B133)
quantile(B133,c(0.025,0.975))
hist(B133, breaks = 30)
dB133<-c(Bmean=mean(B133),quantile(B133,c(0.025,0.975)),n=nrow(CoReES))
dB133
#B134: Reptiles in YSF, Continent
CoReYSF<-subset(data,BGroupContIsla=="ContinentReptilesYSF")
median(CoReYSF$RRS)
R<-10000
B134<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(CoReYSF$RRS, replace = TRUE)
  B134[kk] <- mean(boot.sample)
}
boxplot(B134)
quantile(B134,c(0.025,0.975))
hist(B134, breaks = 30)
dB134<-c(Bmean=mean(B134),quantile(B134,c(0.025,0.975)),n=nrow(CoReYSF))
dB134
#B135: Reptiles in MSF, Continental
CoReMSF<-subset(data,BGroupContIsla=="ContinentReptilesMSF")
median(CoReMSF$RRS)
R<-10000
B135<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(CoReMSF$RRS, replace = TRUE)
  B135[kk] <- mean(boot.sample)
}
boxplot(B135)
quantile(B135,c(0.025,0.975))
hist(B135, breaks = 30)
dB135<-c(Bmean=mean(B135),quantile(B135,c(0.025,0.975)),n=nrow(CoReMSF))
dB135
#B136: Reptiles in OSF, Continent (there is no data)
B136<-rep(NaN,10000)
dB136<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB136
#B137: Birds in ES, Continent 
CoBiES<-subset(data,BGroupContIsla=="ContinentBirdsES")
median(CoBiES$RRS)
R<-10000
B137<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(CoBiES$RRS, replace = TRUE)
  B137[kk] <- mean(boot.sample)
}
boxplot(B137)
quantile(B137,c(0.025,0.975))
hist(B137, breaks = 30)
dB137<-c(Bmean=mean(B137),quantile(B137,c(0.025,0.975)),n=nrow(CoBiES))
dB137
#B138: Birds in YSF, Continent
CoBiYSF<-subset(data,BGroupContIsla=="ContinentBirdsYSF")
median(CoBiYSF$RRS)
nrow(CoBiYSF)
R<-10000
B138<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(CoBiYSF$RRS, replace = TRUE)
  B138[kk] <- mean(boot.sample)
}
boxplot(B138)
quantile(B138,c(0.025,0.975))
hist(B138, breaks = 30)
dB138<-c(Bmean=mean(B138),quantile(B138,c(0.025,0.975)),n=nrow(CoBiYSF))
dB138
#B139: Birds in MSF, Continent
CoBiMSF<-subset(data,BGroupContIsla=="ContinentBirdsMSF")
median(CoBiMSF$RRS)
nrow(CoBiMSF)
R<-10000
B139<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(CoBiMSF$RRS, replace = TRUE)
  B139[kk] <- mean(boot.sample)
}
boxplot(B139)
quantile(B139,c(0.025,0.975))
hist(B139, breaks = 30)
dB139<-c(Bmean=mean(B139),quantile(B139,c(0.025,0.975)),n=nrow(CoBiMSF))
dB139
#B140: Birds in OSF, Continent
B140<-rep(NaN,10000)
dB140<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB140
#B141: Mammals in ES,  Continent
CoMaOSF<-subset(data,BGroupContIsla=="ContinentMammalsOSF")
median(CoMaOSF$RRS)
nrow(CoMaOSF)
R<-10000
B141<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(CoMaOSF$RRS, replace = TRUE)
  B141[kk] <- mean(boot.sample)
}
boxplot(B141)
quantile(B141,c(0.025,0.975))
hist(B141, breaks = 30)
dB141<-c(Bmean=mean(B141),quantile(B141,c(0.025,0.975)),n=nrow(CoMaOSF))
dB141
#B142: Mammals in YSF, Continent
CoMaYSF<-subset(data,BGroupContIsla=="ContinentMammalsYSF")
median(CoMaYSF$RRS)
nrow(CoMaYSF)
R<-10000
B142<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(CoMaYSF$RRS, replace = TRUE)
  B142[kk] <- mean(boot.sample)
}
boxplot(B142)
quantile(B142,c(0.025,0.975))
hist(B142, breaks = 30)
dB142<-c(Bmean=mean(B142),quantile(B142,c(0.025,0.975)),n=nrow(CoMaYSF))
dB142
#B143: Mammals in MSF, Continent
CoMaMSF<-subset(data,BGroupContIsla=="ContinentMammalsMSF")
median(CoMaMSF$RRS)
nrow(CoMaMSF)
R<-10000
B143<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(CoMaMSF$RRS, replace = TRUE)
  B143[kk] <- mean(boot.sample)
}
boxplot(B143)
quantile(B143,c(0.025,0.975))
hist(B143, breaks = 30)
dB143<-c(Bmean=mean(B143),quantile(B143,c(0.025,0.975)),n=nrow(CoMaMSF))
dB143
#B144: Mammals in OSF, Continent
CoMaOSF<-subset(data,BGroupContIsla=="ContinentMammalsOSF")
median(CoMaOSF$RRS)
nrow(CoMaOSF)
R<-10000
B144<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(CoMaOSF$RRS, replace = TRUE)
  B144[kk] <- mean(boot.sample)
}
boxplot(B144)
quantile(B144,c(0.025,0.975))
hist(B144, breaks = 30)
dB144<-c(Bmean=mean(B144),quantile(B144,c(0.025,0.975)),n=nrow(CoMaOSF))
dB144

##~~ Islands
#B145: Amphibians in Early Succession (ES), Islands
IsAmES<-subset(data,BGroupContIsla=="IslandAmphibiansES")
median(IsAmES$RRS)
nrow(IsAmES)
R<-10000
B145<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(IsAmES$RRS, replace = TRUE)
  B145[kk] <- mean(boot.sample)
}
boxplot(B145)
quantile(B145,c(0.025,0.975))
hist(B145, breaks = 30)
dB145<-c(Bmean=mean(B145),quantile(B145,c(0.025,0.975)),n=nrow(IsAmES))
dB145
#B146: Amphibians in Young Secondary Forest (YSF), Islands
IsAmYSF<-subset(data,BGroupContIsla=="IslandAmphibiansYSF")
median(IsAmYSF$RRS)
R<-10000
B146<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(IsAmYSF$RRS, replace = TRUE)
  B146[kk] <- mean(boot.sample)
}
boxplot(B146)
quantile(B146,c(0.025,0.975))
hist(B146, breaks = 30)
dB146<-c(Bmean=mean(B146),quantile(B146,c(0.025,0.975)),n=nrow(IsAmYSF))
dB146
#B147: Amphibians in Mid-successional Secondary Forest (MSF), Island
IsAmMSF<-subset(data,BGroupContIsla=="IslandAmphibiansMSF")
median(IsAmMSF$RRS)
R<-10000
B147<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(IsAmMSF$RRS, replace = TRUE)
  B147[kk] <- mean(boot.sample)
}
boxplot(B147)
quantile(B147,c(0.025,0.975))
hist(B147, breaks = 30)
dB147<-c(Bmean=mean(B147),quantile(B147,c(0.025,0.975)),n=nrow(IsAmMSF))
dB147
#B148: Amphibians in Old Secondary Forest (OSF), Island
B148<-rep(NaN,10000)
dB148<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB148
#B149: Reptiles in ES, Island
IsReES<-subset(data,BGroupContIsla=="IslandReptilesES")
median(IsReES$RRS)
R<-10000
B149<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(IsReES$RRS, replace = TRUE)
  B149[kk] <- mean(boot.sample)
}
boxplot(B149)
quantile(B149,c(0.025,0.975))
hist(B149, breaks = 30)
dB149<-c(Bmean=mean(B149),quantile(B149,c(0.025,0.975)),n=nrow(IsReES))
dB149
#B150: Reptiles in YSF, Island
IsReYSF<-subset(data,BGroupContIsla=="IslandReptilesYSF")
median(IsReYSF$RRS)
R<-10000
B150<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(IsReYSF$RRS, replace = TRUE)
  B150[kk] <- mean(boot.sample)
}
boxplot(B150)
quantile(B150,c(0.025,0.975))
hist(B150, breaks = 30)
dB150<-c(Bmean=mean(B150),quantile(B150,c(0.025,0.975)),n=nrow(IsReYSF))
dB150
#B151: Reptiles in MSF, Island
IsReMSF<-subset(data,BGroupContIsla=="IslandReptilesMSF")
median(IsReMSF$RRS)
R<-10000
B151<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(IsReMSF$RRS, replace = TRUE)
  B151[kk] <- mean(boot.sample)
}
boxplot(B151)
quantile(B151,c(0.025,0.975))
hist(B151, breaks = 30)
dB151<-c(Bmean=mean(B151),quantile(B151,c(0.025,0.975)),n=nrow(IsReMSF))
dB151
#B152: Reptiles in OSF, Island
IsReOSF<-subset(data,BGroupContIsla=="IslandReptilesOSF")
median(IsReOSF$RRS)
R<-10000
B152<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(IsReOSF$RRS, replace = TRUE)
  B152[kk] <- mean(boot.sample)
}
boxplot(B152)
quantile(B152,c(0.025,0.975))
hist(B152, breaks = 30)
dB152<-c(Bmean=mean(B152),quantile(B152,c(0.025,0.975)),n=nrow(IsReOSF))
dB152
#B153: Birds in ES, Island
IsBiES<-subset(data,BGroupContIsla=="IslandBirdsES")
median(IsBiES$RRS)
nrow(IsBiES)
R<-10000
B153<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(IsBiES$RRS, replace = TRUE)
  B153[kk] <- mean(boot.sample)
}
boxplot(B153)
quantile(B153,c(0.025,0.975))
hist(B153, breaks = 30)
dB153<-c(Bmean=mean(B153),quantile(B153,c(0.025,0.975)),n=nrow(IsBiES))
dB153
#B154: Birds in YSF, Island
IsBiYSF<-subset(data,BGroupContIsla=="IslandBirdsYSF")
median(IsBiYSF$RRS)
nrow(IsBiYSF)
R<-10000
B154<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(IsBiYSF$RRS, replace = TRUE)
  B154[kk] <- mean(boot.sample)
}
boxplot(B154)
quantile(B154,c(0.025,0.975))
hist(B154, breaks = 30)
dB154<-c(Bmean=mean(B154),quantile(B154,c(0.025,0.975)),n=nrow(IsBiYSF))
dB154
#B155: Birds in MSF, Island)
IsBiMSF<-subset(data,BGroupContIsla=="IslandBirdsMSF")
median(IsBiMSF$RRS)
nrow(IsBiMSF)
R<-10000
B155<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(IsBiMSF$RRS, replace = TRUE)
  B155[kk] <- mean(boot.sample)
}
boxplot(B155)
quantile(B155,c(0.025,0.975))
hist(B155, breaks = 30)
dB155<-c(Bmean=mean(B155),quantile(B155,c(0.025,0.975)),n=nrow(IsBiMSF))
dB155
#B156: Birds in OSF, Islands
B156<-rep(NaN,10000)
dB156<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB156
#B157: Mammals in ES,  Islands
B157<-rep(NaN,10000)
dB157<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB157
#B158: Mammals in YSF, Island
IsMaYSF<-subset(data,BGroupContIsla=="IslandMammalsYSF")
median(IsMaYSF$RRS)
nrow(IsMaYSF)
R<-10000
B158<-rep(0,R);#change the name of the Bootstrap object
for(kk in 1:R){
  boot.sample <- sample(IsMaYSF$RRS, replace = TRUE)
  B158[kk] <- mean(boot.sample)
}
boxplot(B158)
quantile(B158,c(0.025,0.975))
hist(B158, breaks = 30)
dB158<-c(Bmean=mean(B158),quantile(B158,c(0.025,0.975)),n=nrow(IsMaYSF))
dB158
#B159: Mammals in MSF, Island (no data)
B159<-rep(NaN,10000)
dB159<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB159
#B160: Mammals in OSF, Island (no data)
B160<-rep(NaN,10000)
dB160<-c(Bmean=NaN,"2.5%" = NaN,"97.5%" = NaN,n=0)
dB160

#Join results
#Resume data of Bootstrap (Bootstrap mean, confidence limits[2.5%-97.5%],n)
dbR=data.frame(dB1,dB2,dB3,dB4,dB5,dB6,dB7,dB8,dB9,dB10, dB11,dB12,dB13,dB14,dB15,dB16,dB17,dB18,dB19,dB20, 
               dB21,dB22,dB23,dB24,dB25,dB26,dB27,dB28,dB29,dB30, dB31,dB32,dB33,dB34,dB35,dB36,dB37,dB38,dB39,dB40,
               dB41,dB42,dB43,dB44,dB45,dB46,dB47,dB48,dB49,dB50, dB51,dB52,dB53,dB54,dB55,dB56,dB57,dB58,dB59,dB60,
               dB61,dB62,dB63,dB64,dB65,dB66,dB67,dB68,dB69,dB70, dB71,dB72,dB73,dB74,dB75,dB76,dB77,dB78,dB79,dB80,
               dB81,dB82,dB83,dB84,dB85,dB86,dB87,dB88,dB89,dB90, dB91,dB92,dB93,dB94,dB95,dB96,dB97,dB98,dB99,dB100,
               dB101,dB102,dB103,dB104,dB105,dB106,dB107,dB108,dB109,dB110, dB111,dB112,dB113,dB114,dB115,dB116,dB117,dB118,dB119,dB120, 
               dB121,dB122,dB123,dB124,dB125,dB126,dB127,dB128,dB129,dB130, dB131,dB132,dB133,dB134,dB135,dB136,dB137,dB138,dB139, dB140,
               dB141,dB142,dB143,dB144,dB145,dB146,dB147,dB148,dB149,dB150, dB151,dB152,dB153,dB154,dB155,dB156,dB157,dB158,dB159,dB160)
t.dbR<-t(dbR)
datboots<-as.data.frame(t.dbR)
head(datboots)
datboots$FigurePart<-c(rep("Overall",16),rep("Neotropic",16),rep("Australasia",16),rep("IndoMalay",16),rep("Afrotropic",16),
                       rep("TropMo",16),rep("TropDry",16),rep("TropGrass",16),rep("Continent",16),rep("Island",16))
taxdb<-c(rep("Amphibians",4),rep("Reptiles",4),rep("Birds",4),rep("Mammals",4))
for.covdb<-c(rep("ES",1),rep("YSF",1),rep("MSF",1),rep("OSF",1))
datboots$Taxa<-c(rep(taxdb,10))
datboots$Forest.cov<-c(rep(for.covdb,4))
head(datboots)
tail(datboots)
write.table(datboots, "BootstrappSpSimil.txt",quote=F, sep="\t")

position_n<-datboots[,3]+0.3

#Data to graph the Bootstrap results
bR=data.frame(B1,B2,B3,B4,B5,B6,B7,B8,B9,B10, B11,B12,B13,B14,B15,B16,B17,B18,B19,B20, 
              B21,B22,B23,B24,B25,B26,B27,B28,B29,B30, B31,B32,B33,B34,B35,B36,B37,B38,B39,B40,
              B41,B42,B43,B44,B45,B46,B47,B48,B49,B50, B51,B52,B53,B54,B55,B56,B57,B58,B59,B60,
              B61,B62,B63,B64,B65,B66,B67,B68,B69,B70, B71,B72,B73,B74,B75,B76,B77,B78,B79,B80,
              B81,B82,B83,B84,B85,B86,B87,B88,B89,B90, B91,B92,B93,B94,B95,B96,B97,B98,B99,B100,
              B101,B102,B103,B104,B105,B106,B107,B108,B109, B110, B111,B112,B113,B114,B115,B116,B117,B118,B119,B120, 
              B121,B122,B123,B124,B125,B126,B127,B128, B129, B130, B131,B132,B133,B134,B135,B136,B137,B138, B139, B140,
              B141,B142,B143,B144,B145,B146,B147,B148, B149, B150, B151,B152,B153,B154,B155,B156,B157,B158, B159, B160)
head(bR)
bR$IDB<-seq(1:10000)
tail(bR)
boots<-melt(as.data.frame(bR), id=c("IDB"))
summary(boots)
boots$FigurePart<-c(rep("Overall",160000),rep("Neotropic",160000),rep("Australasia",160000),rep("IndoMalay",160000),rep("Afrotropic",160000),
                    rep("TropMo",160000),rep("TropDry",160000),rep("TropGrass",160000),rep("Continent",160000),rep("Island",160000))
tax<-c(rep("Amphibians",40000),rep("Reptiles",40000),rep("Birds",40000),rep("Mammals",40000))
for.cov<-c(rep("ES",10000),rep("YSF",10000),rep("MSF",10000),rep("OSF",10000))
boots$Taxa<-c(rep(tax,10))
boots$Forest.cov<-c(rep(for.cov,40))
head(boots)
tail(boots)
colnames(boots)<-c("IDB","Code","ResponseR","FigurePart","Taxa","Forest.cov")
overall<-subset(boots,FigurePart=="Overall")
tail(overall)
neot<-subset(boots,FigurePart=="Neotropic")
tail(neot)
austral<-subset(boots,FigurePart=="Australasia")
tail(austral)
indomalay<-subset(boots,FigurePart=="IndoMalay")
tail(indomalay)
afro<-subset(boots,FigurePart=="Afrotropic")
tail(afro)
moist<-subset(boots,FigurePart=="TropMo")
tail(moist)
dry<-subset(boots,FigurePart=="TropDry")
tail(dry)
grass<-subset(boots,FigurePart=="TropGrass")
tail(grass)
continent<-subset(boots,FigurePart=="Continent")
tail(continent)
island<-subset(boots,FigurePart=="Island")
tail(island)
#To graph the mean and quartiles (CI)
mean_quantile<-function(x){
  out<-quantile(x, probs = c(0.025, 0.975))
  names(out) <- c("ymin", "ymax")
  return(out)
}
#Lets graph all this work
overall$Taxa <- factor(overall$Taxa,levels=c('Amphibians','Reptiles','Birds', 'Mammals'))
overall$Forest.cov <- factor(overall$Forest.cov,levels=c('ES','YSF','MSF', 'OSF'))
#theme(legend.position="none")+
a=ggplot(overall, aes(x=Forest.cov, y=ResponseR, shape=Forest.cov)) + 
  #geom_boxplot(notch = TRUE, outlier.colour = "white", coef=0)+
  geom_hline(yintercept = 0,linetype="solid", color="gray", size=0.5)+
  stat_summary(fun.y = "mean", geom="point", 
               position = position_dodge(width = 1), fill="black", size=2)+
  stat_summary(fun.data = mean_quantile, geom="errorbar", width=0.2, size=0.5)+
  ggtitle('a) Overall')+
  ylab("")+  
  xlab("")+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = overall[1, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = overall[90000, ])+
  scale_y_continuous(breaks = seq(-2,0.4, by=0.5),
                     labels = seq(-2,0.4, by=0.5),
                     limits = c(-2,0.4), expand = c(0, 0))+
  scale_shape_manual(name="Forest.cov",
                     values = c("ES" = 21, "YSF"=22, 
                                "MSF"=23,"OSF"=24))+
  facet_wrap(~Taxa,ncol=4,nrow=1)+
  theme_bw()+
  theme(legend.position="none",axis.text.x=element_blank(), axis.text.y=element_text(size=11), 
        plot.margin=unit(c(1,1,-2,1),"mm"),panel.margin = unit(0, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(strip.background = element_blank(),plot.title=element_text(hjust=-0.1, size = 13, face = "bold"))+
  theme(strip.text = element_text(size=11))

Fig3a<-a+
  annotate("text",x=1,y=-0.07,label=c("8","","",""), size=3.5)+
  annotate("text",x=2,y=-0.09,label=c("20","","",""), size=3.5)+
  annotate("text",x=3,y=0.02,label=c("19","","",""), size=3.5)+
  annotate("text",x=1,y=-0.04,label=c("","6","",""), size=3.5)+
  annotate("text",x=2,y=-0.14,label=c("","19","",""), size=3.5)+
  annotate("text",x=3,y=0.06,label=c("","11","",""), size=3.5)+
  annotate("text",x=4,y=0.23,label=c("","1","",""), size=3.5)+
  annotate("text",x=1,y=-0.79,label=c("","","5",""), size=3.5)+
  annotate("text",x=2,y=-0.05,label=c("","","25",""), size=3.5)+
  annotate("text",x=3,y=0.18,label=c("","","13",""), size=3.5)+
  annotate("text",x=1,y=-0.11,label=c("","","","6"), size=3.5)+
  annotate("text",x=2,y=0.05,label=c("","","","30"), size=3.5)+
  annotate("text",x=3,y=0.05,label=c("","","","5"), size=3.5)+
  annotate("text",x=4,y=0.08,label=c("","","","2"), size=3.5)  
Fig3a #Figure of the bootstrapped effect size of species composition similarity Overall
#Now the Neotropic biogeographic realm
neot$Taxa <- factor(neot$Taxa,levels=c('Amphibians','Reptiles','Birds', 'Mammals'))
neot$Forest.cov <- factor(neot$Forest.cov,levels=c('ES','YSF','MSF', 'OSF'))
#theme(legend.position="none")+
b=ggplot(neot, aes(x=Forest.cov, y=ResponseR, shape=Forest.cov)) + 
  #geom_boxplot(notch = TRUE, outlier.colour = "white", coef=0)+
  geom_hline(yintercept = 0,linetype="solid", color="gray", size=0.5)+
  stat_summary(fun.y = "mean", geom="point", 
               position = position_dodge(width = 1), fill="black", size=2)+
  stat_summary(fun.data = mean_quantile, geom="errorbar", width=0.2, size=0.5)+
  ggtitle('b) Neotropic')+
  ylab("")+  
  xlab("")+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = neot[1, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = neot[90000, ])+
  scale_y_continuous(breaks = seq(-2,0.4, by=0.5),
                     labels = seq(-2,0.4, by=0.5),
                     limits = c(-2,0.4), expand = c(0, 0))+
  scale_shape_manual(name="Forest.cov",
                     values = c("ES" = 21, "YSF"=22, 
                                "MSF"=23,"OSF"=24))+
  facet_wrap(~Taxa,ncol=4,nrow=1)+
  theme_bw()+
  theme(legend.position="none",axis.text.x=element_blank(), axis.text.y=element_text(size=11),  
        plot.margin=unit(c(2,1,1,1),"mm"),panel.margin = unit(0, "lines"), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+ 
  theme(strip.background = element_blank(), plot.title=element_text(hjust=-0.1, size = 13, face = "bold"))+
  theme(strip.text = element_blank())
b
Fig3b<-b+
  annotate("text",x=1,y=-0.16,label=c("7","","",""), size=3.5)+
  annotate("text",x=2,y=0.02,label=c("16","","",""), size=3.5)+
  annotate("text",x=3,y=0.03,label=c("15","","",""), size=3.5)+
  annotate("text",x=1,y=-0.06,label=c("","3","",""), size=3.5)+
  annotate("text",x=2,y=-0.01,label=c("","12","",""), size=3.5)+
  annotate("text",x=3,y=0.04,label=c("","8","",""), size=3.5)+
  annotate("text",x=4,y=0.23,label=c("","1","",""), size=3.5)+
  annotate("text",x=1,y=-0.97,label=c("","","1",""), size=3.5)+
  annotate("text",x=2,y=0.03,label=c("","","17",""), size=3.5)+
  annotate("text",x=3,y=0.26,label=c("","","5",""), size=3.5)+
  annotate("text",x=1,y=-0.07,label=c("","","","5"), size=3.5)+
  annotate("text",x=2,y=0.19,label=c("","","","17"), size=3.5)+
  annotate("text",x=3,y=0.08,label=c("","","","3"), size=3.5)+
  annotate("text",x=4,y=0.08,label=c("","","","2"), size=3.5)  
Fig3b #Figure of the bootstrapped effect size of species composition similarity Neotropical realm
#Now the Australasia Biogeographic realm
austral$Taxa <- factor(austral$Taxa,levels=c('Amphibians','Reptiles','Birds', 'Mammals'))
austral$Forest.cov <- factor(austral$Forest.cov,levels=c('ES','YSF','MSF', 'OSF'))
#theme(legend.position="none")+
c=ggplot(austral, aes(x=Forest.cov, y=ResponseR, shape=Forest.cov)) + 
  #geom_boxplot(notch = TRUE, outlier.colour = "white", coef=0)+
  geom_hline(yintercept = 0,linetype="solid", color="gray", size=0.5)+
  stat_summary(fun.y = "mean", geom="point", 
               position = position_dodge(width = 1), fill="black", size=2)+
  stat_summary(fun.data = mean_quantile, geom="errorbar", width=0.2, size=0.5)+
  ggtitle('c) Australasia')+
  ylab("")+  
  xlab("")+
  geom_rect(aes(xmin=1.6, xmax=2.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = austral[1, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = austral[1, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = austral[80000, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = austral[90001, ])+
  geom_rect(aes(xmin=1.6, xmax=2.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = austral[150000, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = austral[150000, ])+
  scale_y_continuous(breaks = seq(-2,0.4, by=0.5),
                     labels = seq(-2,0.4, by=0.5),
                     limits = c(-2,0.4), expand = c(0, 0))+
  scale_x_discrete(labels = c("ES \n", 
                              "YSF \n",
                              "MSF \n",
                              "OSF \n"))+
  scale_shape_manual(name="Forest.cov",
                     values = c("ES" = 21, "YSF"=22, 
                                "MSF"=23,"OSF"=24))+
  facet_wrap(~Taxa,ncol=4,nrow=1)+
  theme_bw()+
  theme(legend.position="none", axis.text.x=element_blank(),axis.text.y=element_text(size=11),  
        plot.margin=unit(c(2,1,1,1),"mm"),panel.margin = unit(0, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(strip.background = element_blank(), plot.title=element_text(hjust=-0.1, size = 13, face = "bold"))+
  theme(strip.text = element_blank())
c
Fig3c<-c+
  annotate("text",x=1,y=0.15,label=c("1","","",""), size=3.5)+
  annotate("text",x=3,y=-0.02,label=c("1","","",""), size=3.5)+
  annotate("text",x=1,y=-0.02,label=c("","3","",""), size=3.5)+
  annotate("text",x=2,y=0.17,label=c("","3","",""), size=3.5)+
  annotate("text",x=3,y=0.13,label=c("","2","",""), size=3.5)+
  annotate("text",x=1,y=-0.54,label=c("","","3",""), size=3.5)+
  annotate("text",x=2,y=-0.17,label=c("","","2",""), size=3.5)+
  annotate("text",x=3,y=0.19,label=c("","","5",""), size=3.5)+
  annotate("text",x=1,y=-0.33,label=c("","","","1"), size=3.5)+
  annotate("text",x=3,y=0.11,label=c("","","","1"), size=3.5)  
Fig3c 
#Now the Indo Malay realm
indomalay$Taxa <- factor(indomalay$Taxa,levels=c('Amphibians','Reptiles','Birds', 'Mammals'))
indomalay$Forest.cov <- factor(indomalay$Forest.cov,levels=c('ES','YSF','MSF', 'OSF'))
#theme(legend.position="none")+
d=ggplot(indomalay, aes(x=Forest.cov, y=ResponseR, shape=Forest.cov)) + 
  #geom_boxplot(notch = TRUE, outlier.colour = "white", coef=0)+
  geom_hline(yintercept = 0,linetype="solid", color="gray", size=0.5)+
  stat_summary(fun.y = "mean", geom="point", 
               position = position_dodge(width = 1), fill="black", size=2)+
  stat_summary(fun.data = mean_quantile, geom="errorbar", width=0.2, size=0.5)+
  ggtitle('d) Indo Malay')+
  ylab("")+  
  xlab("")+
  geom_rect(aes(xmin=0.6, xmax=1.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = indomalay[1, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = indomalay[1, ])+
  geom_rect(aes(xmin=0.6, xmax=1.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = indomalay[60000, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = indomalay[60000, ])+
  geom_rect(aes(xmin=0.6, xmax=1.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = indomalay[90000, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = indomalay[90000, ])+
  geom_rect(aes(xmin=0.6, xmax=1.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = indomalay[130000, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = indomalay[130000, ])+
  scale_y_continuous(breaks = seq(-2,0.4, by=0.5),
                     labels = seq(-2,0.4, by=0.5),
                     limits = c(-2,0.4), expand = c(0, 0))+
  scale_shape_manual(name="Forest.cov",
                     values = c("ES" = 21, "YSF"=22, 
                                "MSF"=23,"OSF"=24))+
  facet_wrap(~Taxa,ncol=4,nrow=1)+
  theme_bw()+
  theme(legend.position="none",axis.text.x=element_blank(), axis.text.y=element_text(size=11),  
        plot.margin=unit(c(2,1,1,1),"mm"),panel.margin = unit(0, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(strip.background = element_blank(), plot.title=element_text(hjust=-.1, size = 13, face = "bold"))+
  theme(strip.text = element_blank())
d
Fig3d<-d+
  annotate("text",x=2,y=-1.33,label=c("3","","",""), size=3.5)+
  annotate("text",x=3,y=-0.15,label=c("1","","",""), size=3.5)+
  annotate("text",x=2,y=-0.53,label=c("","4","",""), size=3.5)+
  annotate("text",x=3,y=0.16,label=c("","1","",""), size=3.5)+
  annotate("text",x=2,y=-0.09,label=c("","","4",""), size=3.5)+
  annotate("text",x=3,y=0.16,label=c("","","3",""), size=3.5)+
  annotate("text",x=2,y=-0.20,label=c("","","","10"), size=3.5)+
  annotate("text",x=3,y=0.10,label=c("","","","1"), size=3.5)  
Fig3d #Figure of the bootstrapped effect size of richeness Tropical Moist forest
#Now the Afrotropic
afro$Taxa <- factor(afro$Taxa,levels=c('Amphibians','Reptiles','Birds', 'Mammals'))
afro$Forest.cov <- factor(afro$Forest.cov,levels=c('ES','YSF','MSF', 'OSF'))
#theme(legend.position="none")+
e=ggplot(afro, aes(x=Forest.cov, y=ResponseR, shape=Forest.cov)) + 
  #geom_boxplot(notch = TRUE, outlier.colour = "white", coef=0)+
  geom_hline(yintercept = 0,linetype="solid", color="gray", size=0.5)+
  stat_summary(fun.y = "mean", geom="point", 
               position = position_dodge(width = 1), fill="black", size=2)+
  stat_summary(fun.data = mean_quantile, geom="errorbar", width=0.2, size=0.5)+
  ggtitle('e) Afrotropic')+
  ylab("")+  
  xlab("")+
  geom_rect(aes(xmin=0.6, xmax=1.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = afro[20000, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = afro[40000, ])+
  geom_rect(aes(xmin=0.6, xmax=1.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = afro[60000, ])+
  geom_rect(aes(xmin=1.6, xmax=2.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = afro[60000, ])+
  geom_rect(aes(xmin=2.6, xmax=3.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = afro[60000, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = afro[60000, ])+
  geom_rect(aes(xmin=2.6, xmax=3.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = afro[100000, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = afro[100000, ])+
  geom_rect(aes(xmin=0.6, xmax=1.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = afro[160000, ])+
  geom_rect(aes(xmin=2.6, xmax=3.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = afro[160000, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = afro[160000, ])+
  scale_y_continuous(breaks = seq(-2,0.4, by=0.5),
                     labels = seq(-2,0.4, by=0.5),
                     limits = c(-2,0.4), expand = c(0, 0))+
  scale_x_discrete(labels = c("ES", 
                              "YSF",
                              "MSF",
                              "OSF"))+
  scale_shape_manual(name="Forest.cov",
                     values = c("ES" = 21, "YSF"=22, 
                                "MSF"=23,"OSF"=24))+
  facet_wrap(~Taxa,ncol=4,nrow=1)+
  theme_bw()+
  theme(legend.position="none",axis.text.y=element_text(size=11),axis.text.x=element_text(size=11, angle = 60, hjust = 1),  
        plot.margin=unit(c(2,1,-6,1),"mm"),panel.margin = unit(0, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(strip.background = element_blank(), plot.title=element_text(hjust=-.1, size = 13, face = "bold"))+
  theme(strip.text = element_blank())
e
Fig3e<-e+
  annotate("text",x=2,y=-0.43,label=c("1","","",""), size=3.5)+
  annotate("text",x=3,y=0.06,label=c("2","","",""), size=3.5)+
  annotate("text",x=1,y=-1.60,label=c("","","1",""), size=3.5)+
  annotate("text",x=2,y=0.01,label=c("","","2",""), size=3.5)+
   annotate("text",x=2,y=0.25,label=c("","","","3"), size=3.5)
Fig3e #Figure of the bootstrapped effect size of richeness Tropical Dry forest
#Now the Tropical Moist Forest Biome
moist$Taxa <- factor(moist$Taxa,levels=c('Amphibians','Reptiles','Birds', 'Mammals'))
moist$Forest.cov <- factor(moist$Forest.cov,levels=c('ES','YSF','MSF', 'OSF'))
#theme(legend.position="none")+
f=ggplot(moist, aes(x=Forest.cov, y=ResponseR, shape=Forest.cov)) + 
  #geom_boxplot(notch = TRUE, outlier.colour = "white", coef=0)+
  geom_hline(yintercept = 0,linetype="solid", color="gray", size=0.5)+
  stat_summary(fun.y = "mean", geom="point", 
               position = position_dodge(width = 1), fill="black", size=2)+
  stat_summary(fun.data = mean_quantile, geom="errorbar", width=0.2, size=0.5)+
  ggtitle('f) Tropical moist forest')+
  ylab("")+  
  xlab("")+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = moist[1, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = moist[100000, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = moist[150000, ])+
  scale_y_continuous(breaks = seq(-2,0.4, by=0.5),
                     labels = seq(-2,0.4, by=0.5),
                     limits = c(-2,0.4), expand = c(0, 0))+
  scale_shape_manual(name="Forest.cov",
                     values = c("ES" = 21, "YSF"=22, 
                                "MSF"=23,"OSF"=24))+
  facet_wrap(~Taxa,ncol=4,nrow=1)+
  theme_bw()+
  theme(legend.position="none",axis.text.x=element_blank(), axis.text.y=element_blank(),  
        plot.margin=unit(c(1,1,-2,1),"mm"),panel.margin = unit(0, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(strip.background = element_blank(), plot.title=element_text(hjust=-.1, size = 13, face = "bold"))+
  theme(strip.text = element_text(size=11))
f
Fig3f<-f+
  annotate("text",x=1,y=-0.12,label=c("6","","",""), size=3.5)+
  annotate("text",x=2,y=-0.10,label=c("18","","",""), size=3.5)+
  annotate("text",x=3,y=0.03,label=c("17","","",""), size=3.5)+
  annotate("text",x=1,y=-0.02,label=c("","4","",""), size=3.5)+
  annotate("text",x=2,y=-0.13,label=c("","17","",""), size=3.5)+
  annotate("text",x=3,y=0.08,label=c("","9","",""), size=3.5)+
  annotate("text",x=4,y=0.23,label=c("","1","",""), size=3.5)+
  annotate("text",x=1,y=-1.02,label=c("","","4",""), size=3.5)+
  annotate("text",x=2,y=-0.11,label=c("","","21",""), size=3.5)+
  annotate("text",x=3,y=0.17,label=c("","","8",""), size=3.5)+
  annotate("text",x=1,y=-0.13,label=c("","","","1"), size=3.5)+
  annotate("text",x=2,y=-0.02,label=c("","","","24"), size=3.5)+
  annotate("text",x=3,y=0.1,label=c("","","","2"), size=3.5)  
Fig3f #Figure of the bootstrapped effect size of richeness Tropical Moist forest
#Now the Tropical Dry Forest Biome
dry$Taxa <- factor(dry$Taxa,levels=c('Amphibians','Reptiles','Birds', 'Mammals'))
dry$Forest.cov <- factor(dry$Forest.cov,levels=c('ES','YSF','MSF', 'OSF'))
#theme(legend.position="none")+
g=ggplot(dry, aes(x=Forest.cov, y=ResponseR, shape=Forest.cov)) + 
  #geom_boxplot(notch = TRUE, outlier.colour = "white", coef=0)+
  geom_hline(yintercept = 0,linetype="solid", color="gray", size=0.5)+
  stat_summary(fun.y = "mean", geom="point", 
               position = position_dodge(width = 1), fill="black", size=2)+
  stat_summary(fun.data = mean_quantile, geom="errorbar", width=0.2, size=0.5)+
  ggtitle('g) Tropical dry forest')+
  ylab("")+  
  xlab("")+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = dry[1, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = dry[80000, ])+
  geom_rect(aes(xmin=0.6, xmax=1.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = dry[100000, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = dry[100000, ])+
  scale_y_continuous(breaks = seq(-2,0.4, by=0.5),
                     labels = seq(-2,0.4, by=0.5),
                     limits = c(-2,0.4), expand = c(0, 0))+
  scale_x_discrete(labels = c("ES \n", 
                              "YSF \n",
                              "MSF \n",
                              "OSF \n"))+
  scale_shape_manual(name="Forest.cov",
                     values = c("ES" = 21, "YSF"=22, 
                                "MSF"=23,"OSF"=24))+
  facet_wrap(~Taxa,ncol=4,nrow=1)+
  theme_bw()+
  theme(legend.position="none", axis.text.x=element_blank(),axis.text.y=element_blank(),  
        plot.margin=unit(c(2,1,1,1),"mm"),panel.margin = unit(0, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(strip.background = element_blank(), plot.title=element_text(hjust=0.01, size = 13, face = "bold"))+
  theme(strip.text = element_blank())
g
Fig3g<-g+
  annotate("text",x=1,y=-0.51,label=c("1","","",""), size=3.5)+
  annotate("text",x=2,y=-0.05,label=c("2","","",""), size=3.5)+
  annotate("text",x=3,y=-0.17,label=c("1","","",""), size=3.5)+
  annotate("text",x=1,y=-0.06,label=c("","1","",""), size=3.5)+
  annotate("text",x=2,y=-0.11,label=c("","2","",""), size=3.5)+
  annotate("text",x=3,y=0.03,label=c("","1","",""), size=3.5)+
  annotate("text",x=2,y=0.24,label=c("","","4",""), size=3.5)+
  annotate("text",x=3,y=0.26,label=c("","","4",""), size=3.5)+
  annotate("text",x=1,y=-0.30,label=c("","","","3"), size=3.5)+
  annotate("text",x=2,y=0.27,label=c("","","","6"), size=3.5)+
  annotate("text",x=3,y=0.08,label=c("","","","2"), size=3.5)+
  annotate("text",x=4,y=0.08,label=c("","","","2"), size=3.5)  
Fig3g #Figure of the bootstrapped effect size of richeness Tropical Dry forest
#Now the Tropical Savannas Biome
grass$Taxa <- factor(dry$Taxa,levels=c('Amphibians','Reptiles','Birds', 'Mammals'))
grass$Forest.cov <- factor(dry$Forest.cov,levels=c('ES','YSF','MSF', 'OSF'))
#theme(legend.position="none")+
h=ggplot(grass, aes(x=Forest.cov, y=ResponseR, shape=Forest.cov)) + 
  #geom_boxplot(notch = TRUE, outlier.colour = "white", coef=0)+
  geom_hline(yintercept = 0,linetype="solid", color="gray", size=0.5)+
  stat_summary(fun.y = "mean", geom="point", 
               position = position_dodge(width = 1), fill="black", size=2)+
  stat_summary(fun.data = mean_quantile, geom="errorbar", width=0.2, size=0.5)+
  ggtitle('h) Tropical savannahs')+
  ylab("")+  
  xlab("")+
  geom_rect(aes(xmin=1.6, xmax=2.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = grass[20000, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = grass[40000, ])+
  geom_rect(aes(xmin=1.6, xmax=2.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = grass[60000, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = grass[80000, ])+
  geom_rect(aes(xmin=1.6, xmax=2.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = grass[100000, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = grass[120000, ])+
  geom_rect(aes(xmin=1.6, xmax=2.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = grass[140000, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = grass[160000, ])+
  scale_y_continuous(breaks = seq(-2,0.4, by=0.5),
                     labels = seq(-2,0.4, by=0.5),
                     limits = c(-2,0.4), expand = c(0, 0))+
  scale_x_discrete(labels = c("ES \n", 
                              "YSF \n",
                              "MSF \n",
                              "OSF \n"))+
  scale_shape_manual(name="Forest.cov",
                     values = c("ES" = 21, "YSF"=22, 
                                "MSF"=23,"OSF"=24))+
  facet_wrap(~Taxa,ncol=4,nrow=1)+
  theme_bw()+
  theme(legend.position="none",axis.text.y=element_blank(),axis.text.x=element_blank(),  
        plot.margin=unit(c(2,1,1,1),"mm"),panel.margin = unit(0, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(strip.background = element_blank(), plot.title=element_text(hjust=0.01, size = 13, face = "bold"))+
  theme(strip.text = element_blank())
h
Fig3h<-h+
  annotate("text",x=1,y=0.15,label=c("1","","",""), size=3.5)+
  annotate("text",x=3,y=-0.02,label=c("1","","",""), size=3.5)+
  annotate("text",x=1,y=-0.31,label=c("","1","",""), size=3.5)+
  annotate("text",x=3,y=-0.01,label=c("","1","",""), size=3.5)+
  annotate("text",x=1,y=-0.54,label=c("","","1",""), size=3.5)+
  annotate("text",x=3,y=-0.01,label=c("","","1",""), size=3.5)+
  annotate("text",x=1,y=0.03,label=c("","","","2"), size=3.5)+
  annotate("text",x=3,y=-0.10,label=c("","","","1"), size=3.5)
Fig3h #Figure of the bootstrapped effect size of richeness Tropical Dry forest
#Now the Continent
continent$Taxa <- factor(continent$Taxa,levels=c('Amphibians','Reptiles','Birds', 'Mammals'))
continent$Forest.cov <- factor(continent$Forest.cov,levels=c('ES','YSF','MSF', 'OSF'))
#theme(legend.position="none")+
i=ggplot(continent, aes(x=Forest.cov, y=ResponseR, shape=Forest.cov)) + 
  #geom_boxplot(notch = TRUE, outlier.colour = "white", coef=0)+
  geom_hline(yintercept = 0,linetype="solid", color="gray", size=0.5)+
  stat_summary(fun.y = "mean", geom="point", 
               position = position_dodge(width = 1), fill="black", size=2)+
  stat_summary(fun.data = mean_quantile, geom="errorbar", width=0.2, size=0.5)+
  ggtitle('i) Continent')+
  ylab("")+  
  xlab("")+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = afro[40000, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = afro[60000, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = afro[100000, ])+
  scale_y_continuous(breaks = seq(-2,0.4, by=0.5),
                     labels = seq(-2,0.4, by=0.5),
                     limits = c(-2,0.4), expand = c(0, 0))+
  scale_x_discrete(labels = c("ES \n", 
                              "YSF \n",
                              "MSF \n",
                              "OSF \n"))+
  scale_shape_manual(name="Forest.cov",
                     values = c("ES" = 21, "YSF"=22, 
                                "MSF"=23,"OSF"=24))+
  facet_wrap(~Taxa,ncol=4,nrow=1)+
  theme_bw()+
  theme(legend.position="none",axis.text.y=element_blank(),axis.text.x=element_blank(),  
        plot.margin=unit(c(2,1,1,1),"mm"),panel.margin = unit(0, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(strip.background = element_blank(), plot.title=element_text(hjust=0.01, size = 13, face = "bold"))+
  theme(strip.text = element_blank())
i
Fig3i<-i+
  annotate("text",x=1,y=-0.07,label=c("6","","",""), size=3.5)+
  annotate("text",x=2,y=-0.1,label=c("18","","",""), size=3.5)+
  annotate("text",x=3,y=0.01,label=c("16","","",""), size=3.5)+
  annotate("text",x=1,y=-0.06,label=c("","2","",""), size=3.5)+
  annotate("text",x=2,y=-0.35,label=c("","13","",""), size=3.5)+
  annotate("text",x=3,y=-0.01,label=c("","8","",""), size=3.5)+
  annotate("text",x=1,y=-0.54,label=c("","","3",""), size=3.5)+
  annotate("text",x=2,y=-0.02,label=c("","","22",""), size=3.5)+
  annotate("text",x=3,y=0.22,label=c("","","7",""), size=3.5)+
  annotate("text",x=1,y=0.08,label=c("","","","2"), size=3.5)+
  annotate("text",x=2,y=0.07,label=c("","","","29"), size=3.5)+
  annotate("text",x=3,y=0.05,label=c("","","","5"), size=3.5)+
  annotate("text",x=4,y=0.08,label=c("","","","2"), size=3.5)  
Fig3i #Figure of the bootstrapped effect size of richeness in Continent
#Now the Islands
island$Taxa <- factor(island$Taxa,levels=c('Amphibians','Reptiles','Birds', 'Mammals'))
island$Forest.cov <- factor(island$Forest.cov,levels=c('ES','YSF','MSF', 'OSF'))
#theme(legend.position="none")+
j=ggplot(island, aes(x=Forest.cov, y=ResponseR, shape=Forest.cov)) + 
  #geom_boxplot(notch = TRUE, outlier.colour = "white", coef=0)+
  geom_hline(yintercept = 0,linetype="solid", color="gray", size=0.5)+
  stat_summary(fun.y = "mean", geom="point", 
               position = position_dodge(width = 1), fill="black", size=2)+
  stat_summary(fun.data = mean_quantile, geom="errorbar", width=0.2, size=0.5)+
  ggtitle('j) Island')+
  ylab("")+  
  xlab("")+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = island[40000, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = island[120000, ])+
  geom_rect(aes(xmin=0.6, xmax=1.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = island[160000, ])+
  geom_rect(aes(xmin=2.6, xmax=3.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = island[160000, ])+
  geom_rect(aes(xmin=3.6, xmax=4.4, ymin=-Inf, ymax=Inf),alpha=0.5,fill='gray90', data = island[160000, ])+
  scale_y_continuous(breaks = seq(-2,0.4, by=0.5),
                     labels = seq(-2,0.4, by=0.5),
                     limits = c(-2,0.4), expand = c(0, 0))+
  scale_x_discrete(labels = c("ES", 
                              "YSF",
                              "MSF",
                              "OSF"))+
  scale_shape_manual(name="Forest.cov",
                     values = c("ES" = 21, "YSF"=22, 
                                "MSF"=23,"OSF"=24))+
  facet_wrap(~Taxa,ncol=4,nrow=1)+
  theme_bw()+
  theme(legend.position="none",axis.text.y=element_blank(),axis.text.x=element_text(size=11, angle = 60, hjust = 1),  
        plot.margin=unit(c(2,1,-6,1),"mm"),panel.margin = unit(0, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  theme(strip.background = element_blank(), plot.title=element_text(hjust=0.01, size = 13, face = "bold"))+
  theme(strip.text = element_blank())
j
Fig3j<-j+
  annotate("text",x=1,y=0.01,label=c("2","","",""), size=3.5)+
  annotate("text",x=2,y=0.3,label=c("2","","",""), size=3.5)+
  annotate("text",x=3,y=0.18,label=c("3","","",""), size=3.5)+
  annotate("text",x=1,y=-0.02,label=c("","4","",""), size=3.5)+
  annotate("text",x=2,y=0.19,label=c("","6","",""), size=3.5)+
  annotate("text",x=3,y=0.20,label=c("","3","",""), size=3.5)+
  annotate("text",x=4,y=0.23,label=c("","1","",""), size=3.5)+
  annotate("text",x=1,y=-1.17,label=c("","","2",""), size=3.5)+
  annotate("text",x=2,y=-0.17,label=c("","","3",""), size=3.5)+
  annotate("text",x=3,y=0.19,label=c("","","6",""), size=3.5)+
  annotate("text",x=2,y=-0.55,label=c("","","","1"), size=3.5)
Fig3j #Figure of the bootstrapped effect size of richeness in Continent


Fig3<-grid.arrange(arrangeGrob(arrangeGrob(Fig3a, Fig3b, Fig3c, Fig3d, Fig3e, nrow=5),
                               arrangeGrob(Fig3f, Fig3g, Fig3h, Fig3i, Fig3j, nrow=5),
                               widths=c(1,1)),
                               left=textGrob("Log Response ratio of species composition similarity \n (bootstrapped effect size)", 
                               				 rot=90, vjust=1),
                               top=textGrob(""),
                               right=textGrob(""),
                               bottom=textGrob("\nSuccessional stages"))

