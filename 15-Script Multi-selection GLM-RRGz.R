#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Multimodel selection and key predictors
#Clear memory
rm(list=ls())
#Working in the office
setwd("//Users/orlando/Dropbox/Vertebrate Recovery - Nature/data")
#Call the file MainData.csv
mydata <- read.csv("//Users/orlando/Dropbox/Vertebrate Recovery - Nature/data/MainData.csv")

#Working in my house
setwd("//Users/mac/Dropbox/Vertebrate Recovery - Nature/data")
#Call the file MainData.csv
mydata <- read.csv("//Users/mac/Dropbox/Vertebrate Recovery - Nature/data/MainData.csv")

#~~~~~~~~~~~~~~No birds data
#Mammals
taxa<-subset(mydata,Taxa=="Mammals")
#Exclude the with species richness 0 (the RRGz can't be calculated)
mam<-subset(taxa,RRGz!="NA")
nrow(mam)#21
summary(mam)
#Absolute Response Ratio
mam$RRGz2<-abs(mam$RRGz)
#Remove outlyers
outlyersmam<-boxplot(mam$RRGz2)
sort(outlyersmam$out) # >1.7
mam2<-subset(mam, mam$RRGz2<1.7)
boxplot(mam2$RRGz2)
mam2$RRGz2<-mam2$RRGz2+0.5

#Scale variables
library(dplyr)
mam3<- mam2 %>% mutate_each_(funs(scale(.) %>% as.vector),
							vars=c("Mean.Age","Mean_Elevation","Mean_Rainfall",
							"PatchPtoA5","PatchPtoA10","PatchPtoA25","PatchPtoA50","PatchPtoA75","PatchPtoA100",
							"PatchArea5","PatchArea10","PatchArea25","PatchArea50","PatchArea75","PatchArea100",
							"nPatches5","nPatches10","nPatches25","nPatches50","nPatches75","nPatches100",
							"p.Forest5km","p.Forest10km","p.Forest25km","p.Forest50km","p.Forest75km","p.Forest100km"))
summary(mam3)

### Testing correlation among covariates. From Ada Sanchez-Mercado.peronal communication (and then Marconi Cmamos)
#Being rs2 your stack of predictive variables
require(Hmisc)
names(mam3)
corrmam=mam3[,c(6,60:83,15,16,84)]
corrmam.t=data.matrix(corrmam)
clustmam = varclus(corrmam.t)
clustmam
plot(clustmam)
abline(h=0.3,lty=2,col=2)
#Chose variables with no Correlation

#First overall Selection of models, using glmulti
library(glmulti)
mam4 <- mam3[!apply(mam3[,c("SuccStage","Mean.Age","PreviousLandUse", 
				"PatchArea10","PatchArea100",
				"nPatches5","nPatches25",
				"PatchPtoA5","PatchPtoA10","PatchPtoA25","PatchPtoA100","PatchPtoA50",
				"p.Forest50km",
				"Continent_Island","Zoogeographical_Realm_Holt_etal_2013",
				"Mean_Elevation","Mean_Rainfall")], 1, anyNA),]
nrow(mam4)#n=18
## Good of fitness and normality
global.model<-glm(RRGz2 ~1+SuccStage+Mean.Age+PreviousLandUse+
				PatchArea10+PatchArea100+
				nPatches5+nPatches25+
				PatchPtoA5+PatchPtoA10+PatchPtoA25+PatchPtoA100+PatchPtoA50+
				p.Forest50km+
				Continent_Island+Zoogeographical_Realm_Holt_etal_2013+
				Mean_Elevation+Mean_Rainfall,data=mam4, Gamma("identity"))
plot(global.model)
hist(resid(global.model), breaks=20)
#1- Model selection
model <- glmulti(RRGz2 ~1+SuccStage+Mean.Age+PreviousLandUse+
				PatchArea10+PatchArea100+
				nPatches5+nPatches25+
				PatchPtoA5+PatchPtoA10+PatchPtoA25+PatchPtoA100+PatchPtoA50+
				p.Forest50km+
				Continent_Island+Zoogeographical_Realm_Holt_etal_2013+
				Mean_Elevation+Mean_Rainfall,data=mam4,
               level=1, family=Gamma (link="identity"), crit="aicc")
ranks<-weightable(model)
ranks
plot(model,type="s")
ranks7<-ranks[ranks$aicc <= min(ranks$aicc) + 7,]
ranks7
#20 models + null

#To select one comparison for each study
randomRows = function(df,n){
  return(df[sample(nrow(df),n),])
}

##2 - Model selection - AICc
m.AICc <- function(modelos,data,n){
  LL <- sapply(modelos,logLik)
  k <- sapply(lapply(modelos,logLik),attr,"df")
  AIC <- -2*LL+2*k
  AICc <- AIC+((2*k*(k+1))/(n-k-1))
  d.AICc <- AICc-min(AICc)
  w.AICc <- (exp(-0.5*d.AICc))/sum(exp(-0.5*d.AICc))
  data<-data
  data.frame(n.par=k,log.lik=LL,AICc=AICc,AIC=AIC,delta.AICc=d.AICc,w.AICc=w.AICc,name=data,row.names=names(LL))[order(d.AICc),]
}

library(plyr)
## Run models via bootstrap
top_ranked<-list()
null_model_AICc<-list()
null_model_position<-list()
#coef<-list()
#r2<-list()
for (i in 1:10000){
## sample 1 data per study
mam5<-ddply(mam4,.(Site),randomRows,1)
m1<-glm(RRGz2~1, data=mam5, Gamma('identity'))
m2<-glm(RRGz2~1+Continent_Island+PatchPtoA5, data=mam5, Gamma('identity'))
m3<-glm(RRGz2~1+Continent_Island+Mean.Age+PatchPtoA5, data=mam5, Gamma('identity'))
m4<-glm(RRGz2~1+Continent_Island+PatchArea10+PatchPtoA5, data=mam5, Gamma('identity'))
m5<-glm(RRGz2~1+Continent_Island+PatchPtoA5+Mean_Rainfall, data=mam5, Gamma('identity'))
m6<-glm(RRGz2~1+Continent_Island+PatchPtoA5+PatchPtoA10, data=mam5, Gamma('identity'))
m7<-glm(RRGz2~1+Continent_Island+nPatches25, data=mam5, Gamma('identity'))
m8<-glm(RRGz2~1+Continent_Island+PatchPtoA5+p.Forest50km, data=mam5, Gamma('identity'))
m9<-glm(RRGz2~1+Continent_Island+PatchPtoA5+PatchPtoA100, data=mam5, Gamma('identity'))
m10<-glm(RRGz2~1+Continent_Island+PatchArea100+PatchPtoA5, data=mam5, Gamma('identity'))
m11<-glm(RRGz2~1+Continent_Island+PatchPtoA5+PatchPtoA50, data=mam5, Gamma('identity'))
m12<-glm(RRGz2~1+Continent_Island+nPatches25+PatchPtoA5, data=mam5, Gamma('identity'))
m13<-glm(RRGz2~1+Continent_Island+nPatches5+PatchPtoA5, data=mam5, Gamma('identity'))
m14<-glm(RRGz2~1+Continent_Island+PatchPtoA5+Mean_Elevation, data=mam5, Gamma('identity'))
m15<-glm(RRGz2~1+Continent_Island+PatchPtoA5+PatchPtoA25, data=mam5, Gamma('identity'))
m16<-glm(RRGz2~1+Continent_Island+p.Forest50km, data=mam5, Gamma('identity'))
m17<-glm(RRGz2~1+Continent_Island+nPatches25+Mean_Rainfall, data=mam5, Gamma('identity'))
m18<-glm(RRGz2~1+p.Forest50km+Mean_Elevation, data=mam5, Gamma('identity'))
m19<-glm(RRGz2~1+Continent_Island+PatchArea10, data=mam5, Gamma('identity'))
m20<-glm(RRGz2~1+Continent_Island+nPatches5, data=mam5, Gamma('identity'))
m21<-glm(RRGz2~1+Continent_Island, data=mam5, Gamma('identity'))

ranks<-m.AICc(list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21),
	as.integer(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21")),
	length(mam5$RRGz2))

top_ranked[[i]]<-ranks[1,]
null_model_AICc[[i]]<-ranks[ranks$name=="1",]
null_model_position[[i]]<-which(ranks$name=="1",)
#coef[[i]]<-list(summary(m))
#r2[[i]]<-r.squaredLR(m, null = m1)
}

## Transforming lists in data.frames
df_top_ranked<-data.frame(matrix(unlist(top_ranked), nrow=10000, byrow=T))
df_top_ranked<-data.frame(matrix(unlist(top_ranked), nrow=10000, byrow=T))
df_null_model_AICc<-data.frame(matrix(unlist(null_model_AICc), nrow=10000, byrow=T))
df_null_model_position<-data.frame(matrix(unlist(null_model_position), nrow=10000, byrow=T))

## Saving data.frames
#write.table(df_top_ranked, "Top_ranked.txt",quote=F)
#write.table(df_null_model_AICc, "Null_model_AICc.txt",quote=F)
#write.table(df_null_model_position, "Null_model_position.txt",quote=F)
#sink("Coefficients.txt") 
#lapply(coef, print) 
#sink()
#sink("R2.txt") 
#lapply(r2, print) 
#sink() 

## Top-model frequence
plausible_mod_frequency<-as.data.frame(ftable(df_top_ranked[,7]))
plausible_mod_frequency
#m2 (RRGz2~1+Continent_Island+PatchPtoA5)
m2b<-lm(RRGz2~1+Continent_Island+PatchPtoA5, data=mam4)
summary(m2b)
Am2b<-anova(m2b)
Am2bss<-Am2b$"Sum Sq"
print(cbind(Am2b,PctExp=Am2bss/sum(Am2bss)*100)) #PctExp: % explain variation
