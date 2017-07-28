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

#~~~~~~~~~~~~~~
#Birds
taxa<-subset(mydata,Taxa=="Birds")
#Exclude the with species richness 0 (the RRIn can't be calculated)
bir<-subset(taxa,RRIn!="NA")
nrow(bir)#45
summary(bir)
#Absolute Response Ratio
bir$RRIn2<-abs(bir$RRIn)
#Remove outlyers
outlyersbir<-boxplot(bir$RRIn2)
sort(outlyersbir$out) # >0.8
bir2<-subset(bir, bir$RRIn2<0.8)
boxplot(bir2$RRIn2)
bir2$RRIn2<-bir2$RRIn2+0.5

#Scale variables
library(dplyr)
bir3<- bir2 %>% mutate_each_(funs(scale(.) %>% as.vector),
							vars=c("Mean.Age","Mean_Elevation","Mean_Rainfall",
							"PatchPtoA5","PatchPtoA10","PatchPtoA25","PatchPtoA50","PatchPtoA75","PatchPtoA100",
							"PatchArea5","PatchArea10","PatchArea25","PatchArea50","PatchArea75","PatchArea100",
							"nPatches5","nPatches10","nPatches25","nPatches50","nPatches75","nPatches100",
							"p.Forest5km","p.Forest10km","p.Forest25km","p.Forest50km","p.Forest75km","p.Forest100km"))
summary(bir3)

### Testing correlation among covariates. From Ada Sanchez-Mercado.peronal communication (and then Marconi Cbiros)
#Being rs2 your stack of predictive variables
require(Hmisc)
names(bir3)
corrbir=bir3[,c(6,60:83,15,16,84)]
corrbir.t=data.matrix(corrbir)
clustbir = varclus(corrbir.t)
clustbir
plot(clustbir)
abline(h=0.3,lty=2,col=2)
#Chose variables with no Correlation

#First overall Selection of modesl, using glmulti
library(glmulti)
bir4 <- bir3[!apply(bir3[,c("SuccStage","Mean.Age","PreviousLandUse", 
				"PatchPtoA75","PatchPtoA25","PatchPtoA10",
				"nPatches5","nPatches50","nPatches75",
				"PatchArea10","PatchArea50",
				"p.Forest10km","p.Forest5km","p.Forest25km",
				"Continent_Island","Zoogeographical_Realm_Holt_etal_2013",
				"Mean_Elevation","Mean_Rainfall")], 1, anyNA),]
nrow(bir4)#n=34
## Good of fitness and normality
global.model<-glm(RRIn2 ~1+SuccStage+Mean.Age+PreviousLandUse+
				PatchPtoA75+PatchPtoA25+PatchPtoA10+
				nPatches5+nPatches50+nPatches75+
				PatchArea10+PatchArea50+
				p.Forest10km+p.Forest5km+p.Forest25km+
				Continent_Island+Zoogeographical_Realm_Holt_etal_2013+
				Mean_Elevation+Mean_Rainfall,data=bir4, Gamma("identity"))
plot(global.model)
hist(resid(global.model), breaks=20)
#1- Model selection
model <- glmulti(RRIn2 ~1+SuccStage+Mean.Age+PreviousLandUse+
				PatchPtoA75+PatchPtoA25+PatchPtoA10+
				nPatches5+nPatches50+nPatches75+
				PatchArea10+PatchArea50+
				p.Forest10km+p.Forest5km+p.Forest25km+
				Continent_Island+Zoogeographical_Realm_Holt_etal_2013+
				Mean_Elevation+Mean_Rainfall,data=bir4,
               level=1, family=Gamma (link="identity"), crit="aicc")
ranks<-weightable(model)
ranks
plot(model,type="s")
ranks7<-ranks[ranks$aicc <= min(ranks$aicc) + 7,]
ranks7
#11 models + null

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
bir5<-ddply(bir4,.(Site),randomRows,1)
m1<-glm(RRIn2~1, data=bir5, Gamma('identity'))
m2<-glm(RRIn2~1+SuccStage+PreviousLandUse+Continent_Island+Mean.Age+nPatches75+PatchArea50+p.Forest25km, data=bir5, Gamma('identity'))
m3<-glm(RRIn2~1+SuccStage+PreviousLandUse+Continent_Island+Mean.Age+nPatches75+PatchArea10+PatchArea50+p.Forest25km, data=bir5, Gamma('identity'))
m4<-glm(RRIn2~1+SuccStage+PreviousLandUse+Continent_Island+Mean.Age+PatchPtoA25+nPatches75+PatchArea50+p.Forest25km, data=bir5, Gamma('identity'))
m5<-glm(RRIn2~1+SuccStage+PreviousLandUse+Continent_Island+Mean.Age+PatchPtoA10+nPatches75+PatchArea50+p.Forest25km, data=bir5, Gamma('identity'))
m6<-glm(RRIn2~1+SuccStage+PreviousLandUse+Continent_Island+Mean.Age+PatchPtoA75+nPatches75+PatchArea50+p.Forest25km, data=bir5, Gamma('identity'))
m7<-glm(RRIn2~1+SuccStage+PreviousLandUse+Continent_Island+Mean.Age+nPatches5+nPatches75+PatchArea50+p.Forest25km, data=bir5, Gamma('identity'))
m8<-glm(RRIn2~1+SuccStage+PreviousLandUse+Continent_Island+Mean.Age+nPatches75+PatchArea50+p.Forest25km+Mean_Elevation, data=bir5, Gamma('identity'))
m9<-glm(RRIn2~1+SuccStage+PreviousLandUse+Continent_Island+Mean.Age+nPatches75+PatchArea50+p.Forest5km+p.Forest25km, data=bir5, Gamma('identity'))
m10<-glm(RRIn2~1+SuccStage+PreviousLandUse+Continent_Island+Mean.Age+nPatches50+nPatches75+PatchArea50+p.Forest25km, data=bir5, Gamma('identity'))
m11<-glm(RRIn2~1+SuccStage+PreviousLandUse+Continent_Island+Mean.Age+nPatches75+PatchArea50+p.Forest10km+p.Forest25km, data=bir5, Gamma('identity'))
m12<-glm(RRIn2~1+SuccStage+PreviousLandUse+Continent_Island+Mean.Age+nPatches75+PatchArea50+p.Forest25km+Mean_Rainfall, data=bir5, Gamma('identity'))

ranks<-m.AICc(list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12),
	as.integer(c("1","2","3","4","5","6","7","8","9","10","11","12")),
	length(bir5$RRIn2))

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
#m3 (Null)

#~~~~~~~~~~~~~~
#Mammals
taxa<-subset(mydata,Taxa=="Mammals")
#Exclude the with species richness 0 (the RRIn can't be calculated)
mam<-subset(taxa,RRIn!="NA")
nrow(mam)#25
summary(mam)
#Absolute Response Ratio
mam$RRIn2<-abs(mam$RRIn)
#Remove outlyers
outlyersmam<-boxplot(mam$RRIn2)
sort(outlyersmam$out) # >1.6
mam2<-mam
boxplot(mam2$RRIn2)
mam2$RRIn2<-mam2$RRIn2+0.5

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
				"PatchPtoA10","PatchPtoA100","PatchPtoA25","PatchPtoA50",
				"nPatches50","nPatches10",
				"PatchArea10","PatchArea50",
				"p.Forest5km",
				"Continent_Island","Zoogeographical_Realm_Holt_etal_2013",
				"Mean_Elevation","Mean_Rainfall")], 1, anyNA),]
nrow(mam4)#n=17
## Good of fitness and normality
global.model<-glm(RRIn2 ~1+SuccStage+Mean.Age+PreviousLandUse+
				PatchPtoA10+PatchPtoA100+PatchPtoA25+PatchPtoA50+
				nPatches50+nPatches10+
				PatchArea10+PatchArea50+
				p.Forest5km+
				Continent_Island+Zoogeographical_Realm_Holt_etal_2013+
				Mean_Elevation+Mean_Rainfall,data=mam4, Gamma("identity"))
plot(global.model)
hist(resid(global.model), breaks=20)
#1- Model selection
model <- glmulti(RRIn2 ~1+SuccStage+Mean.Age+PreviousLandUse+
				PatchPtoA10+PatchPtoA100+PatchPtoA25+PatchPtoA50+
				nPatches50+nPatches10+
				PatchArea10+PatchArea50+
				p.Forest5km+
				Continent_Island+Zoogeographical_Realm_Holt_etal_2013+
				Mean_Elevation+Mean_Rainfall,data=mam4,
               level=1, family=Gamma (link="identity"), crit="aicc")
ranks<-weightable(model)
ranks
plot(model,type="s")
ranks7<-ranks[ranks$aicc <= min(ranks$aicc) + 7,]
ranks7
#14 models + null

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
m1<-glm(RRIn2~1, data=mam5, Gamma('identity'))
m2<-glm(RRIn2~1+SuccStage+PatchPtoA25+PatchPtoA50+nPatches10, data=mam5, Gamma('identity'))
m3<-glm(RRIn2~1+SuccStage+PatchPtoA25+PatchPtoA50+nPatches50, data=mam5, Gamma('identity'))
m4<-glm(RRIn2~1+SuccStage+PreviousLandUse+PatchPtoA25+PatchPtoA50, data=mam5, Gamma('identity'))
m5<-glm(RRIn2~1+SuccStage+PatchPtoA100+PatchPtoA25, data=mam5, Gamma('identity'))
m6<-glm(RRIn2~1+SuccStage+PatchPtoA25+PatchPtoA50, data=mam5, Gamma('identity'))
m7<-glm(RRIn2~1+PatchPtoA25+nPatches50+nPatches10, data=mam5, Gamma('identity'))
m8<-glm(RRIn2~1+PatchPtoA25+PatchPtoA50+Mean_Elevation, data=mam5, Gamma('identity'))
m9<-glm(RRIn2~1+PatchPtoA25+PatchPtoA50, data=mam5, Gamma('identity'))
m10<-glm(RRIn2~1+SuccStage+PatchPtoA25+PatchPtoA50+PatchArea50, data=mam5, Gamma('identity'))
m11<-glm(RRIn2~1+SuccStage+PatchPtoA10+PatchPtoA25+PatchPtoA50+nPatches10, data=mam5, Gamma('identity'))
m12<-glm(RRIn2~1+SuccStage+PatchPtoA100+PatchPtoA25+nPatches10+PatchArea10, data=mam5, Gamma('identity'))
m13<-glm(RRIn2~1+SuccStage+PatchPtoA25+PatchPtoA50+p.Forest5km, data=mam5, Gamma('identity'))
m14<-glm(RRIn2~1+Mean.Age+PatchPtoA25+PatchPtoA50, data=mam5, Gamma('identity'))
m15<-glm(RRIn2~1+Mean.Age+PatchPtoA25+PatchPtoA50+Mean_Elevation, data=mam5, Gamma('identity'))

ranks<-m.AICc(list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15),
	as.integer(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15")),length(mam5$RRIn2))

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
#m4 (Null)
