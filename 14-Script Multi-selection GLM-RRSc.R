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
#Exclude the with species richness 0 (the RRSc can't be calculated)
bir<-subset(taxa,RRSc!="NA")
nrow(bir)#8
summary(bir)
#Absolute Response Ratio
bir$RRSc2<-abs(bir$RRSc)
#Remove outlyers
outlyersbir<-boxplot(bir$RRSc2)
sort(outlyersbir$out) # no
bir2<-bir
boxplot(bir2$RRSc2)
bir2$RRSc2<-bir2$RRSc2+0.5

#Scale variables
library(dplyr)
bir3<- bir2 %>% mutate_each_(funs(scale(.) %>% as.vector),
							vars=c("Mean.Age","Mean_Elevation","Mean_Rainfall",
							"PatchPtoA5","PatchPtoA10","PatchPtoA25","PatchPtoA50","PatchPtoA75","PatchPtoA100",
							"PatchArea5","PatchArea10","PatchArea25","PatchArea50","PatchArea75","PatchArea100",
							"nPatches5","nPatches10","nPatches25","nPatches50","nPatches75","nPatches100",
							"p.Forest5km","p.Forest10km","p.Forest25km","p.Forest50km","p.Forest75km","p.Forest100km"))
summary(bir3)

#### Testing correlation among covariates. From Ada Sanchez-Mercado.peronal communication (and then Marconi Cbiros)
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
				"nPatches50","nPatches75","nPatches5",
				"p.Forest25km","p.Forest10km",
				"PatchArea100","PatchArea5",
				"PatchPtoA5","PatchPtoA10","PatchPtoA100",
				"Continent_Island","Zoogeographical_Realm_Holt_etal_2013",
				"Mean_Elevation")], 1, anyNA),]
nrow(bir4)#n=5
## Good of fitness and normality
global.model<-glm(RRSc2 ~1+SuccStage+Mean.Age+PreviousLandUse+
				nPatches50+nPatches75+nPatches5+
				p.Forest25km+p.Forest10km+
				PatchArea100+PatchArea5+
				PatchPtoA5+PatchPtoA10+PatchPtoA100+
				Continent_Island+Zoogeographical_Realm_Holt_etal_2013+
				Mean_Elevation,data=bir4, Gamma("identity"))
plot(global.model)
hist(resid(global.model), breaks=20)
##Noooooooooo
#1- Model selection
model <- glmulti(RRSc2 ~1+SuccStage+Mean.Age+PreviousLandUse+
				p.Forest10km+
				PatchArea10+PatchArea100+
				nPatches25+
				PatchPtoA10+PatchPtoA100+PatchPtoA75+PatchPtoA25+
				Continent_Island+Zoogeographical_Realm_Holt_etal_2013+
				Mean_Elevation+Mean_Rainfall,data=bir4,
               level=1, family=Gamma (link="identity"), crit="aicc")
ranks<-weightable(model)
ranks
plot(model,type="s")
ranks7<-ranks[ranks$aicc <= min(ranks$aicc) + 7,]
ranks7
#70 models including null

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
m1<-glm(RRSc2~1, data=bir5, Gamma('identity'))
m2<-glm(RRSc2~1+p.Forest10km+PatchArea10+nPatches25+PatchPtoA25, data=bir5, Gamma('identity'))
m3<-glm(RRSc2~1+PreviousLandUse+Continent_Island+Mean_Rainfall, data=bir5, Gamma('identity'))
m4<-glm(RRSc2~1+p.Forest10km+PatchArea10+nPatches25+PatchPtoA25+Mean_Rainfall, data=bir5, Gamma('identity'))
m5<-glm(RRSc2~1+p.Forest10km+PatchArea10+PatchPtoA25, data=bir5, Gamma('identity'))
m6<-glm(RRSc2~1+PreviousLandUse+Continent_Island+Mean_Elevation+Mean_Rainfall, data=bir5, Gamma('identity'))
m7<-glm(RRSc2~1+p.Forest10km+PatchArea10+PatchPtoA100, data=bir5, Gamma('identity'))
m8<-glm(RRSc2~1+p.Forest10km+PatchArea10+nPatches25, data=bir5, Gamma('identity'))
m9<-glm(RRSc2~1+PreviousLandUse+Continent_Island+PatchPtoA75+Mean_Rainfall, data=bir5, Gamma('identity'))
m10<-glm(RRSc2~1+PreviousLandUse+Mean_Rainfall, data=bir5, Gamma('identity'))
m11<-glm(RRSc2~1+p.Forest10km+PatchArea10+PatchPtoA100+PatchPtoA25, data=bir5, Gamma('identity'))
m12<-glm(RRSc2~1+PreviousLandUse+Continent_Island, data=bir5, Gamma('identity'))
m13<-glm(RRSc2~1+p.Forest10km+PatchArea10+nPatches25+PatchPtoA10, data=bir5, Gamma('identity'))
m14<-glm(RRSc2~1+PreviousLandUse, data=bir5, Gamma('identity'))
m15<-glm(RRSc2~1+Continent_Island+p.Forest10km+PatchArea10+PatchArea100+nPatches25+Mean_Rainfall, data=bir5, Gamma('identity'))
m16<-glm(RRSc2~1+p.Forest10km+PatchArea10+nPatches25+Mean_Rainfall, data=bir5, Gamma('identity'))
m17<-glm(RRSc2~1+p.Forest10km+PatchArea10+PatchArea100+nPatches25+Mean_Rainfall, data=bir5, Gamma('identity'))
m18<-glm(RRSc2~1+Continent_Island+nPatches25, data=bir5, Gamma('identity'))
m19<-glm(RRSc2~1+p.Forest10km+PatchArea10+nPatches25+PatchPtoA100, data=bir5, Gamma('identity'))
m20<-glm(RRSc2~1+p.Forest10km+PatchArea10+nPatches25+PatchPtoA100+PatchPtoA25, data=bir5, Gamma('identity'))
m21<-glm(RRSc2~1+nPatches25, data=bir5, Gamma('identity'))
m22<-glm(RRSc2~1+Continent_Island+p.Forest10km+PatchArea10+nPatches25+PatchPtoA10+Mean_Rainfall, data=bir5, Gamma('identity'))
m23<-glm(RRSc2~1+PatchArea10, data=bir5, Gamma('identity'))
m24<-glm(RRSc2~1+Continent_Island+p.Forest10km+PatchArea10+nPatches25+PatchPtoA25, data=bir5, Gamma('identity'))
m25<-glm(RRSc2~1+Continent_Island+p.Forest10km+PatchArea10+nPatches25, data=bir5, Gamma('identity'))
m26<-glm(RRSc2~1+Continent_Island+p.Forest10km+PatchArea10+nPatches25+Mean_Rainfall, data=bir5, Gamma('identity'))
m27<-glm(RRSc2~1+p.Forest10km+PatchArea10+nPatches25+PatchPtoA10+Mean_Rainfall, data=bir5, Gamma('identity'))
m28<-glm(RRSc2~1+p.Forest10km+PatchArea10+nPatches25+PatchPtoA10+PatchPtoA25, data=bir5, Gamma('identity'))
m29<-glm(RRSc2~1+Mean.Age+p.Forest10km+PatchArea10+nPatches25+PatchPtoA25, data=bir5, Gamma('identity'))
m30<-glm(RRSc2~1+PreviousLandUse+Continent_Island+PatchArea100+Mean_Rainfall, data=bir5, Gamma('identity'))
m31<-glm(RRSc2~1+Continent_Island+p.Forest10km+PatchArea10+nPatches25+PatchPtoA10, data=bir5, Gamma('identity'))
m32<-glm(RRSc2~1+p.Forest10km+PatchArea10+nPatches25+PatchPtoA100+Mean_Rainfall, data=bir5, Gamma('identity'))
m33<-glm(RRSc2~1+p.Forest10km+PatchArea10+nPatches25+PatchPtoA25+Mean_Elevation, data=bir5, Gamma('identity'))
m34<-glm(RRSc2~1+p.Forest10km+PatchArea10+nPatches25+PatchPtoA75+PatchPtoA25, data=bir5, Gamma('identity'))
m35<-glm(RRSc2~1+p.Forest10km+PatchArea10+PatchArea100+nPatches25+PatchPtoA25, data=bir5, Gamma('identity'))
m36<-glm(RRSc2~1+p.Forest10km+PatchArea10+PatchPtoA10+PatchPtoA25, data=bir5, Gamma('identity'))
m37<-glm(RRSc2~1+p.Forest10km+PatchArea10+PatchPtoA25+Mean_Rainfall, data=bir5, Gamma('identity'))
m38<-glm(RRSc2~1+Continent_Island+p.Forest10km+PatchArea10+nPatches25+PatchPtoA25+Mean_Rainfall, data=bir5, Gamma('identity'))
m39<-glm(RRSc2~1+p.Forest10km+PatchArea10+PatchPtoA100+Mean_Rainfall, data=bir5, Gamma('identity'))
m40<-glm(RRSc2~1+Mean.Age+p.Forest10km+PatchArea10+PatchPtoA25, data=bir5, Gamma('identity'))
m41<-glm(RRSc2~1+Continent_Island+nPatches25+Mean_Rainfall, data=bir5, Gamma('identity'))
m42<-glm(RRSc2~1+PreviousLandUse+Continent_Island+nPatches25+Mean_Rainfall, data=bir5, Gamma('identity'))
m43<-glm(RRSc2~1+p.Forest10km+PatchArea10+PatchArea100+PatchPtoA25, data=bir5, Gamma('identity'))
m44<-glm(RRSc2~1+Continent_Island+p.Forest10km+PatchArea10+PatchPtoA25, data=bir5, Gamma('identity'))
m45<-glm(RRSc2~1+Continent_Island+PatchArea10, data=bir5, Gamma('identity'))
m46<-glm(RRSc2~1+p.Forest10km+PatchArea10+PatchPtoA75+PatchPtoA25, data=bir5, Gamma('identity'))
m47<-glm(RRSc2~1+p.Forest10km+PatchArea10+PatchPtoA25+Mean_Elevation, data=bir5, Gamma('identity'))
m48<-glm(RRSc2~1+p.Forest10km+PatchArea10+PatchArea100+nPatches25+PatchPtoA25+Mean_Rainfall, data=bir5, Gamma('identity'))
m49<-glm(RRSc2~1+p.Forest10km+nPatches25, data=bir5, Gamma('identity'))
m50<-glm(RRSc2~1+p.Forest10km+PatchArea10+PatchPtoA100+Mean_Elevation, data=bir5, Gamma('identity'))
m51<-glm(RRSc2~1+PatchArea10+PatchPtoA100, data=bir5, Gamma('identity'))
m52<-glm(RRSc2~1+p.Forest10km+PatchArea10, data=bir5, Gamma('identity'))
m53<-glm(RRSc2~1+p.Forest10km+PatchArea10+PatchPtoA100+PatchPtoA75, data=bir5, Gamma('identity'))
m54<-glm(RRSc2~1+p.Forest10km+PatchArea10+nPatches25+PatchPtoA100+PatchPtoA25+Mean_Rainfall, data=bir5, Gamma('identity'))
m55<-glm(RRSc2~1+PatchPtoA10, data=bir5, Gamma('identity'))
m56<-glm(RRSc2~1+p.Forest10km+PatchArea10+PatchPtoA10+PatchPtoA100, data=bir5, Gamma('identity'))
m57<-glm(RRSc2~1+Continent_Island+p.Forest10km+nPatches25, data=bir5, Gamma('identity'))
m58<-glm(RRSc2~1+PreviousLandUse+Continent_Island+PatchPtoA100+Mean_Rainfall, data=bir5, Gamma('identity'))
m59<-glm(RRSc2~1+PreviousLandUse+Continent_Island+p.Forest10km+Mean_Rainfall, data=bir5, Gamma('identity'))
m60<-glm(RRSc2~1+p.Forest10km+PatchArea10+PatchArea100+nPatches25, data=bir5, Gamma('identity'))
m61<-glm(RRSc2~1+Continent_Island+Zoogeographical_Realm_Holt_etal_2013+nPatches25+PatchPtoA75, data=bir5, Gamma('identity'))
m62<-glm(RRSc2~1+nPatches25+Mean_Rainfall, data=bir5, Gamma('identity'))
m63<-glm(RRSc2~1+Continent_Island+p.Forest10km+nPatches25+Mean_Rainfall, data=bir5, Gamma('identity'))
m64<-glm(RRSc2~1+PreviousLandUse+p.Forest10km, data=bir5, Gamma('identity'))
m65<-glm(RRSc2~1+PatchArea10+PatchPtoA25, data=bir5, Gamma('identity'))
m66<-glm(RRSc2~1+PatchPtoA25, data=bir5, Gamma('identity'))
m67<-glm(RRSc2~1+Continent_Island+p.Forest10km+PatchArea10+PatchPtoA100, data=bir5, Gamma('identity'))
m68<-glm(RRSc2~1+PatchArea10+nPatches25, data=bir5, Gamma('identity'))
m69<-glm(RRSc2~1+PatchPtoA100+PatchPtoA25, data=bir5, Gamma('identity'))
m70<-glm(RRSc2~1+PatchArea100, data=bir5, Gamma('identity'))


ranks<-m.AICc(list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,
					m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39,m40,
					m41,m42,m43,m44,m45,m46,m47,m48,m49,m50,m51,m52,m53,m54,m55,m56,m57,m58,m59,m60,
					m61,m62,m63,m64,m65,m66,m67,m68,m69,m70),
					as.integer(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20",
					"21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40",
					"41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60",
					"61","62","63","64","65","66","67","68","69","70")),
					length(bir5$RRSc2))

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
#m3 (null)

#~~~~~~~~~~~~~~
#Mammals
taxa<-subset(mydata,Taxa=="Mammals")
#Exclude the with species richness 0 (the RRSc can't be calculated)
mam<-subset(taxa,RRSc!="NA")
nrow(mam)#6
summary(mam)
#Absolute Response Ratio
mam$RRSc2<-abs(mam$RRSc)
#Remove outlyers
outlyersmam<-boxplot(mam$RRSc2)
sort(outlyersmam$out) #no
mam2<-mam
boxplot(mam2$RRSc2)
mam2$RRSc2<-mam2$RRSc2+0.5

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
				"nPatches100","nPatches10","nPatches75",
				"PatchPtoA100",
				"p.Forest5km","p.Forest75km","p.Forest10km",
				"PatchArea10",
				"Zoogeographical_Realm_Holt_etal_2013",
				"Mean_Elevation","Mean_Rainfall")], 1, anyNA),]
nrow(mam4)#n=5
## Good of fitness and normality
global.model<-glm(RRSc2 ~1+SuccStage+Mean.Age+PreviousLandUse+
				nPatches100+nPatches10+nPatches75+
				PatchPtoA100+
				p.Forest5km+p.Forest75km+p.Forest10km+
				PatchArea10+
				Mean_Elevation+Mean_Rainfall,data=mam4, Gamma("identity"))
plot(global.model)
hist(resid(global.model), breaks=20)
##No Genera residuales, ahhhhh

#1- Model selection
model <- glmulti(RRSc2 ~1+SuccStage+Mean.Age+PreviousLandUse+
				nPatches5+nPatches100+
				p.Forest75km+
				PatchArea10+PatchArea25+
				PatchPtoA50+PatchPtoA75+PatchPtoA10+PatchPtoA25+
				Continent_Island+Zoogeographical_Realm_Holt_etal_2013+
				Mean_Elevation+Mean_Rainfall,data=mam4,
               level=1, family=Gamma (link="identity"), crit="aicc")
ranks<-weightable(model)
ranks
plot(model,type="s")
ranks7<-ranks[ranks$aicc <= min(ranks$aicc) + 7,]
ranks7
#34 models + null

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
m1<-glm(RRSc2~1, data=mam5, Gamma('identity'))
m2<-glm(RRSc2~1+PreviousLandUse, data=mam5, Gamma('identity'))
m3<-glm(RRSc2~1, data=mam5, Gamma('identity'))
m4<-glm(RRSc2~1+PatchPtoA25, data=mam5, Gamma('identity'))
m5<-glm(RRSc2~1+PatchPtoA10, data=mam5, Gamma('identity'))
m6<-glm(RRSc2~1+Mean_Elevation, data=mam5, Gamma('identity'))
m7<-glm(RRSc2~1+nPatches100+Mean_Elevation, data=mam5, Gamma('identity'))
m8<-glm(RRSc2~1+PatchPtoA75, data=mam5, Gamma('identity'))
m9<-glm(RRSc2~1+PatchPtoA50, data=mam5, Gamma('identity'))
m10<-glm(RRSc2~1+nPatches100, data=mam5, Gamma('identity'))
m11<-glm(RRSc2~1+PatchArea25, data=mam5, Gamma('identity'))
m12<-glm(RRSc2~1+nPatches5, data=mam5, Gamma('identity'))
m13<-glm(RRSc2~1+Mean_Rainfall, data=mam5, Gamma('identity'))
m14<-glm(RRSc2~1+PatchArea10, data=mam5, Gamma('identity'))
m15<-glm(RRSc2~1+p.Forest75km, data=mam5, Gamma('identity'))
m16<-glm(RRSc2~1+Continent_Island, data=mam5, Gamma('identity'))
m17<-glm(RRSc2~1+Mean.Age, data=mam5, Gamma('identity'))
m18<-glm(RRSc2~1+nPatches100+p.Forest75km, data=mam5, Gamma('identity'))
m19<-glm(RRSc2~1+nPatches100+PatchPtoA25, data=mam5, Gamma('identity'))
m20<-glm(RRSc2~1+nPatches100+Mean_Rainfall, data=mam5, Gamma('identity'))
m21<-glm(RRSc2~1+PatchPtoA75+PatchPtoA25, data=mam5, Gamma('identity'))
m22<-glm(RRSc2~1+PatchArea10+PatchPtoA25, data=mam5, Gamma('identity'))
m23<-glm(RRSc2~1+PatchPtoA50+PatchPtoA25, data=mam5, Gamma('identity'))
m24<-glm(RRSc2~1+Continent_Island+PatchPtoA25, data=mam5, Gamma('identity'))
m25<-glm(RRSc2~1+PatchArea10+PatchPtoA10, data=mam5, Gamma('identity'))
m26<-glm(RRSc2~1+Continent_Island+Mean_Elevation, data=mam5, Gamma('identity'))
m27<-glm(RRSc2~1+PreviousLandUse+nPatches5, data=mam5, Gamma('identity'))
m28<-glm(RRSc2~1+nPatches100+PatchArea10, data=mam5, Gamma('identity'))
m29<-glm(RRSc2~1+nPatches5+PatchPtoA25, data=mam5, Gamma('identity'))
m30<-glm(RRSc2~1+PatchPtoA10+PatchPtoA25, data=mam5, Gamma('identity'))
m31<-glm(RRSc2~1+p.Forest75km+PatchPtoA25, data=mam5, Gamma('identity'))
m32<-glm(RRSc2~1+Mean.Age+PatchPtoA25, data=mam5, Gamma('identity'))
m33<-glm(RRSc2~1+PatchPtoA25+Mean_Rainfall, data=mam5, Gamma('identity'))
m34<-glm(RRSc2~1+PatchArea25+PatchPtoA25, data=mam5, Gamma('identity'))
m35<-glm(RRSc2~1+PatchPtoA25+Mean_Elevation, data=mam5, Gamma('identity'))

ranks<-m.AICc(list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,
					m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35),
					as.integer(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20",
					"21","22","23","24","25","26","27","28","29","30","31","32","33","34","35")),
					length(mam5$RRSc2))

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
m2b<-lm(RRSc2 ~ 1 + Continent_Island + PatchPtoA10 + Mean_Elevation, data=mam4)
summary(m2b)
m2b<-anova(m2b)
m2bss<-m2b$"Sum Sq"
print(cbind(m2b,PctExp=m2bss/sum(m2bss)*100)) #PctExp: % explain variation
