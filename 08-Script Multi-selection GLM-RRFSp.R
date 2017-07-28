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
#Amphibians
taxa<-subset(mydata,Taxa=="Amphibians")
#Exclude the with species richness 0 (the RRFSp can't be calculated)
amp<-subset(taxa,RRFSp!="NA")
nrow(amp)#37
summary(amp)
#Absolute Response Ratio
amp$RRFSp2<-abs(amp$RRFSp)
#Remove outlyers
outlyersamp<-boxplot(amp$RRFSp2)
sort(outlyersamp$out) # >3
amp2<-subset(amp, amp$RRFSp2<3)
boxplot(amp2$RRFSp2)
amp2$RRFSp2<-amp2$RRFSp2+0.5

#Scale variables
library(dplyr)
amp3<- amp2 %>% mutate_each_(funs(scale(.) %>% as.vector),
							vars=c("Mean.Age","Mean_Elevation","Mean_Rainfall",
							"PatchPtoA5","PatchPtoA10","PatchPtoA25","PatchPtoA50","PatchPtoA75","PatchPtoA100",
							"PatchArea5","PatchArea10","PatchArea25","PatchArea50","PatchArea75","PatchArea100",
							"nPatches5","nPatches10","nPatches25","nPatches50","nPatches75","nPatches100",
							"p.Forest5km","p.Forest10km","p.Forest25km","p.Forest50km","p.Forest75km","p.Forest100km"))
summary(amp3)

### Testing correlation among covariates. From Ada Sanchez-Mercado.peronal communication (and then Marconi Campos)
#Being rs2 your stack of predictive variables
require(Hmisc)
names(amp3)
corramp=amp3[,c(6,60:83,15,16,84)]
corramp.t=data.matrix(corramp)
clustamp = varclus(corramp.t)
clustamp
plot(clustamp)
abline(h=0.3,lty=2,col=2)
#Chose variables with no Correlation
summary(amp3)
#First overall Selection of modesl, using glmulti
library(glmulti)
amp4 <- amp3[!apply(amp3[,c("SuccStage","Mean.Age","PreviousLandUse", 
				"nPatches100","nPatches5",
				"PatchPtoA10","PatchPtoA25",
				"p.Forest10km",
				"PatchArea100","PatchArea25",
				"Continent_Island","Zoogeographical_Realm_Holt_etal_2013",
				"Mean_Elevation","Mean_Rainfall")], 1, anyNA),]
nrow(amp4)#n=26
## Good of fitness and normality
global.model<-glm(RRFSp2 ~1+SuccStage+Mean.Age+PreviousLandUse+
				nPatches100+nPatches5+
				PatchPtoA10+PatchPtoA25+
				PatchArea25+PatchArea100+
				p.Forest10km+
				Continent_Island+Zoogeographical_Realm_Holt_etal_2013+
				Mean_Elevation+Mean_Rainfall,data=amp4, Gamma("identity"))
plot(global.model)
hist(resid(global.model), breaks=20)
#1- Model selection
model <- glmulti(RRFSp2 ~1+SuccStage+Mean.Age+PreviousLandUse+
				nPatches100+nPatches5+
				PatchPtoA10+PatchPtoA25+
				PatchArea25+PatchArea100+
				p.Forest10km+
				Continent_Island+Zoogeographical_Realm_Holt_etal_2013+
				Mean_Elevation+Mean_Rainfall,data=amp4,
               level=1, family=Gamma (link="identity"), crit="aicc")
ranks<-weightable(model)
ranks
plot(model,type="s")
ranks7<-ranks[ranks$aicc <= min(ranks$aicc) + 7,]
ranks7
#16 models + null

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
## Sample 1 data per study
amp5<-ddply(amp4,.(Site),randomRows,1)
m1<-glm(RRFSp2~1, data=amp5, Gamma('identity'))
m2<-glm(RRFSp2~1+Continent_Island+Mean.Age+PatchPtoA25+PatchArea25+p.Forest100km+Mean_Rainfall, data=amp5, Gamma('identity'))
m3<-glm(RRFSp2~1+Continent_Island+PatchPtoA25+PatchArea25+p.Forest100km+Mean_Rainfall, data=amp5, Gamma('identity'))
m4<-glm(RRFSp2~1+Continent_Island+nPatches5+PatchPtoA25+PatchArea25+p.Forest100km+Mean_Rainfall, data=amp5, Gamma('identity'))
m5<-glm(RRFSp2~1+Continent_Island+nPatches100+PatchPtoA25+PatchArea25+p.Forest100km+Mean_Rainfall, data=amp5, Gamma('identity'))
m6<-glm(RRFSp2~1+Continent_Island+PatchPtoA25+PatchArea25+p.Forest100km+Mean_Elevation+Mean_Rainfall, data=amp5, Gamma('identity'))
m7<-glm(RRFSp2~1+Continent_Island+Mean.Age+nPatches5+PatchPtoA25+PatchArea25+p.Forest100km+Mean_Rainfall, data=amp5, Gamma('identity'))
m8<-glm(RRFSp2~1+Continent_Island+Mean.Age+PatchPtoA10+PatchPtoA25+PatchArea25+p.Forest100km+Mean_Rainfall, data=amp5, Gamma('identity'))
m9<-glm(RRFSp2~1+Continent_Island+Mean.Age+nPatches100+PatchPtoA25+PatchArea25+p.Forest100km+Mean_Rainfall, data=amp5, Gamma('identity'))
m10<-glm(RRFSp2~1+Continent_Island+PatchPtoA10+PatchPtoA25+PatchArea25+p.Forest100km+Mean_Rainfall, data=amp5, Gamma('identity'))
m11<-glm(RRFSp2~1+Continent_Island+Mean.Age+PatchPtoA25+PatchArea25+p.Forest100km+Mean_Elevation+Mean_Rainfall, data=amp5, Gamma('identity'))
m12<-glm(RRFSp2~1+Continent_Island+nPatches100+PatchPtoA25+PatchArea25+p.Forest100km+Mean_Elevation+Mean_Rainfall, data=amp5, Gamma('identity'))
m13<-glm(RRFSp2~1+Continent_Island+nPatches100+nPatches5+PatchPtoA25+PatchArea25+p.Forest100km+Mean_Rainfall, data=amp5, Gamma('identity'))
m14<-glm(RRFSp2~1+Continent_Island+nPatches5+PatchPtoA25+PatchArea25+p.Forest100km+Mean_Elevation+Mean_Rainfall, data=amp5, Gamma('identity'))
m15<-glm(RRFSp2~1+Continent_Island+Mean.Age+PatchPtoA10+PatchPtoA25+PatchArea25+Mean_Rainfall, data=amp5, Gamma('identity'))
m16<-glm(RRFSp2~1+Continent_Island+nPatches5+PatchPtoA10+PatchPtoA25+PatchArea25+p.Forest100km+Mean_Rainfall, data=amp5, Gamma('identity'))
m17<-glm(RRFSp2~1+Continent_Island+nPatches100+PatchPtoA10+PatchPtoA25+PatchArea25+p.Forest100km+Mean_Rainfall, data=amp5, Gamma('identity'))


ranks<-m.AICc(list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17),
				as.integer(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")),
				length(amp5$RRFSp2))

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
#m1 (Null)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Reptiles
taxa<-subset(mydata,Taxa=="Reptiles")
#Exclude the with species richness 0 (the RRFSp can't be calculated)
rep<-subset(taxa,RRFSp!="NA")
nrow(rep)#24
summary(rep)
#Absolute Response Ratio
rep$RRFSp2<-abs(rep$RRFSp)
#Remove outlyers
outlyersrep<-boxplot(rep$RRFSp2)
sort(outlyersrep$out) # 2.1
rep2<-subset(rep, rep$RRFSp2<2.1)
boxplot(rep2$RRFSp2)
rep2$RRFSp2<-rep2$RRFSp2+0.5

#Scale variables
library(dplyr)
rep3<- rep2 %>% mutate_each_(funs(scale(.) %>% as.vector),
							vars=c("Mean.Age","Mean_Elevation","Mean_Rainfall",
							"PatchPtoA5","PatchPtoA10","PatchPtoA25","PatchPtoA50","PatchPtoA75","PatchPtoA100",
							"PatchArea5","PatchArea10","PatchArea25","PatchArea50","PatchArea75","PatchArea100",
							"nPatches5","nPatches10","nPatches25","nPatches50","nPatches75","nPatches100",
							"p.Forest5km","p.Forest10km","p.Forest25km","p.Forest50km","p.Forest75km","p.Forest100km"))
summary(rep3)

### Testing correlation among covariates. From Ada Sanchez-Mercado.peronal communication (and then Marconi Crepos)
#Being rs2 your stack of predictive variables
require(Hmisc)
names(rep3)
coRRFSpep=rep3[,c(6,60:83,15,16,84)]
coRRFSpep.t=data.matrix(coRRFSpep)
clustrep = varclus(coRRFSpep.t)
clustrep
plot(clustrep)
abline(h=0.3,lty=2,col=2)
#Chose variables with no Correlation

#First overall Selection of modesl, using glmulti
library(glmulti)
rep4 <- rep3[!apply(rep3[,c("SuccStage","Mean.Age","PreviousLandUse", 
				"nPatches100","nPatches50",
				"PatchPtoA10","PatchPtoA5","PatchPtoA100","PatchPtoA50","PatchPtoA25","PatchPtoA75",
				"PatchArea75","PatchArea50","PatchArea25",
				"Continent_Island","Zoogeographical_Realm_Holt_etal_2013",
				"Mean_Elevation","Mean_Rainfall")], 1, anyNA),]
nrow(rep4)#n=20
## Good of fitness and normality
global.model<-glm(RRFSp2 ~1+SuccStage+Mean.Age+PreviousLandUse+
				nPatches100+nPatches50+
				PatchPtoA10+PatchPtoA5+PatchPtoA100+PatchPtoA50+PatchPtoA25+PatchPtoA75+
				PatchArea75+PatchArea50+PatchArea25+
				Continent_Island+Zoogeographical_Realm_Holt_etal_2013+
				Mean_Elevation+Mean_Rainfall,data=rep4, Gamma("identity"))
plot(global.model)
hist(resid(global.model), breaks=20)
#1-Model selection
model <- glmulti(RRFSp2 ~1+SuccStage+Mean.Age+PreviousLandUse+
				nPatches100+nPatches50+
				PatchPtoA10+PatchPtoA5+PatchPtoA100+PatchPtoA50+PatchPtoA25+PatchPtoA75+
				PatchArea75+PatchArea50+PatchArea25+
				Continent_Island+Zoogeographical_Realm_Holt_etal_2013+
				Mean_Elevation+Mean_Rainfall,data=rep4,
               level=1, family=Gamma (link="identity"), crit="aicc")
ranks<-weightable(model)
ranks
plot(model,type="s")
ranks7<-ranks[ranks$aicc <= min(ranks$aicc) + 7,]
ranks7
#76 models + null

randomRows = function(df,n){
  return(df[sample(nrow(df),n),])
}
## Model selection - AICc
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
## Sample 1 data per study
rep5<-ddply(rep4,.(Site),randomRows,1)

m1<-glm(RRFSp2~1, data=rep5, Gamma('identity'))
m2<-glm(RRFSp2~1+SuccStage+PatchArea75, data=rep5, Gamma('identity'))
m3<-glm(RRFSp2~1+SuccStage+PatchPtoA100+PatchArea50, data=rep5, Gamma('identity'))
m4<-glm(RRFSp2~1+SuccStage+PatchPtoA50+PatchArea75+PatchArea50, data=rep5, Gamma('identity'))
m5<-glm(RRFSp2~1+SuccStage+PatchPtoA75+Mean_Elevation, data=rep5, Gamma('identity'))
m6<-glm(RRFSp2~1+SuccStage+PatchPtoA100+PatchArea75, data=rep5, Gamma('identity'))
m7<-glm(RRFSp2~1+SuccStage+PatchPtoA50+PatchArea75, data=rep5, Gamma('identity'))
m8<-glm(RRFSp2~1+SuccStage+nPatches50+PatchArea75, data=rep5, Gamma('identity'))
m9<-glm(RRFSp2~1+SuccStage+PatchPtoA100+PatchPtoA50+PatchArea75, data=rep5, Gamma('identity'))
m10<-glm(RRFSp2~1+SuccStage+PatchPtoA75+PatchArea50, data=rep5, Gamma('identity'))
m11<-glm(RRFSp2~1+SuccStage+nPatches100+PatchArea75, data=rep5, Gamma('identity'))
m12<-glm(RRFSp2~1+SuccStage+PatchArea75+Mean_Elevation, data=rep5, Gamma('identity'))
m13<-glm(RRFSp2~1+SuccStage+Continent_Island+PatchArea75, data=rep5, Gamma('identity'))
m14<-glm(RRFSp2~1+SuccStage+PatchPtoA50+PatchPtoA75+PatchArea75, data=rep5, Gamma('identity'))
m15<-glm(RRFSp2~1+SuccStage+PatchPtoA25+PatchArea75, data=rep5, Gamma('identity'))
m16<-glm(RRFSp2~1+SuccStage+PatchPtoA75+PatchArea75, data=rep5, Gamma('identity'))
m17<-glm(RRFSp2~1+SuccStage+nPatches50, data=rep5, Gamma('identity'))
m18<-glm(RRFSp2~1+SuccStage+nPatches100, data=rep5, Gamma('identity'))
m19<-glm(RRFSp2~1+SuccStage+PatchArea50, data=rep5, Gamma('identity'))
m20<-glm(RRFSp2~1+SuccStage+PatchArea75+PatchArea25, data=rep5, Gamma('identity'))
m21<-glm(RRFSp2~1+SuccStage+Continent_Island+PatchArea50, data=rep5, Gamma('identity'))
m22<-glm(RRFSp2~1+SuccStage+PatchPtoA75+Mean_Elevation+Mean_Rainfall, data=rep5, Gamma('identity'))
m23<-glm(RRFSp2~1+SuccStage+PatchPtoA25+PatchPtoA75+PatchArea50, data=rep5, Gamma('identity'))
m24<-glm(RRFSp2~1+SuccStage+PatchArea75+PatchArea50, data=rep5, Gamma('identity'))
m25<-glm(RRFSp2~1+SuccStage+PatchPtoA10+PatchArea75, data=rep5, Gamma('identity'))
m26<-glm(RRFSp2~1+SuccStage+PatchArea75+Mean_Rainfall, data=rep5, Gamma('identity'))
m27<-glm(RRFSp2~1+SuccStage+PatchPtoA5+PatchArea75, data=rep5, Gamma('identity'))
m28<-glm(RRFSp2~1+SuccStage+Mean.Age+PatchArea75, data=rep5, Gamma('identity'))
m29<-glm(RRFSp2~1+SuccStage+Continent_Island+PatchPtoA50+PatchArea75, data=rep5, Gamma('identity'))
m30<-glm(RRFSp2~1+SuccStage+PatchPtoA50+PatchPtoA75+PatchArea50, data=rep5, Gamma('identity'))
m31<-glm(RRFSp2~1+SuccStage+nPatches100+PatchArea50, data=rep5, Gamma('identity'))
m32<-glm(RRFSp2~1+SuccStage+Mean_Elevation, data=rep5, Gamma('identity'))
m33<-glm(RRFSp2~1+SuccStage+PatchPtoA75+PatchArea75+Mean_Elevation, data=rep5, Gamma('identity'))
m34<-glm(RRFSp2~1+SuccStage+PatchPtoA100+PatchPtoA75+Mean_Elevation, data=rep5, Gamma('identity'))
m35<-glm(RRFSp2~1+SuccStage+PatchPtoA75+PatchArea50+Mean_Elevation, data=rep5, Gamma('identity'))
m36<-glm(RRFSp2~1+SuccStage+nPatches50+PatchArea50, data=rep5, Gamma('identity'))
m37<-glm(RRFSp2~1+SuccStage+PatchPtoA5+PatchPtoA75, data=rep5, Gamma('identity'))
m38<-glm(RRFSp2~1+SuccStage+PatchPtoA100+PatchArea75+Mean_Elevation, data=rep5, Gamma('identity'))
m39<-glm(RRFSp2~1+SuccStage+nPatches50+PatchPtoA75, data=rep5, Gamma('identity'))
m40<-glm(RRFSp2~1+SuccStage+PatchArea75+PatchArea50+Mean_Elevation, data=rep5, Gamma('identity'))
m41<-glm(RRFSp2~1+SuccStage+PatchPtoA25+PatchArea75+PatchArea25, data=rep5, Gamma('identity'))
m42<-glm(RRFSp2~1+SuccStage+PatchPtoA100+PatchArea75+PatchArea50, data=rep5, Gamma('identity'))
m43<-glm(RRFSp2~1+SuccStage+Mean_Elevation+Mean_Rainfall, data=rep5, Gamma('identity'))
m44<-glm(RRFSp2~1+SuccStage+nPatches50+PatchPtoA25+PatchArea75, data=rep5, Gamma('identity'))
m45<-glm(RRFSp2~1+SuccStage+PatchPtoA100+PatchArea50+Mean_Rainfall, data=rep5, Gamma('identity'))
m46<-glm(RRFSp2~1+SuccStage+Continent_Island+PatchPtoA100+PatchArea50, data=rep5, Gamma('identity'))
m47<-glm(RRFSp2~1+SuccStage+PatchPtoA100+PatchPtoA75+PatchArea50, data=rep5, Gamma('identity'))
m48<-glm(RRFSp2~1+SuccStage+PatchPtoA75+PatchArea50+PatchArea25, data=rep5, Gamma('identity'))
m49<-glm(RRFSp2~1+SuccStage+PatchPtoA75+PatchArea25+Mean_Elevation, data=rep5, Gamma('identity'))
m50<-glm(RRFSp2~1+SuccStage+PatchPtoA25+Mean_Elevation, data=rep5, Gamma('identity'))
m51<-glm(RRFSp2~1+SuccStage+Mean.Age+PatchPtoA100+PatchArea50, data=rep5, Gamma('identity'))
m52<-glm(RRFSp2~1+SuccStage+PatchPtoA25+PatchPtoA75+Mean_Elevation, data=rep5, Gamma('identity'))
m53<-glm(RRFSp2~1+SuccStage+PatchPtoA5+PatchPtoA100+PatchArea75, data=rep5, Gamma('identity'))
m54<-glm(RRFSp2~1+SuccStage+PatchPtoA50+PatchPtoA25+PatchArea75, data=rep5, Gamma('identity'))
m55<-glm(RRFSp2~1+SuccStage+nPatches100+PatchArea25, data=rep5, Gamma('identity'))
m56<-glm(RRFSp2~1+SuccStage+nPatches50+PatchPtoA100+PatchArea50, data=rep5, Gamma('identity'))
m57<-glm(RRFSp2~1+SuccStage+PatchArea50+Mean_Elevation, data=rep5, Gamma('identity'))
m58<-glm(RRFSp2~1+SuccStage+PatchPtoA10+PatchPtoA100+PatchArea50, data=rep5, Gamma('identity'))
m59<-glm(RRFSp2~1+SuccStage+PatchPtoA50+PatchArea75+PatchArea25, data=rep5, Gamma('identity'))
m60<-glm(RRFSp2~1+SuccStage+PatchPtoA75, data=rep5, Gamma('identity'))
m61<-glm(RRFSp2~1+SuccStage+Continent_Island+PatchArea75+Mean_Elevation, data=rep5, Gamma('identity'))
m62<-glm(RRFSp2~1+SuccStage+nPatches100+PatchPtoA25, data=rep5, Gamma('identity'))
m63<-glm(RRFSp2~1+SuccStage+nPatches100+PatchPtoA75, data=rep5, Gamma('identity'))
m64<-glm(RRFSp2~1+SuccStage+nPatches50+PatchPtoA100+PatchArea75, data=rep5, Gamma('identity'))
m65<-glm(RRFSp2~1+SuccStage+nPatches100+PatchPtoA50+PatchArea75, data=rep5, Gamma('identity'))
m66<-glm(RRFSp2~1+SuccStage+PatchPtoA100+PatchPtoA50+PatchArea50, data=rep5, Gamma('identity'))
m67<-glm(RRFSp2~1+SuccStage+PatchPtoA100+PatchArea50+Mean_Elevation, data=rep5, Gamma('identity'))
m68<-glm(RRFSp2~1+SuccStage+PatchPtoA5+PatchPtoA100+PatchArea50, data=rep5, Gamma('identity'))
m69<-glm(RRFSp2~1+SuccStage+nPatches100+PatchPtoA100+PatchArea50, data=rep5, Gamma('identity'))
m70<-glm(RRFSp2~1+SuccStage+PatchPtoA100+PatchArea75+Mean_Rainfall, data=rep5, Gamma('identity'))
m71<-glm(RRFSp2~1+SuccStage+Mean.Age+PatchPtoA75+Mean_Elevation, data=rep5, Gamma('identity'))
m72<-glm(RRFSp2~1+SuccStage+PatchPtoA10+PatchPtoA50+PatchArea75, data=rep5, Gamma('identity'))
m73<-glm(RRFSp2~1+SuccStage+PatchArea25+Mean_Elevation, data=rep5, Gamma('identity'))
m74<-glm(RRFSp2~1+SuccStage+PatchPtoA100+PatchPtoA25+PatchArea50, data=rep5, Gamma('identity'))
m75<-glm(RRFSp2~1+SuccStage+PatchPtoA100+PatchArea50+PatchArea25, data=rep5, Gamma('identity'))
m76<-glm(RRFSp2~1+SuccStage+PatchPtoA50+Mean_Elevation, data=rep5, Gamma('identity'))
m77<-glm(RRFSp2~1+SuccStage+Continent_Island+PatchPtoA75+PatchArea50, data=rep5, Gamma('identity'))


ranks<-m.AICc(list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,
					m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39,m40,
					m41,m42,m43,m44,m45,m46,m47,m48,m49,m50,m51,m52,m53,m54,m55,m56,m57,m58,m59,m60,
					m61,m62,m63,m64,m65,m66,m67,m68,m69,m70,m71,m72,m73,m74,m75,m76,m77),
					as.integer(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20",
					"21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40",
					"41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60",
					"61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77")),
					length(rep5$RRFSp2))

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
#Null

#~~~~~~~~~~~~~~
#Birds
taxa<-subset(mydata,Taxa=="Birds")
#Exclude the with species richness 0 (the RRFSp can't be calculated)
bir<-subset(taxa,RRFSp!="NA")
nrow(bir)#44
summary(bir)
#Absolute Response Ratio
bir$RRFSp2<-abs(bir$RRFSp)
#Remove outlyers
outlyersbir<-boxplot(bir$RRFSp2)
sort(outlyersbir$out) # >2.1
bir2<-subset(bir, bir$RRFSp2<2.1)
boxplot(bir2$RRFSp2)
bir2$RRFSp2<-bir2$RRFSp2+0.5

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
				"PatchArea25","PatchArea10",
				"PatchPtoA10","PatchPtoA75","PatchPtoA25",
				"nPatches5","nPatches75","nPatches50",
				"p.Forest5km","p.Forest10km", "p.Forest25km",
				"Continent_Island","Zoogeographical_Realm_Holt_etal_2013",
				"Mean_Elevation","Mean_Rainfall")], 1, anyNA),]
nrow(bir4)#n=33
## Good of fitness and normality
global.model<-glm(RRFSp2 ~1+SuccStage+Mean.Age+PreviousLandUse+
				PatchArea25+PatchArea10+
				PatchPtoA10+PatchPtoA75+PatchPtoA25+
				nPatches5+nPatches75+nPatches50+
				p.Forest5km+p.Forest10km+p.Forest25km+
				Continent_Island+Zoogeographical_Realm_Holt_etal_2013+
				Mean_Elevation+Mean_Rainfall,data=bir4, Gamma("identity"))
plot(global.model)
hist(resid(global.model), breaks=20)
#1- Model selection
model <- glmulti(RRFSp2 ~1+SuccStage+Mean.Age+PreviousLandUse+
				PatchArea25+PatchArea10+
				PatchPtoA10+PatchPtoA75+PatchPtoA25+
				nPatches5+nPatches75+nPatches50+
				p.Forest5km+p.Forest10km+p.Forest25km+
				Continent_Island+Zoogeographical_Realm_Holt_etal_2013+
				Mean_Elevation+Mean_Rainfall,data=bir4,
               level=1, family=Gamma (link="identity"), crit="aicc")
ranks<-weightable(model)
ranks
plot(model,type="s")
ranks7<-ranks[ranks$aicc <= min(ranks$aicc) + 7,]
ranks7
#100 models + null

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
m1<-glm(RRFSp2~1, data=bir5, Gamma('identity'))
m2<-glm(RRFSp2~1+Mean.Age+PatchPtoA25+Mean_Elevation, data=bir5, Gamma('identity'))
m3<-glm(RRFSp2~1+Mean.Age+PatchArea25+Mean_Elevation, data=bir5, Gamma('identity'))
m4<-glm(RRFSp2~1+Mean.Age+Mean_Elevation, data=bir5, Gamma('identity'))
m5<-glm(RRFSp2~1+Mean.Age+PatchPtoA10+Mean_Elevation, data=bir5, Gamma('identity'))
m6<-glm(RRFSp2~1+Mean.Age, data=bir5, Gamma('identity'))
m7<-glm(RRFSp2~1+Mean.Age+PatchPtoA10, data=bir5, Gamma('identity'))
m8<-glm(RRFSp2~1+Mean.Age+PatchPtoA25, data=bir5, Gamma('identity'))
m9<-glm(RRFSp2~1+Mean.Age+PatchArea25+PatchPtoA75+Mean_Elevation, data=bir5, Gamma('identity'))
m10<-glm(RRFSp2~1+Mean.Age+PatchArea10+Mean_Elevation, data=bir5, Gamma('identity'))
m11<-glm(RRFSp2~1+Mean.Age+nPatches50+Mean_Elevation, data=bir5, Gamma('identity'))
m12<-glm(RRFSp2~1+Mean.Age+PatchArea25+PatchArea10+Mean_Elevation, data=bir5, Gamma('identity'))
m13<-glm(RRFSp2~1+Mean.Age+nPatches75+Mean_Elevation, data=bir5, Gamma('identity'))
m14<-glm(RRFSp2~1+Mean.Age+nPatches50, data=bir5, Gamma('identity'))
m15<-glm(RRFSp2~1+Mean.Age+PatchArea25, data=bir5, Gamma('identity'))
m16<-glm(RRFSp2~1+Mean.Age+nPatches75, data=bir5, Gamma('identity'))
m17<-glm(RRFSp2~1+Mean.Age+PatchPtoA75+PatchPtoA25+Mean_Elevation, data=bir5, Gamma('identity'))
m18<-glm(RRFSp2~1+Mean.Age+PatchArea25+p.Forest25km+Mean_Elevation, data=bir5, Gamma('identity'))
m19<-glm(RRFSp2~1+Mean.Age+PatchPtoA10+PatchPtoA75, data=bir5, Gamma('identity'))
m20<-glm(RRFSp2~1+Mean.Age+PatchPtoA75, data=bir5, Gamma('identity'))
m21<-glm(RRFSp2~1+Mean.Age+PatchArea25+nPatches5+Mean_Elevation, data=bir5, Gamma('identity'))
m22<-glm(RRFSp2~1+Mean.Age+PatchPtoA25+Mean_Elevation+Mean_Rainfall, data=bir5, Gamma('identity'))
m23<-glm(RRFSp2~1+SuccStage, data=bir5, Gamma('identity'))
m24<-glm(RRFSp2~1+Mean.Age+PatchPtoA75+Mean_Elevation, data=bir5, Gamma('identity'))
m25<-glm(RRFSp2~1+Mean.Age+PatchPtoA25+nPatches5+Mean_Elevation, data=bir5, Gamma('identity'))
m26<-glm(RRFSp2~1+Mean.Age+PatchArea25+p.Forest10km+Mean_Elevation, data=bir5, Gamma('identity'))
m27<-glm(RRFSp2~1+Mean.Age+PatchArea10, data=bir5, Gamma('identity'))
m28<-glm(RRFSp2~1+Mean.Age+PatchArea10+PatchPtoA25+Mean_Elevation, data=bir5, Gamma('identity'))
m29<-glm(RRFSp2~1+Mean.Age+PatchArea25+Mean_Elevation+Mean_Rainfall, data=bir5, Gamma('identity'))
m30<-glm(RRFSp2~1+Mean.Age+PatchPtoA10+PatchPtoA25+Mean_Elevation, data=bir5, Gamma('identity'))
m31<-glm(RRFSp2~1+Mean.Age+PatchArea25+PatchPtoA25+Mean_Elevation, data=bir5, Gamma('identity'))
m32<-glm(RRFSp2~1+Mean.Age+PatchArea25+nPatches75+Mean_Elevation, data=bir5, Gamma('identity'))
m33<-glm(RRFSp2~1+Mean.Age+PatchPtoA75+PatchPtoA25, data=bir5, Gamma('identity'))
m34<-glm(RRFSp2~1+Mean.Age+PatchArea25+nPatches50+Mean_Elevation, data=bir5, Gamma('identity'))
m35<-glm(RRFSp2~1+Mean.Age+PatchPtoA25+nPatches50+Mean_Elevation, data=bir5, Gamma('identity'))
m36<-glm(RRFSp2~1+Continent_Island+Mean.Age+PatchPtoA25+Mean_Elevation, data=bir5, Gamma('identity'))
m37<-glm(RRFSp2~1+Mean.Age+PatchPtoA25+p.Forest25km+Mean_Elevation, data=bir5, Gamma('identity'))
m38<-glm(RRFSp2~1+Continent_Island+Mean.Age+PatchArea25+Mean_Elevation, data=bir5, Gamma('identity'))
m39<-glm(RRFSp2~1+Mean.Age+PatchPtoA25+nPatches75+Mean_Elevation, data=bir5, Gamma('identity'))
m40<-glm(RRFSp2~1+Mean.Age+PatchPtoA25+p.Forest5km+Mean_Elevation, data=bir5, Gamma('identity'))
m41<-glm(RRFSp2~1+Mean.Age+PatchPtoA25+p.Forest10km+Mean_Elevation, data=bir5, Gamma('identity'))
m42<-glm(RRFSp2~1+Mean.Age+PatchArea25+p.Forest5km+Mean_Elevation, data=bir5, Gamma('identity'))
m43<-glm(RRFSp2~1+Mean.Age+PatchArea25+PatchPtoA10+Mean_Elevation, data=bir5, Gamma('identity'))
m44<-glm(RRFSp2~1+Mean.Age+PatchPtoA25+p.Forest25km, data=bir5, Gamma('identity'))
m45<-glm(RRFSp2~1+Mean.Age+PatchPtoA10+p.Forest25km, data=bir5, Gamma('identity'))
m46<-glm(RRFSp2~1+Mean.Age+PatchPtoA10+PatchPtoA75+Mean_Elevation, data=bir5, Gamma('identity'))
m47<-glm(RRFSp2~1+Mean.Age+PatchArea10+PatchPtoA10, data=bir5, Gamma('identity'))
m48<-glm(RRFSp2~1+Continent_Island+Mean.Age+Mean_Elevation, data=bir5, Gamma('identity'))
m49<-glm(RRFSp2~1+Mean.Age+p.Forest10km+Mean_Elevation, data=bir5, Gamma('identity'))
m50<-glm(RRFSp2~1+Mean.Age+PatchArea25+PatchPtoA75, data=bir5, Gamma('identity'))
m51<-glm(RRFSp2~1+Mean.Age+PatchArea10+PatchPtoA75+Mean_Elevation, data=bir5, Gamma('identity'))
m52<-glm(RRFSp2~1+Mean.Age+p.Forest5km+Mean_Elevation, data=bir5, Gamma('identity'))
m53<-glm(RRFSp2~1+Mean.Age+PatchArea25+p.Forest25km, data=bir5, Gamma('identity'))
m54<-glm(RRFSp2~1+Mean.Age+nPatches50+p.Forest25km, data=bir5, Gamma('identity'))
m55<-glm(RRFSp2~1+Mean.Age+PatchPtoA25+Mean_Rainfall, data=bir5, Gamma('identity'))
m56<-glm(RRFSp2~1+Mean.Age+PatchPtoA10+nPatches5+Mean_Elevation, data=bir5, Gamma('identity'))
m57<-glm(RRFSp2~1+Mean.Age+Mean_Elevation+Mean_Rainfall, data=bir5, Gamma('identity'))
m58<-glm(RRFSp2~1+Mean.Age+PatchArea25+PatchPtoA25, data=bir5, Gamma('identity'))
m59<-glm(RRFSp2~1+Mean.Age+p.Forest25km+Mean_Elevation, data=bir5, Gamma('identity'))
m60<-glm(RRFSp2~1+Mean.Age+p.Forest10km, data=bir5, Gamma('identity'))
m61<-glm(RRFSp2~1+Mean.Age+nPatches5+Mean_Elevation, data=bir5, Gamma('identity'))
m62<-glm(RRFSp2~1+Mean.Age+PatchArea10+PatchPtoA75, data=bir5, Gamma('identity'))
m63<-glm(RRFSp2~1+Continent_Island+Mean.Age, data=bir5, Gamma('identity'))
m64<-glm(RRFSp2~1+SuccStage+Mean_Elevation, data=bir5, Gamma('identity'))
m65<-glm(RRFSp2~1+Mean.Age+nPatches5, data=bir5, Gamma('identity'))
m66<-glm(RRFSp2~1+Mean.Age+PatchPtoA10+nPatches50+Mean_Elevation, data=bir5, Gamma('identity'))
m67<-glm(RRFSp2~1+Mean.Age+p.Forest25km, data=bir5, Gamma('identity'))
m68<-glm(RRFSp2~1+Mean.Age+p.Forest5km, data=bir5, Gamma('identity'))
m69<-glm(RRFSp2~1+Mean.Age+Mean_Rainfall, data=bir5, Gamma('identity'))
m70<-glm(RRFSp2~1+Mean.Age+PatchPtoA10+nPatches75+Mean_Elevation, data=bir5, Gamma('identity'))
m71<-glm(RRFSp2~1+Mean.Age+PatchPtoA10+Mean_Rainfall, data=bir5, Gamma('identity'))
m72<-glm(RRFSp2~1+Mean.Age+PatchArea25+PatchArea10+PatchPtoA75+Mean_Elevation, data=bir5, Gamma('identity'))
m73<-glm(RRFSp2~1+Continent_Island+Mean.Age+PatchPtoA10+Mean_Elevation, data=bir5, Gamma('identity'))
m74<-glm(RRFSp2~1+Mean.Age+PatchArea10+PatchPtoA25, data=bir5, Gamma('identity'))
m75<-glm(RRFSp2~1+Mean.Age+PatchPtoA10+nPatches75, data=bir5, Gamma('identity'))
m76<-glm(RRFSp2~1+Mean.Age+PatchPtoA10+nPatches50, data=bir5, Gamma('identity'))
m77<-glm(RRFSp2~1+Mean.Age+PatchPtoA25+nPatches50, data=bir5, Gamma('identity'))
m78<-glm(RRFSp2~1+SuccStage+Mean.Age, data=bir5, Gamma('identity'))
m79<-glm(RRFSp2~1+Mean.Age+nPatches5+nPatches50+Mean_Elevation, data=bir5, Gamma('identity'))
m80<-glm(RRFSp2~1+Mean.Age+PatchPtoA10+p.Forest10km, data=bir5, Gamma('identity'))
m81<-glm(RRFSp2~1+Mean.Age+PatchPtoA10+Mean_Elevation+Mean_Rainfall, data=bir5, Gamma('identity'))
m82<-glm(RRFSp2~1+Mean.Age+PatchPtoA75+nPatches50, data=bir5, Gamma('identity'))
m83<-glm(RRFSp2~1+Mean.Age+PatchPtoA25+nPatches75, data=bir5, Gamma('identity'))
m84<-glm(RRFSp2~1+Mean.Age+PatchPtoA10+p.Forest5km+Mean_Elevation, data=bir5, Gamma('identity'))
m85<-glm(RRFSp2~1+Mean.Age+PatchPtoA10+p.Forest25km+Mean_Elevation, data=bir5, Gamma('identity'))
m86<-glm(RRFSp2~1+Mean.Age+PatchArea10+PatchPtoA10+Mean_Elevation, data=bir5, Gamma('identity'))
m87<-glm(RRFSp2~1+Mean.Age+PatchPtoA10+p.Forest10km+Mean_Elevation, data=bir5, Gamma('identity'))
m88<-glm(RRFSp2~1+Mean.Age+PatchPtoA25+p.Forest10km, data=bir5, Gamma('identity'))
m89<-glm(RRFSp2~1+Mean.Age+PatchArea10+nPatches75+Mean_Elevation, data=bir5, Gamma('identity'))
m90<-glm(RRFSp2~1+Mean.Age+PatchPtoA10+PatchPtoA25, data=bir5, Gamma('identity'))
m91<-glm(RRFSp2~1+Mean.Age+PatchPtoA10+nPatches5, data=bir5, Gamma('identity'))
m92<-glm(RRFSp2~1+Mean.Age+PatchPtoA10+p.Forest5km, data=bir5, Gamma('identity'))
m93<-glm(RRFSp2~1+Mean.Age+PatchArea25+nPatches50+p.Forest25km, data=bir5, Gamma('identity'))
m94<-glm(RRFSp2~1+Mean.Age+PatchArea10+nPatches50+Mean_Elevation, data=bir5, Gamma('identity'))
m95<-glm(RRFSp2~1+Mean.Age+PatchArea10+nPatches5+Mean_Elevation, data=bir5, Gamma('identity'))
m96<-glm(RRFSp2~1+Mean.Age+PatchArea10+p.Forest25km, data=bir5, Gamma('identity'))
m97<-glm(RRFSp2~1+Mean.Age+PatchPtoA25+p.Forest5km, data=bir5, Gamma('identity'))
m98<-glm(RRFSp2~1+Continent_Island+Mean.Age+nPatches50+Mean_Elevation, data=bir5, Gamma('identity'))
m99<-glm(RRFSp2~1+Continent_Island+Mean.Age+PatchPtoA10, data=bir5, Gamma('identity'))
m100<-glm(RRFSp2~1+Mean.Age+PatchPtoA75+nPatches75, data=bir5, Gamma('identity'))
m101<-glm(RRFSp2~1+Continent_Island+Mean.Age+nPatches75+Mean_Elevation, data=bir5, Gamma('identity'))


ranks<-m.AICc(list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,
					m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39,m40,
					m41,m42,m43,m44,m45,m46,m47,m48,m49,m50,m51,m52,m53,m54,m55,m56,m57,m58,m59,m60,
					m61,m62,m63,m64,m65,m66,m67,m68,m69,m70,m71,m72,m73,m74,m75,m76,m77,m78,m79,m80,
					m81,m82,m83,m84,m85,m86,m87,m88,m89,m90,m91,m92,m93,m94,m95,m96,m97,m98,m99,m100,m101),
					as.integer(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20",
					"21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40",
					"41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60",
					"61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80",
					"81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101")),
					length(bir5$RRFSp2))

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
#m3 (RRFSp2 ~ 1 + Mean.Age)
m2b<-lm(RRFSp2 ~ 1 + Mean.Age, data=bir4)
summary(m2b)
Am2b<-anova(m2b)
Am2bss<-Am2b$"Sum Sq"
print(cbind(Am2b,PctExp=Am2bss/sum(Am2bss)*100)) #PctExp: % explain variation

#~~~~~~~~~~~~~~
#Mammals
taxa<-subset(mydata,Taxa=="Mammals")
#Exclude the with species richness 0 (the RRFSp can't be calculated)
mam<-subset(taxa,RRFSp!="NA")
nrow(mam)#31
summary(mam)
#Absolute Response Ratio
mam$RRFSp2<-abs(mam$RRFSp)
#Remove outlyers
outlyersmam<-boxplot(mam$RRFSp2)
sort(outlyersmam$out) # >1.6
mam2<-subset(mam, mam$RRFSp2<1.6)
boxplot(mam2$RRFSp2)
mam2$RRFSp2<-mam2$RRFSp2+0.5

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
				"nPatches5","nPatches25","nPatches50",
				"PatchArea10","PatchArea75",
				"p.Forest10km",
				"PatchPtoA5","PatchPtoA10","PatchPtoA25","PatchPtoA100","PatchPtoA50",
				"Continent_Island","Zoogeographical_Realm_Holt_etal_2013",
				"Mean_Elevation","Mean_Rainfall")], 1, anyNA),]
nrow(mam4)#n=22
## Good of fitness and normality
global.model<-glm(RRFSp2 ~1+SuccStage+Mean.Age+PreviousLandUse+
				nPatches5+nPatches25+nPatches50+
				PatchArea10+PatchArea75+
				p.Forest10km+
				PatchPtoA5+PatchPtoA10+PatchPtoA25+PatchPtoA100+PatchPtoA50+
				Continent_Island+Zoogeographical_Realm_Holt_etal_2013+
				Mean_Elevation+Mean_Rainfall,data=mam4, Gamma("identity"))
plot(global.model)
hist(resid(global.model), breaks=20)
#1- Model selection
model <- glmulti(RRFSp2 ~1+SuccStage+Mean.Age+PreviousLandUse+
				nPatches5+nPatches25+nPatches50+
				PatchArea10+PatchArea75+
				p.Forest10km+
				PatchPtoA5+PatchPtoA10+PatchPtoA25+PatchPtoA100+PatchPtoA50+
				Continent_Island+Zoogeographical_Realm_Holt_etal_2013+
				Mean_Elevation+Mean_Rainfall,data=mam4,
               level=1, family=Gamma (link="identity"), crit="aicc")
ranks<-weightable(model)
ranks
plot(model,type="s")
ranks7<-ranks[ranks$aicc <= min(ranks$aicc) + 7,]
ranks7
#100 models + null

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
m1<-glm(RRFSp2~1, data=mam5, Gamma('identity'))
m2<-glm(RRFSp2~1+Mean.Age+p.Forest10km+PatchPtoA10+PatchPtoA50, data=mam5, Gamma('identity'))
m3<-glm(RRFSp2~1+p.Forest10km+PatchPtoA10+PatchPtoA50, data=mam5, Gamma('identity'))
m4<-glm(RRFSp2~1+Mean.Age+p.Forest10km+PatchPtoA10+PatchPtoA25, data=mam5, Gamma('identity'))
m5<-glm(RRFSp2~1+p.Forest10km+PatchPtoA10+PatchPtoA25, data=mam5, Gamma('identity'))
m6<-glm(RRFSp2~1+Continent_Island+p.Forest10km+PatchPtoA10+PatchPtoA50, data=mam5, Gamma('identity'))
m7<-glm(RRFSp2~1+Mean_Rainfall, data=mam5, Gamma('identity'))
m8<-glm(RRFSp2~1+PatchArea10+p.Forest10km+PatchPtoA10+PatchPtoA25, data=mam5, Gamma('identity'))
m9<-glm(RRFSp2~1+Mean.Age+p.Forest10km+PatchPtoA10+PatchPtoA100+PatchPtoA50, data=mam5, Gamma('identity'))
m10<-glm(RRFSp2~1+Continent_Island+PatchArea10+PatchPtoA5+PatchPtoA25+PatchPtoA50, data=mam5, Gamma('identity'))
m11<-glm(RRFSp2~1+Continent_Island+Mean.Age+p.Forest10km+PatchPtoA10+PatchPtoA50, data=mam5, Gamma('identity'))
m12<-glm(RRFSp2~1+Mean.Age+Mean_Rainfall, data=mam5, Gamma('identity'))
m13<-glm(RRFSp2~1+p.Forest10km+PatchPtoA10+PatchPtoA100+PatchPtoA50, data=mam5, Gamma('identity'))
m14<-glm(RRFSp2~1+Mean.Age+PatchArea75+p.Forest10km+PatchPtoA10+PatchPtoA50, data=mam5, Gamma('identity'))
m15<-glm(RRFSp2~1+p.Forest10km+PatchPtoA10+PatchPtoA25+PatchPtoA50, data=mam5, Gamma('identity'))
m16<-glm(RRFSp2~1+Continent_Island+PatchPtoA5+PatchPtoA10+PatchPtoA25+PatchPtoA50, data=mam5, Gamma('identity'))
m17<-glm(RRFSp2~1+nPatches25+p.Forest10km+PatchPtoA10+PatchPtoA50, data=mam5, Gamma('identity'))
m18<-glm(RRFSp2~1+Mean.Age+p.Forest10km+PatchPtoA10+PatchPtoA25+PatchPtoA50, data=mam5, Gamma('identity'))
m19<-glm(RRFSp2~1+p.Forest10km+PatchPtoA5+PatchPtoA10+PatchPtoA50, data=mam5, Gamma('identity'))
m20<-glm(RRFSp2~1+PatchArea75+p.Forest10km+PatchPtoA10+PatchPtoA50, data=mam5, Gamma('identity'))
m21<-glm(RRFSp2~1+Continent_Island+Mean_Rainfall, data=mam5, Gamma('identity'))
m22<-glm(RRFSp2~1+PatchArea10+p.Forest10km+PatchPtoA10+PatchPtoA50, data=mam5, Gamma('identity'))
m23<-glm(RRFSp2~1+nPatches50+p.Forest10km+PatchPtoA10+PatchPtoA50, data=mam5, Gamma('identity'))
m24<-glm(RRFSp2~1+p.Forest10km+PatchPtoA10+PatchPtoA50+Mean_Elevation, data=mam5, Gamma('identity'))
m25<-glm(RRFSp2~1+Mean.Age+nPatches5+p.Forest10km+PatchPtoA10+PatchPtoA50, data=mam5, Gamma('identity'))
m26<-glm(RRFSp2~1+Continent_Island+p.Forest10km+PatchPtoA10+PatchPtoA100, data=mam5, Gamma('identity'))
m27<-glm(RRFSp2~1+nPatches5+p.Forest10km+PatchPtoA10+PatchPtoA50, data=mam5, Gamma('identity'))
m28<-glm(RRFSp2~1+p.Forest10km+PatchPtoA10+PatchPtoA50+Mean_Rainfall, data=mam5, Gamma('identity'))
m29<-glm(RRFSp2~1+PatchPtoA10+Mean_Rainfall, data=mam5, Gamma('identity'))
m30<-glm(RRFSp2~1+Mean.Age+PatchArea10+p.Forest10km+PatchPtoA10+PatchPtoA25, data=mam5, Gamma('identity'))
m31<-glm(RRFSp2~1+p.Forest10km+PatchPtoA10+PatchPtoA25+PatchPtoA100, data=mam5, Gamma('identity'))
m32<-glm(RRFSp2~1+Mean.Age+PatchArea10+p.Forest10km+PatchPtoA10+PatchPtoA50, data=mam5, Gamma('identity'))
m33<-glm(RRFSp2~1+PatchPtoA5+PatchPtoA10+PatchPtoA50, data=mam5, Gamma('identity'))
m34<-glm(RRFSp2~1+PatchArea75+p.Forest10km+PatchPtoA10+PatchPtoA25, data=mam5, Gamma('identity'))
m35<-glm(RRFSp2~1+Mean.Age+p.Forest10km+PatchPtoA10+PatchPtoA50+Mean_Rainfall, data=mam5, Gamma('identity'))
m36<-glm(RRFSp2~1+p.Forest10km+PatchPtoA10+PatchPtoA100, data=mam5, Gamma('identity'))
m37<-glm(RRFSp2~1+Mean.Age+p.Forest10km+PatchPtoA5+PatchPtoA10+PatchPtoA50, data=mam5, Gamma('identity'))
m38<-glm(RRFSp2~1+Mean.Age+nPatches25+p.Forest10km+PatchPtoA10+PatchPtoA50, data=mam5, Gamma('identity'))
m39<-glm(RRFSp2~1+Mean.Age+nPatches50+p.Forest10km+PatchPtoA10+PatchPtoA50, data=mam5, Gamma('identity'))
m40<-glm(RRFSp2~1+PatchArea75+p.Forest10km+PatchPtoA10, data=mam5, Gamma('identity'))
m41<-glm(RRFSp2~1+Mean.Age+p.Forest10km+PatchPtoA10+PatchPtoA50+Mean_Elevation, data=mam5, Gamma('identity'))
m42<-glm(RRFSp2~1+p.Forest10km+PatchPtoA10+PatchPtoA25+Mean_Rainfall, data=mam5, Gamma('identity'))
m43<-glm(RRFSp2~1+PatchPtoA5+Mean_Rainfall, data=mam5, Gamma('identity'))
m44<-glm(RRFSp2~1+Continent_Island+PatchArea10+p.Forest10km+PatchPtoA10+PatchPtoA50, data=mam5, Gamma('identity'))
m45<-glm(RRFSp2~1+PatchPtoA100+Mean_Rainfall, data=mam5, Gamma('identity'))
m46<-glm(RRFSp2~1+Mean.Age+PatchPtoA5+PatchPtoA10+PatchPtoA50, data=mam5, Gamma('identity'))
m47<-glm(RRFSp2~1+Mean_Elevation+Mean_Rainfall, data=mam5, Gamma('identity'))
m48<-glm(RRFSp2~1+PatchPtoA50+Mean_Rainfall, data=mam5, Gamma('identity'))
m49<-glm(RRFSp2~1+nPatches5+Mean_Rainfall, data=mam5, Gamma('identity'))
m50<-glm(RRFSp2~1+PatchArea75+Mean_Rainfall, data=mam5, Gamma('identity'))
m51<-glm(RRFSp2~1+PatchPtoA25+Mean_Rainfall, data=mam5, Gamma('identity'))
m52<-glm(RRFSp2~1+Mean.Age+p.Forest10km+PatchPtoA10+PatchPtoA25+PatchPtoA100, data=mam5, Gamma('identity'))
m53<-glm(RRFSp2~1+nPatches50+Mean_Rainfall, data=mam5, Gamma('identity'))
m54<-glm(RRFSp2~1+p.Forest10km+Mean_Rainfall, data=mam5, Gamma('identity'))
m55<-glm(RRFSp2~1+nPatches25+Mean_Rainfall, data=mam5, Gamma('identity'))
m56<-glm(RRFSp2~1+PatchArea10+Mean_Rainfall, data=mam5, Gamma('identity'))
m57<-glm(RRFSp2~1+PatchArea10+PatchPtoA10+Mean_Rainfall, data=mam5, Gamma('identity'))
m58<-glm(RRFSp2~1+Continent_Island+PatchArea10+p.Forest10km+PatchPtoA10+PatchPtoA25, data=mam5, Gamma('identity'))
m59<-glm(RRFSp2~1+Mean.Age+p.Forest10km+PatchPtoA5+PatchPtoA10+PatchPtoA25, data=mam5, Gamma('identity'))
m60<-glm(RRFSp2~1+Continent_Island+p.Forest10km+PatchPtoA10+PatchPtoA25, data=mam5, Gamma('identity'))
m61<-glm(RRFSp2~1+PatchArea10+PatchArea75+PatchPtoA10+PatchPtoA25, data=mam5, Gamma('identity'))
m62<-glm(RRFSp2~1+Mean.Age+nPatches5+p.Forest10km+PatchPtoA10+PatchPtoA25, data=mam5, Gamma('identity'))
m63<-glm(RRFSp2~1+p.Forest10km+PatchPtoA5+PatchPtoA10+PatchPtoA25, data=mam5, Gamma('identity'))
m64<-glm(RRFSp2~1+Mean.Age+PatchArea75+p.Forest10km+PatchPtoA10+PatchPtoA25, data=mam5, Gamma('identity'))
m65<-glm(RRFSp2~1+Mean.Age+PatchPtoA10+PatchPtoA100+PatchPtoA50, data=mam5, Gamma('identity'))
m66<-glm(RRFSp2~1+Mean.Age+p.Forest10km+PatchPtoA10+PatchPtoA25+Mean_Elevation, data=mam5, Gamma('identity'))
m67<-glm(RRFSp2~1+nPatches5+p.Forest10km+PatchPtoA10+PatchPtoA25, data=mam5, Gamma('identity'))
m68<-glm(RRFSp2~1+Mean.Age+p.Forest10km+PatchPtoA10+PatchPtoA100, data=mam5, Gamma('identity'))
m69<-glm(RRFSp2~1+Continent_Island+Mean.Age+Mean_Rainfall, data=mam5, Gamma('identity'))
m70<-glm(RRFSp2~1+Continent_Island+p.Forest10km+PatchPtoA5+PatchPtoA10+PatchPtoA50, data=mam5, Gamma('identity'))
m71<-glm(RRFSp2~1+Continent_Island+PatchArea75+p.Forest10km+PatchPtoA10, data=mam5, Gamma('identity'))
m72<-glm(RRFSp2~1+Continent_Island+nPatches25+p.Forest10km+PatchPtoA10+PatchPtoA50, data=mam5, Gamma('identity'))
m73<-glm(RRFSp2~1+nPatches25+p.Forest10km+PatchPtoA10+PatchPtoA25, data=mam5, Gamma('identity'))
m74<-glm(RRFSp2~1+Mean.Age+p.Forest10km+PatchPtoA10+PatchPtoA25+Mean_Rainfall, data=mam5, Gamma('identity'))
m75<-glm(RRFSp2~1+nPatches50+p.Forest10km+PatchPtoA10+PatchPtoA25, data=mam5, Gamma('identity'))
m76<-glm(RRFSp2~1+p.Forest10km+PatchPtoA10+PatchPtoA25+Mean_Elevation, data=mam5, Gamma('identity'))
m77<-glm(RRFSp2~1+Mean.Age+PatchPtoA10+PatchPtoA100+PatchPtoA50+Mean_Rainfall, data=mam5, Gamma('identity'))
m78<-glm(RRFSp2~1+Mean.Age+PatchPtoA10+Mean_Rainfall, data=mam5, Gamma('identity'))
m79<-glm(RRFSp2~1+PatchPtoA10+PatchPtoA25+Mean_Rainfall, data=mam5, Gamma('identity'))
m80<-glm(RRFSp2~1+Continent_Island+PatchPtoA5+PatchPtoA25+PatchPtoA50, data=mam5, Gamma('identity'))
m81<-glm(RRFSp2~1+PatchPtoA5+PatchPtoA10+PatchPtoA50+Mean_Rainfall, data=mam5, Gamma('identity'))
m82<-glm(RRFSp2~1+Mean.Age+PatchPtoA10+PatchPtoA25+Mean_Rainfall, data=mam5, Gamma('identity'))
m83<-glm(RRFSp2~1+nPatches5+PatchArea10+p.Forest10km+PatchPtoA10+PatchPtoA25, data=mam5, Gamma('identity'))
m84<-glm(RRFSp2~1+Mean.Age+nPatches50+p.Forest10km+PatchPtoA10+PatchPtoA25, data=mam5, Gamma('identity'))
m85<-glm(RRFSp2~1+Continent_Island+Mean.Age+p.Forest10km+PatchPtoA10+PatchPtoA25, data=mam5, Gamma('identity'))
m86<-glm(RRFSp2~1+Mean.Age+nPatches25+p.Forest10km+PatchPtoA10+PatchPtoA25, data=mam5, Gamma('identity'))
m87<-glm(RRFSp2~1+Continent_Island+p.Forest10km+PatchPtoA10+PatchPtoA100+PatchPtoA50, data=mam5, Gamma('identity'))
m88<-glm(RRFSp2~1+Mean.Age+nPatches5+p.Forest10km+PatchPtoA10+PatchPtoA100+PatchPtoA50, data=mam5, Gamma('identity'))
m89<-glm(RRFSp2~1+PatchPtoA10+PatchPtoA50+Mean_Rainfall, data=mam5, Gamma('identity'))
m90<-glm(RRFSp2~1+Mean.Age+nPatches25+Mean_Rainfall, data=mam5, Gamma('identity'))
m91<-glm(RRFSp2~1+Mean.Age+PatchArea75+PatchPtoA10+PatchPtoA50, data=mam5, Gamma('identity'))
m92<-glm(RRFSp2~1+Mean.Age+p.Forest10km+Mean_Rainfall, data=mam5, Gamma('identity'))
m93<-glm(RRFSp2~1+Continent_Island+PatchArea10+PatchPtoA10+Mean_Rainfall, data=mam5, Gamma('identity'))
m94<-glm(RRFSp2~1+Continent_Island+nPatches50+p.Forest10km+PatchPtoA10+PatchPtoA50, data=mam5, Gamma('identity'))
m95<-glm(RRFSp2~1+Mean.Age+nPatches50+Mean_Rainfall, data=mam5, Gamma('identity'))
m96<-glm(RRFSp2~1+Continent_Island+PatchArea75+p.Forest10km+PatchPtoA10+PatchPtoA50, data=mam5, Gamma('identity'))
m97<-glm(RRFSp2~1+PatchArea10+PatchArea75+p.Forest10km+PatchPtoA10+PatchPtoA25, data=mam5, Gamma('identity'))
m98<-glm(RRFSp2~1+Mean.Age+PatchArea10+Mean_Rainfall, data=mam5, Gamma('identity'))
m99<-glm(RRFSp2~1+Continent_Island+p.Forest10km+PatchPtoA10+PatchPtoA25+PatchPtoA50, data=mam5, Gamma('identity'))
m100<-glm(RRFSp2~1+Continent_Island+p.Forest10km+PatchPtoA10+PatchPtoA50+Mean_Rainfall, data=mam5, Gamma('identity'))
m101<-glm(RRFSp2~1+Mean.Age+PatchPtoA100+Mean_Rainfall, data=mam5, Gamma('identity'))


ranks<-m.AICc(list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,
					m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39,m40,
					m41,m42,m43,m44,m45,m46,m47,m48,m49,m50,m51,m52,m53,m54,m55,m56,m57,m58,m59,m60,
					m61,m62,m63,m64,m65,m66,m67,m68,m69,m70,m71,m72,m73,m74,m75,m76,m77,m78,m79,m80,
					m81,m82,m83,m84,m85,m86,m87,m88,m89,m90,m91,m92,m93,m94,m95,m96,m97,m98,m99,m100,m101),
					as.integer(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20",
					"21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40",
					"41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60",
					"61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80",
					"81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101")),
					length(mam5$RRFSp2))

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
#m4 (RRFSp2~1+Mean_Rainfall)
m4b<-lm(RRFSp2~1 + Mean_Rainfall, data=mam4)
summary(m4b)
Am4b<-anova(m4b)
Am4bss<-Am4b$"Sum Sq"
print(cbind(Am4b,PctExp=Am4bss/sum(Am4bss)*100)) #PctExp: 39.27% explain variation

##Extra fig
library(ggplot2)
library(ggrepel)
ggplot(mam4,aes(x=Mean_Rainfall, y=RRFSp2))+
	geom_point(size=5)+
	geom_smooth(method="lm",se=F)+
  geom_text_repel(aes(label=mam4$Name), size=6, segment.color='#cccccc', arrow = arrow(length=unit(0.01,'npc')),force=10)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = 	element_blank())

