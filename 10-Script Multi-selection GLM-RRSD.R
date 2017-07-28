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
#Exclude the with species richness 0 (the RRSD can't be calculated)
bir<-subset(taxa,RRSD!="NA")
nrow(bir)#43
summary(bir)
#Absolute Response Ratio
bir$RRSD2<-abs(bir$RRSD)
#Remove outlyers
outlyersbir<-boxplot(bir$RRSD2)
sort(outlyersbir$out) # >1.2
bir2<-subset(bir, bir$RRSD2<1.2)
boxplot(bir2$RRSD2)
bir2$RRSD2<-bir2$RRSD2+0.5

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
				"PatchArea25","PatchArea10","PatchArea50",
				"p.Forest5km","p.Forest25km",
				"nPatches5","nPatches50",
				"Continent_Island","Zoogeographical_Realm_Holt_etal_2013",
				"Mean_Elevation","Mean_Rainfall")], 1, anyNA),]
nrow(bir4)#n=33
## Good of fitness and normality
global.model<-glm(RRSD2 ~1+SuccStage+Mean.Age+PreviousLandUse+
				PatchPtoA75+PatchPtoA25+PatchPtoA10+
				PatchArea25+PatchArea10+PatchArea50+
				p.Forest5km+p.Forest25km+
				nPatches5+nPatches50+
				Continent_Island+Zoogeographical_Realm_Holt_etal_2013+
				Mean_Elevation+Mean_Rainfall,data=bir4, Gamma("identity"))
plot(global.model)
hist(resid(global.model), breaks=20)
#1- Model selection
model <- glmulti(RRSD2 ~1+SuccStage+Mean.Age+PreviousLandUse+
				PatchPtoA75+PatchPtoA25+PatchPtoA10+
				PatchArea25+PatchArea10+PatchArea50+
				p.Forest5km+p.Forest25km+
				nPatches5+nPatches50+
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
m1<-glm(RRSD2~1, data=bir5, Gamma('identity'))
m2<-glm(RRSD2~1+Mean.Age+PatchArea10+PatchArea50+Mean_Elevation, data=bir5, Gamma('identity'))
m3<-glm(RRSD2~1+Continent_Island+Mean.Age+PatchArea10+PatchArea50+Mean_Elevation, data=bir5, Gamma('identity'))
m4<-glm(RRSD2~1+Mean.Age+PatchPtoA10+PatchArea10+PatchArea50+Mean_Elevation, data=bir5, Gamma('identity'))
m5<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA10+PatchArea10+PatchArea50+Mean_Elevation, data=bir5, Gamma('identity'))
m6<-glm(RRSD2~1+Mean.Age+PatchArea10+PatchArea50+Mean_Elevation+Mean_Rainfall, data=bir5, Gamma('identity'))
m7<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea25+PatchArea10, data=bir5, Gamma('identity'))
m8<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea10+Mean_Elevation, data=bir5, Gamma('identity'))
m9<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchArea10+PatchArea50+Mean_Elevation, data=bir5, Gamma('identity'))
m10<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea25+PatchArea10+Mean_Elevation, data=bir5, Gamma('identity'))
m11<-glm(RRSD2~1+Continent_Island+Mean.Age+PatchPtoA10+PatchArea10+PatchArea50+Mean_Elevation, data=bir5, Gamma('identity'))
m12<-glm(RRSD2~1+Continent_Island+Mean.Age+PatchArea10+PatchArea50+Mean_Elevation+Mean_Rainfall, data=bir5, Gamma('identity'))
m13<-glm(RRSD2~1+Mean.Age+PatchArea10+PatchArea50+p.Forest25km+Mean_Elevation, data=bir5, Gamma('identity'))
m14<-glm(RRSD2~1+Mean.Age+PatchArea10+PatchArea50+nPatches50+Mean_Elevation, data=bir5, Gamma('identity'))
m15<-glm(RRSD2~1+Mean.Age+PatchArea10+PatchArea50+p.Forest5km+Mean_Elevation, data=bir5, Gamma('identity'))
m16<-glm(RRSD2~1+Continent_Island+Mean.Age+PatchArea10+PatchArea50+p.Forest25km+Mean_Elevation, data=bir5, Gamma('identity'))
m17<-glm(RRSD2~1+Mean.Age+PatchPtoA25+PatchArea10+PatchArea50+Mean_Elevation, data=bir5, Gamma('identity'))
m18<-glm(RRSD2~1+Mean.Age+PatchArea10+PatchArea50+nPatches5+Mean_Elevation, data=bir5, Gamma('identity'))
m19<-glm(RRSD2~1+Mean.Age+PatchArea25+PatchArea10+PatchArea50+Mean_Elevation, data=bir5, Gamma('identity'))
m20<-glm(RRSD2~1+Mean.Age+PatchPtoA10+PatchArea10+PatchArea50+nPatches50+Mean_Elevation, data=bir5, Gamma('identity'))
m21<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea10+PatchArea50+Mean_Elevation, data=bir5, Gamma('identity'))
m22<-glm(RRSD2~1+Continent_Island+Mean.Age+PatchArea25+PatchArea10+PatchArea50+Mean_Elevation, data=bir5, Gamma('identity'))
m23<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea10+Mean_Elevation+Mean_Rainfall, data=bir5, Gamma('identity'))
m24<-glm(RRSD2~1+Mean.Age+PatchPtoA10+PatchArea10+PatchArea50+Mean_Elevation+Mean_Rainfall, data=bir5, Gamma('identity'))
m25<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea25+PatchArea10+p.Forest25km, data=bir5, Gamma('identity'))
m26<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea10+p.Forest25km+Mean_Elevation, data=bir5, Gamma('identity'))
m27<-glm(RRSD2~1+Continent_Island+Mean.Age+PatchPtoA25+PatchArea10+PatchArea50+Mean_Elevation, data=bir5, Gamma('identity'))
m28<-glm(RRSD2~1+Continent_Island+Mean.Age+PatchArea10+PatchArea50+nPatches5+Mean_Elevation, data=bir5, Gamma('identity'))
m29<-glm(RRSD2~1+Continent_Island+Mean.Age+PatchArea10+PatchArea50+p.Forest5km+Mean_Elevation, data=bir5, Gamma('identity'))
m30<-glm(RRSD2~1+Continent_Island+Mean.Age+PatchPtoA75+PatchArea10+PatchArea50+Mean_Elevation, data=bir5, Gamma('identity'))
m31<-glm(RRSD2~1+Mean.Age+PatchPtoA25+PatchPtoA10+PatchArea10+PatchArea50+Mean_Elevation, data=bir5, Gamma('identity'))
m32<-glm(RRSD2~1+Mean.Age+PatchPtoA10+PatchArea10+PatchArea50+nPatches5+Mean_Elevation, data=bir5, Gamma('identity'))
m33<-glm(RRSD2~1+Continent_Island+Mean.Age+PatchArea10+PatchArea50+nPatches50+Mean_Elevation, data=bir5, Gamma('identity'))
m34<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea25+PatchArea10+PatchArea50+Mean_Elevation, data=bir5, Gamma('identity'))
m35<-glm(RRSD2~1+Mean.Age+PatchPtoA10+PatchArea25+PatchArea10+PatchArea50+Mean_Elevation, data=bir5, Gamma('identity'))
m36<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea10+p.Forest25km, data=bir5, Gamma('identity'))
m37<-glm(RRSD2~1+Mean.Age+PatchPtoA10+PatchArea10+PatchArea50+p.Forest25km+Mean_Elevation, data=bir5, Gamma('identity'))
m38<-glm(RRSD2~1+Mean.Age+PatchPtoA10+PatchArea10+PatchArea50+p.Forest5km+Mean_Elevation, data=bir5, Gamma('identity'))
m39<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchArea10+PatchArea50+p.Forest25km+Mean_Elevation, data=bir5, Gamma('identity'))
m40<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA10+PatchArea10+PatchArea50+p.Forest25km+Mean_Elevation, data=bir5, Gamma('identity'))
m41<-glm(RRSD2~1+Continent_Island+Mean.Age+PatchArea25+PatchArea50+Mean_Elevation, data=bir5, Gamma('identity'))
m42<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea25+PatchArea10+Mean_Rainfall, data=bir5, Gamma('identity'))
m43<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA10+PatchArea10+PatchArea50+nPatches5+Mean_Elevation, data=bir5, Gamma('identity'))
m44<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea25+PatchArea10+Mean_Elevation+Mean_Rainfall, data=bir5, Gamma('identity'))
m45<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchArea10+PatchArea50+Mean_Elevation+Mean_Rainfall, data=bir5, Gamma('identity'))
m46<-glm(RRSD2~1+Continent_Island+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea25+PatchArea10, data=bir5, Gamma('identity'))
m47<-glm(RRSD2~1+Mean.Age+PatchArea10+PatchArea50+nPatches5+Mean_Elevation+Mean_Rainfall, data=bir5, Gamma('identity'))
m48<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchPtoA10+PatchArea10+Mean_Elevation, data=bir5, Gamma('identity'))
m49<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea25+nPatches5+Mean_Elevation+Mean_Rainfall, data=bir5, Gamma('identity'))
m50<-glm(RRSD2~1+Mean.Age+PatchArea10+PatchArea50+nPatches50+Mean_Elevation+Mean_Rainfall, data=bir5, Gamma('identity'))
m51<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA10+PatchArea10+PatchArea50+Mean_Elevation+Mean_Rainfall, data=bir5, Gamma('identity'))
m52<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea25+PatchArea10+PatchArea50, data=bir5, Gamma('identity'))
m53<-glm(RRSD2~1+Continent_Island+Mean.Age+PatchPtoA75+PatchPtoA10+PatchArea10+PatchArea50+Mean_Elevation, data=bir5, Gamma('identity'))
m54<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA10+PatchArea10+PatchArea50+nPatches50+Mean_Elevation, data=bir5, Gamma('identity'))
m55<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchPtoA10+PatchArea10+PatchArea50+Mean_Elevation, data=bir5, Gamma('identity'))
m56<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchArea10+PatchArea50+p.Forest5km+Mean_Elevation, data=bir5, Gamma('identity'))
m57<-glm(RRSD2~1+Mean.Age+PatchArea10+PatchArea50+p.Forest5km+Mean_Elevation+Mean_Rainfall, data=bir5, Gamma('identity'))
m58<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea25+PatchArea10+p.Forest25km+Mean_Elevation, data=bir5, Gamma('identity'))
m59<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchPtoA10+PatchArea25+PatchArea10, data=bir5, Gamma('identity'))
m60<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA10+PatchArea25+PatchArea10+PatchArea50+Mean_Elevation, data=bir5, Gamma('identity'))
m61<-glm(RRSD2~1+Continent_Island+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea25+PatchArea10+Mean_Elevation, data=bir5, Gamma('identity'))
m62<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea10+p.Forest25km+nPatches50, data=bir5, Gamma('identity'))
m63<-glm(RRSD2~1+Mean.Age+PatchPtoA25+PatchArea10+PatchArea50+Mean_Elevation+Mean_Rainfall, data=bir5, Gamma('identity'))
m64<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA10+PatchArea10+PatchArea50+p.Forest25km, data=bir5, Gamma('identity'))
m65<-glm(RRSD2~1+Mean.Age+PatchArea25+PatchArea10+PatchArea50+Mean_Elevation+Mean_Rainfall, data=bir5, Gamma('identity'))
m66<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea25+PatchArea10+nPatches5+Mean_Elevation, data=bir5, Gamma('identity'))
m67<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea10+p.Forest5km+Mean_Elevation, data=bir5, Gamma('identity'))
m68<-glm(RRSD2~1+Mean.Age+PatchArea10+PatchArea50+p.Forest25km+Mean_Elevation+Mean_Rainfall, data=bir5, Gamma('identity'))
m69<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea25+PatchArea10+p.Forest5km, data=bir5, Gamma('identity'))
m70<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea25, data=bir5, Gamma('identity'))
m71<-glm(RRSD2~1+Continent_Island+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea25, data=bir5, Gamma('identity'))
m72<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea10+nPatches5+Mean_Elevation, data=bir5, Gamma('identity'))
m73<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA10+PatchArea10+PatchArea50+p.Forest5km+Mean_Elevation, data=bir5, Gamma('identity'))
m74<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea25+Mean_Rainfall, data=bir5, Gamma('identity'))
m75<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchArea10+PatchArea50+nPatches5+Mean_Elevation, data=bir5, Gamma('identity'))
m76<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea25+PatchArea10+nPatches5, data=bir5, Gamma('identity'))
m77<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchPtoA10+PatchArea25, data=bir5, Gamma('identity'))
m78<-glm(RRSD2~1+Continent_Island+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea10+Mean_Elevation, data=bir5, Gamma('identity'))
m79<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea25+p.Forest25km, data=bir5, Gamma('identity'))
m80<-glm(RRSD2~1+Continent_Island+Mean.Age+PatchPtoA10+PatchArea25+PatchArea10+PatchArea50+Mean_Elevation, data=bir5, Gamma('identity'))
m81<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea10+nPatches5+Mean_Elevation+Mean_Rainfall, data=bir5, Gamma('identity'))
m82<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea25+PatchArea10+nPatches50, data=bir5, Gamma('identity'))
m83<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea10+PatchArea50+p.Forest25km+Mean_Elevation, data=bir5, Gamma('identity'))
m84<-glm(RRSD2~1+Continent_Island+Mean.Age+PatchPtoA10+PatchArea10+PatchArea50+Mean_Elevation+Mean_Rainfall, data=bir5, Gamma('identity'))
m85<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchPtoA10+PatchArea10+p.Forest25km, data=bir5, Gamma('identity'))
m86<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea10+nPatches50+Mean_Elevation, data=bir5, Gamma('identity'))
m87<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea25+PatchArea10+p.Forest5km+Mean_Elevation, data=bir5, Gamma('identity'))
m88<-glm(RRSD2~1+Continent_Island+Mean.Age+PatchPtoA75+PatchPtoA25+PatchPtoA10+PatchArea25, data=bir5, Gamma('identity'))
m89<-glm(RRSD2~1+Continent_Island+Mean.Age+PatchArea25+PatchArea50+nPatches5+Mean_Elevation, data=bir5, Gamma('identity'))
m90<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchArea10+PatchArea50+nPatches50+Mean_Elevation, data=bir5, Gamma('identity'))
m91<-glm(RRSD2~1+Continent_Island+Mean.Age+PatchArea10+PatchArea50+nPatches5+Mean_Elevation+Mean_Rainfall, data=bir5, Gamma('identity'))
m92<-glm(RRSD2~1+Mean.Age+PatchPtoA75+nPatches5+Mean_Elevation+Mean_Rainfall, data=bir5, Gamma('identity'))
m93<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchArea25+PatchArea10+PatchArea50+Mean_Elevation, data=bir5, Gamma('identity'))
m94<-glm(RRSD2~1+Mean.Age+PatchArea25+PatchArea50+Mean_Elevation, data=bir5, Gamma('identity'))
m95<-glm(RRSD2~1+Continent_Island+Mean.Age+PatchPtoA25+PatchArea25+PatchArea10+PatchArea50+Mean_Elevation, data=bir5, Gamma('identity'))
m96<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchArea10+PatchArea50+Mean_Elevation+Mean_Rainfall, data=bir5, Gamma('identity'))
m97<-glm(RRSD2~1+Mean.Age+PatchArea25+PatchArea50+Mean_Elevation+Mean_Rainfall, data=bir5, Gamma('identity'))
m98<-glm(RRSD2~1+Continent_Island+Mean.Age+PatchPtoA10+PatchArea10+PatchArea50+nPatches5+Mean_Elevation, data=bir5, Gamma('identity'))
m99<-glm(RRSD2~1+PatchArea10+PatchArea50+Mean_Elevation, data=bir5, Gamma('identity'))
m100<-glm(RRSD2~1+Continent_Island+Mean.Age+PatchPtoA10+PatchArea10+PatchArea50+p.Forest25km+Mean_Elevation, data=bir5, Gamma('identity'))
m101<-glm(RRSD2~1+Mean.Age+PatchPtoA75+PatchPtoA25+PatchPtoA10+PatchArea25+PatchArea10+Mean_Elevation, data=bir5, Gamma('identity'))

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
					length(bir5$RRSD2))

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
#m3 (RRSD2~1+Mean.Age+PatchArea10+PatchArea50+Mean_Elevation)
m2b<-lm(RRSD2~1+Mean.Age+PatchArea10+PatchArea50+Mean_Elevation, data=bir4)
summary(m2b)
Am2b<-anova(m2b)
Am2bss<-Am2b$"Sum Sq"
print(cbind(Am2b,PctExp=Am2bss/sum(Am2bss)*100)) #PctExp: % explain variation

#~~~~~~~~~~~~~~
#Mammals
taxa<-subset(mydata,Taxa=="Mammals")
#Exclude the with species richness 0 (the RRSD can't be calculated)
mam<-subset(taxa,RRSD!="NA")
nrow(mam)#31
summary(mam)
#Absolute Response Ratio
mam$RRSD2<-abs(mam$RRSD)
#Remove outlyers
outlyersmam<-boxplot(mam$RRSD2)
sort(outlyersmam$out) # >1.6
mam2<-mam
boxplot(mam2$RRSD2)
mam2$RRSD2<-mam2$RRSD2+0.5

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
				"PatchPtoA5","PatchPtoA10","PatchPtoA25","PatchPtoA100","PatchPtoA50",
				"PatchArea5","PatchArea75",
				"nPatches5","nPatches25","nPatches50",
				"p.Forest10km",
				"Continent_Island","Zoogeographical_Realm_Holt_etal_2013",
				"Mean_Elevation","Mean_Rainfall")], 1, anyNA),]
nrow(mam4)#n=24
## Good of fitness and normality
global.model<-glm(RRSD2 ~1+SuccStage+Mean.Age+PreviousLandUse+
				PatchPtoA5+PatchPtoA10+PatchPtoA25+PatchPtoA100+PatchPtoA50+
				PatchArea5+PatchArea75+
				nPatches5+nPatches25+nPatches50+
				p.Forest10km+
				Continent_Island+Zoogeographical_Realm_Holt_etal_2013+
				Mean_Elevation+Mean_Rainfall,data=mam4, Gamma("identity"))
plot(global.model)
hist(resid(global.model), breaks=20)
#1- Model selection
model <- glmulti(RRSD2 ~1+SuccStage+Mean.Age+PreviousLandUse+
				PatchPtoA5+PatchPtoA10+PatchPtoA25+PatchPtoA100+PatchPtoA50+
				PatchArea5+PatchArea75+
				nPatches5+nPatches25+nPatches50+
				p.Forest10km+
				Continent_Island+Zoogeographical_Realm_Holt_etal_2013+
				Mean_Elevation+Mean_Rainfall,data=mam4,
               level=1, family=Gamma (link="identity"), crit="aicc")
ranks<-weightable(model)
ranks
plot(model,type="s")
ranks7<-ranks[ranks$aicc <= min(ranks$aicc) + 7,]
ranks7
#73 models + null

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
m1<-glm(RRSD2~1, data=mam5, Gamma('identity'))
m2<-glm(RRSD2~1+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Elevation+Mean_Rainfall, data=mam5, Gamma('identity'))
m3<-glm(RRSD2~1+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Rainfall, data=mam5, Gamma('identity'))
m4<-glm(RRSD2~1+PatchPtoA100+PatchPtoA50+nPatches25+nPatches50+p.Forest10km+Mean_Elevation+Mean_Rainfall, data=mam5, Gamma('identity'))
m5<-glm(RRSD2~1+PatchPtoA100+PatchPtoA50+nPatches25+p.Forest10km+Mean_Elevation+Mean_Rainfall, data=mam5, Gamma('identity'))
m6<-glm(RRSD2~1+PatchPtoA25+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Rainfall, data=mam5, Gamma('identity'))
m7<-glm(RRSD2~1+Mean.Age+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Elevation+Mean_Rainfall, data=mam5, Gamma('identity'))
m8<-glm(RRSD2~1+PatchPtoA100+PatchPtoA50+nPatches25+p.Forest10km+Mean_Rainfall, data=mam5, Gamma('identity'))
m9<-glm(RRSD2~1+Mean.Age+PatchPtoA25+PatchPtoA100+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m10<-glm(RRSD2~1+PatchPtoA10+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Elevation+Mean_Rainfall, data=mam5, Gamma('identity'))
m11<-glm(RRSD2~1+SuccStage+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Rainfall, data=mam5, Gamma('identity'))
m12<-glm(RRSD2~1+Continent_Island+PatchPtoA10+PatchArea5+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m13<-glm(RRSD2~1+Mean.Age+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Rainfall, data=mam5, Gamma('identity'))
m14<-glm(RRSD2~1+PatchPtoA25+PatchPtoA100+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m15<-glm(RRSD2~1+Mean.Age+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m16<-glm(RRSD2~1+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m17<-glm(RRSD2~1+PatchPtoA100+PatchPtoA50+nPatches5+p.Forest10km+Mean_Elevation+Mean_Rainfall, data=mam5, Gamma('identity'))
m18<-glm(RRSD2~1+PatchPtoA100+PatchPtoA50+PatchArea75+p.Forest10km+Mean_Rainfall, data=mam5, Gamma('identity'))
m19<-glm(RRSD2~1+PatchPtoA100+PatchPtoA50+nPatches50+p.Forest10km+Mean_Rainfall, data=mam5, Gamma('identity'))
m20<-glm(RRSD2~1+PatchPtoA100+PatchPtoA50+nPatches25+nPatches50+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m21<-glm(RRSD2~1+PatchPtoA100+PatchPtoA50+nPatches50+p.Forest10km+Mean_Elevation+Mean_Rainfall, data=mam5, Gamma('identity'))
m22<-glm(RRSD2~1+PatchPtoA10+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Rainfall, data=mam5, Gamma('identity'))
m23<-glm(RRSD2~1+PatchPtoA5+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Elevation+Mean_Rainfall, data=mam5, Gamma('identity'))
m24<-glm(RRSD2~1+PatchPtoA25+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Elevation+Mean_Rainfall, data=mam5, Gamma('identity'))
m25<-glm(RRSD2~1+PatchPtoA100+PatchPtoA50+nPatches5+p.Forest10km+Mean_Rainfall, data=mam5, Gamma('identity'))
m26<-glm(RRSD2~1+Continent_Island+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Rainfall, data=mam5, Gamma('identity'))
m27<-glm(RRSD2~1+PatchPtoA100+PatchPtoA50+PatchArea75+p.Forest10km+Mean_Elevation+Mean_Rainfall, data=mam5, Gamma('identity'))
m28<-glm(RRSD2~1+Continent_Island+PatchPtoA25+PatchPtoA100+nPatches5+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m29<-glm(RRSD2~1+PatchPtoA100+PatchPtoA50+PatchArea5+p.Forest10km+Mean_Elevation+Mean_Rainfall, data=mam5, Gamma('identity'))
m30<-glm(RRSD2~1+PatchPtoA5+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Rainfall, data=mam5, Gamma('identity'))
m31<-glm(RRSD2~1+PatchPtoA100+PatchPtoA50+PatchArea5+p.Forest10km+Mean_Rainfall, data=mam5, Gamma('identity'))
m32<-glm(RRSD2~1+Continent_Island+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Elevation+Mean_Rainfall, data=mam5, Gamma('identity'))
m33<-glm(RRSD2~1+PatchPtoA25+PatchPtoA100+nPatches5+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m34<-glm(RRSD2~1+PatchPtoA25+PatchPtoA100+nPatches25+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m35<-glm(RRSD2~1+SuccStage+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Elevation+Mean_Rainfall, data=mam5, Gamma('identity'))
m36<-glm(RRSD2~1+PatchPtoA100+PatchPtoA50+nPatches25+nPatches50+p.Forest10km+Mean_Rainfall, data=mam5, Gamma('identity'))
m37<-glm(RRSD2~1+Mean.Age+PatchPtoA10+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Elevation+Mean_Rainfall, data=mam5, Gamma('identity'))
m38<-glm(RRSD2~1+Continent_Island+PatchPtoA100+PatchPtoA50+PatchArea5+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m39<-glm(RRSD2~1+Mean.Age+PatchPtoA25+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m40<-glm(RRSD2~1+PatchPtoA100+PatchPtoA50+PatchArea75+nPatches25+nPatches50+p.Forest10km+Mean_Rainfall, data=mam5, Gamma('identity'))
m41<-glm(RRSD2~1+PatchPtoA100+PatchPtoA50+nPatches25+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m42<-glm(RRSD2~1+Continent_Island+PatchPtoA25+PatchPtoA50+nPatches5+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m43<-glm(RRSD2~1+PreviousLandUse+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Rainfall, data=mam5, Gamma('identity'))
m44<-glm(RRSD2~1+Continent_Island+PatchPtoA25+PatchPtoA100+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m45<-glm(RRSD2~1+PatchPtoA100+PatchPtoA50+nPatches5+nPatches25+nPatches50+p.Forest10km+Mean_Elevation+Mean_Rainfall, data=mam5, Gamma('identity'))
m46<-glm(RRSD2~1+Continent_Island+PatchPtoA100+PatchPtoA50+nPatches25+nPatches50+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m47<-glm(RRSD2~1+Mean.Age+PatchPtoA25+PatchPtoA100+nPatches25+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m48<-glm(RRSD2~1+PatchPtoA10+PatchPtoA25+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Rainfall, data=mam5, Gamma('identity'))
m49<-glm(RRSD2~1+Mean.Age+PatchPtoA25+PatchPtoA100+nPatches5+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m50<-glm(RRSD2~1+Continent_Island+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m51<-glm(RRSD2~1+PatchPtoA100+PatchPtoA50+PatchArea5+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m52<-glm(RRSD2~1+PatchPtoA25+PatchPtoA100+PatchPtoA50+nPatches25+p.Forest10km+Mean_Rainfall, data=mam5, Gamma('identity'))
m53<-glm(RRSD2~1+PatchPtoA25+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m54<-glm(RRSD2~1+PatchPtoA25+PatchPtoA100+nPatches50+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m55<-glm(RRSD2~1+Continent_Island+PatchPtoA10+PatchArea5+nPatches50+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m56<-glm(RRSD2~1+PatchPtoA100+PatchPtoA50+nPatches5+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m57<-glm(RRSD2~1+PatchPtoA10+PatchPtoA100+PatchPtoA50+nPatches25+p.Forest10km+Mean_Elevation+Mean_Rainfall, data=mam5, Gamma('identity'))
m58<-glm(RRSD2~1+Mean.Age+PatchPtoA25+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Elevation+Mean_Rainfall, data=mam5, Gamma('identity'))
m59<-glm(RRSD2~1+Continent_Island+PatchPtoA25+PatchPtoA100+nPatches25+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m60<-glm(RRSD2~1+PatchPtoA25+PatchPtoA100+p.Forest10km+Mean_Elevation+Mean_Rainfall, data=mam5, Gamma('identity'))
m61<-glm(RRSD2~1+Continent_Island+PatchPtoA25+PatchArea5+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m62<-glm(RRSD2~1+PatchPtoA25+PatchPtoA100+PatchArea75+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m63<-glm(RRSD2~1+PatchPtoA25+PatchPtoA100+nPatches5+p.Forest10km+Mean_Elevation+Mean_Rainfall, data=mam5, Gamma('identity'))
m64<-glm(RRSD2~1+Continent_Island+PatchPtoA10+PatchPtoA100+PatchArea5+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m65<-glm(RRSD2~1+Continent_Island+PatchPtoA25+PatchPtoA100+PatchArea5+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m66<-glm(RRSD2~1+Continent_Island+PatchPtoA25+PatchPtoA100+Mean_Elevation, data=mam5, Gamma('identity'))
m67<-glm(RRSD2~1+Mean.Age+PatchPtoA25+PatchPtoA100+nPatches50+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m68<-glm(RRSD2~1+Mean.Age+PatchPtoA25+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Rainfall, data=mam5, Gamma('identity'))
m69<-glm(RRSD2~1+PatchPtoA100+PatchPtoA50+PatchArea5+nPatches25+nPatches50+p.Forest10km+Mean_Elevation+Mean_Rainfall, data=mam5, Gamma('identity'))
m70<-glm(RRSD2~1+PatchPtoA25+PatchPtoA100+PatchPtoA50+nPatches25+p.Forest10km+Mean_Elevation+Mean_Rainfall, data=mam5, Gamma('identity'))
m71<-glm(RRSD2~1+Continent_Island+PatchPtoA10+PatchArea5+nPatches5+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m72<-glm(RRSD2~1+PatchPtoA25+PatchPtoA100+PatchPtoA50+nPatches25+nPatches50+p.Forest10km+Mean_Rainfall, data=mam5, Gamma('identity'))
m73<-glm(RRSD2~1+PatchPtoA25+PatchPtoA50+nPatches5+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))
m74<-glm(RRSD2~1+PatchPtoA100+PatchPtoA50+PatchArea75+p.Forest10km+Mean_Elevation, data=mam5, Gamma('identity'))

ranks<-m.AICc(list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,
					m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39,m40,
					m41,m42,m43,m44,m45,m46,m47,m48,m49,m50,m51,m52,m53,m54,m55,m56,m57,m58,m59,m60,
					m61,m62,m63,m64,m65,m66,m67,m68,m69,m70,m71,m72,m73,m74),
					as.integer(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20",
					"21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40",
					"41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60",
					"61","62","63","64","65","66","67","68","69","70","71","72","73","74")),length(mam5$RRSD2))

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
#m4 (RRSD2~1+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Rainfall)
m2b<-lm(RRSD2~1+PatchPtoA100+PatchPtoA50+p.Forest10km+Mean_Rainfall, data=mam4)
summary(m2b)
m2b<-anova(m2b)
m2bss<-m2b$"Sum Sq"
print(cbind(m2b,PctExp=m2bss/sum(m2bss)*100)) #PctExp: % explain variation

##Extra fig
library(ggplot2)
library(ggrepel)
ggplot(mam4,aes(x=Mean_Rainfall, y=RRSD2))+
	geom_point(size=5)+
	geom_smooth(method="lm",se=F)+
  geom_text_repel(aes(label=mam4$Name), size=6, segment.color='#cccccc', arrow = arrow(length=unit(0.01,'npc')),force=10)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = 	element_blank())
