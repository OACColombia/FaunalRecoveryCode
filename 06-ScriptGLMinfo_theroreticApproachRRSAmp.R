#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Bootstrapping of the Response ratio in Species composition similarity
#Clear memory
rm(list=ls())
library(plyr)  # load library
setwd("c:/Users/Flaco/Dropbox/Tesis MSc/Capitulo 1 - Review/Articulos base/Results/GEB/Submitted/Resubmission")
#Call the file MainData.csv
mydata <- read.csv("c:/Users/Flaco/Dropbox/Tesis MSc/Capitulo 1 - Review/Articulos base/Results/GEB/Submitted/Resubmission/MainData.csv")
#~~~~~~~~~~~~~~
#Amphibians
taxa<-subset(mydata,Taxa=="Amphibians")
#Exclude the forest with no list of species (the RRS can't be calculated), and the OGF data
data1<-subset(taxa,RRS!="NA")
data<-subset(data1,Forest.cov!="OGF")
nrow(data)#47
summary(data)
#Sign inverted
data$RRS2<-data$RRS*(-1)
#Remove outlyers
outlyers<-boxplot(data$RRS2)
sort(outlyers$out) # >1.179
data2<-subset(data, data$RRS2<1.179)
boxplot(data2$RRS2)
## Good of fitness and normality
names(data2)
data2$RRS2<-data2$RRS2+0.5
global.model<-glm(RRS2~Age +Forest.cov +Continent_Island +Biogeographic_realms +Nbiome +Mean_Elevation +Mean_Rainfall +Matrix +Previous_Use, data=data2, Gamma("identity"))
plot(global.model)
hist(resid(global.model), breaks=30)
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
## Run models via bootstrap
top_ranked<-list()
null_model_AICc<-list()
null_model_position<-list()
#coef<-list()
#r2<-list()
for (i in 1:10000){
data3<-ddply(data2,.(Study),randomRows,1)
m1<-glm(RRS2~1, data=data3, Gamma("identity"))
m2<-glm(RRS2~Age, data=data3, Gamma("identity"))
m3<-glm(RRS2~Forest.cov, data=data3, Gamma("identity"))
m4<-glm(RRS2~Continent_Island, data=data3, Gamma("identity"))
m5<-glm(RRS2~Biogeographic_realms, data=data3, Gamma("identity"))
m6<-glm(RRS2~Nbiome, data=data3, Gamma("identity"))
m7<-glm(RRS2~Mean_Elevation, data=data3, Gamma("identity"))
m8<-glm(RRS2~Mean_Rainfall, data=data3, Gamma("identity"))
m9<-glm(RRS2~Matrix, data=data3, Gamma("identity"))
m10<-glm(RRS2~Previous_Use, data=data3, Gamma("identity"))
m11<-glm(RRS2~Age+Forest.cov, data=data3, Gamma("identity"))
m12<-glm(RRS2~Age+Continent_Island, data=data3, Gamma("identity"))

m13<-glm(RRS2~Age+Biogeographic_realms, data=data3, Gamma("identity"))
m14<-glm(RRS2~Age+Nbiome, data=data3, Gamma("identity"))
m15<-glm(RRS2~Age+Mean_Elevation, data=data3, Gamma("identity"))
m16<-glm(RRS2~Age+Mean_Rainfall, data=data3, Gamma("identity"))
m17<-glm(RRS2~Age+Matrix, data=data3, Gamma("identity"))
m18<-glm(RRS2~Age+Previous_Use, data=data3, Gamma("identity"))
m19<-glm(RRS2~Forest.cov+Continent_Island, data=data3, Gamma("identity"))

m20<-glm(RRS2~Forest.cov+Biogeographic_realms, data=data3, Gamma("identity"))
m21<-glm(RRS2~Forest.cov+Nbiome, data=data3, Gamma("identity"))
m22<-glm(RRS2~Forest.cov+Mean_Elevation, data=data3, Gamma("identity"))
m23<-glm(RRS2~Forest.cov+Mean_Rainfall, data=data3, Gamma("identity"))
m24<-glm(RRS2~Forest.cov+Matrix, data=data3, Gamma("identity"))
m25<-glm(RRS2~Forest.cov+Previous_Use, data=data3, Gamma("identity"))

m26<-glm(RRS2~Continent_Island+Biogeographic_realms, data=data3, Gamma("identity"))
m27<-glm(RRS2~Continent_Island+Nbiome, data=data3, Gamma("identity"))
m28<-glm(RRS2~Continent_Island+Mean_Elevation, data=data3, Gamma("identity"))
m29<-glm(RRS2~Continent_Island+Mean_Rainfall, data=data3, Gamma("identity"))
m30<-glm(RRS2~Continent_Island+Matrix, data=data3, Gamma("identity"))
m31<-glm(RRS2~Continent_Island+Previous_Use, data=data3, Gamma("identity"))

m32<-glm(RRS2~Biogeographic_realms+Nbiome, data=data3, Gamma("identity"))
m33<-glm(RRS2~Biogeographic_realms+Mean_Elevation, data=data3, Gamma("identity"))
m34<-glm(RRS2~Biogeographic_realms+Mean_Rainfall, data=data3, Gamma("identity"))
m35<-glm(RRS2~Biogeographic_realms+Matrix, data=data3, Gamma("identity"))
m36<-glm(RRS2~Biogeographic_realms+Previous_Use, data=data3, Gamma("identity"))
m37<-glm(RRS2~Nbiome+Mean_Elevation, data=data3, Gamma("identity"))
m38<-glm(RRS2~Nbiome+Mean_Rainfall, data=data3, Gamma("identity"))
m39<-glm(RRS2~Nbiome+Matrix, data=data3, Gamma("identity"))
m40<-glm(RRS2~Nbiome+Previous_Use, data=data3, Gamma("identity"))
m41<-glm(RRS2~Mean_Elevation+Mean_Rainfall, data=data3, Gamma("identity"))
m42<-glm(RRS2~Mean_Elevation+Matrix, data=data3, Gamma("identity"))
m43<-glm(RRS2~Mean_Elevation+Previous_Use, data=data3, Gamma("identity"))
m44<-glm(RRS2~Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m45<-glm(RRS2~Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m46<-glm(RRS2~Matrix+Previous_Use, data=data3, Gamma("identity"))
m47<-glm(RRS2~Age+Forest.cov+Continent_Island, data=data3, Gamma("identity"))

m48<-glm(RRS2~Age+Forest.cov+Biogeographic_realms, data=data3, Gamma("identity"))
m49<-glm(RRS2~Age+Forest.cov+Nbiome, data=data3, Gamma("identity"))
m50<-glm(RRS2~Age+Forest.cov+Mean_Elevation, data=data3, Gamma("identity"))
m51<-glm(RRS2~Age+Forest.cov+Mean_Rainfall, data=data3, Gamma("identity"))
m52<-glm(RRS2~Age+Forest.cov+Matrix, data=data3, Gamma("identity"))
m53<-glm(RRS2~Age+Forest.cov+Previous_Use, data=data3, Gamma("identity"))

m54<-glm(RRS2~Age+Continent_Island+Biogeographic_realms, data=data3, Gamma("identity"))
m55<-glm(RRS2~Age+Continent_Island+Nbiome, data=data3, Gamma("identity"))
m56<-glm(RRS2~Age+Continent_Island+Mean_Elevation, data=data3, Gamma("identity"))
m57<-glm(RRS2~Age+Continent_Island+Mean_Rainfall, data=data3, Gamma("identity"))
m58<-glm(RRS2~Age+Continent_Island+Matrix, data=data3, Gamma("identity"))
m59<-glm(RRS2~Age+Continent_Island+Previous_Use, data=data3, Gamma("identity"))

m60<-glm(RRS2~Age+Biogeographic_realms+Nbiome, data=data3, Gamma("identity"))
m61<-glm(RRS2~Age+Biogeographic_realms+Mean_Elevation, data=data3, Gamma("identity"))
m62<-glm(RRS2~Age+Biogeographic_realms+Mean_Rainfall, data=data3, Gamma("identity"))
m63<-glm(RRS2~Age+Biogeographic_realms+Matrix, data=data3, Gamma("identity"))
m64<-glm(RRS2~Age+Biogeographic_realms+Previous_Use, data=data3, Gamma("identity"))
m65<-glm(RRS2~Age+Nbiome+Mean_Elevation, data=data3, Gamma("identity"))
m66<-glm(RRS2~Age+Nbiome+Mean_Rainfall, data=data3, Gamma("identity"))
m67<-glm(RRS2~Age+Nbiome+Matrix, data=data3, Gamma("identity"))
m68<-glm(RRS2~Age+Nbiome+Previous_Use, data=data3, Gamma("identity"))
m69<-glm(RRS2~Age+Mean_Elevation+Mean_Rainfall, data=data3, Gamma("identity"))
m70<-glm(RRS2~Age+Mean_Elevation+Matrix, data=data3, Gamma("identity"))
m71<-glm(RRS2~Age+Mean_Elevation+Previous_Use, data=data3, Gamma("identity"))
m72<-glm(RRS2~Age+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m73<-glm(RRS2~Age+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m74<-glm(RRS2~Age+Matrix+Previous_Use, data=data3, Gamma("identity"))

m75<-glm(RRS2~Forest.cov+Continent_Island+Biogeographic_realms, data=data3, Gamma("identity"))
m76<-glm(RRS2~Forest.cov+Continent_Island+Nbiome, data=data3, Gamma("identity"))
m77<-glm(RRS2~Forest.cov+Continent_Island+Mean_Elevation, data=data3, Gamma("identity"))
m78<-glm(RRS2~Forest.cov+Continent_Island+Mean_Rainfall, data=data3, Gamma("identity"))
m79<-glm(RRS2~Forest.cov+Continent_Island+Matrix, data=data3, Gamma("identity"))
m80<-glm(RRS2~Forest.cov+Continent_Island+Previous_Use, data=data3, Gamma("identity"))

m81<-glm(RRS2~Forest.cov+Biogeographic_realms+Nbiome, data=data3, Gamma("identity"))
m82<-glm(RRS2~Forest.cov+Biogeographic_realms+Mean_Elevation, data=data3, Gamma("identity"))
m83<-glm(RRS2~Forest.cov+Biogeographic_realms+Mean_Rainfall, data=data3, Gamma("identity"))
m84<-glm(RRS2~Forest.cov+Biogeographic_realms+Matrix, data=data3, Gamma("identity"))
m85<-glm(RRS2~Forest.cov+Biogeographic_realms+Previous_Use, data=data3, Gamma("identity"))
m86<-glm(RRS2~Forest.cov+Nbiome+Mean_Elevation, data=data3, Gamma("identity"))
m87<-glm(RRS2~Forest.cov+Nbiome+Mean_Rainfall, data=data3, Gamma("identity"))
m88<-glm(RRS2~Forest.cov+Nbiome+Matrix, data=data3, Gamma("identity"))
m89<-glm(RRS2~Forest.cov+Nbiome+Previous_Use, data=data3, Gamma("identity"))
m90<-glm(RRS2~Forest.cov+Mean_Elevation+Mean_Rainfall, data=data3, Gamma("identity"))
m91<-glm(RRS2~Forest.cov+Mean_Elevation+Matrix, data=data3, Gamma("identity"))
m92<-glm(RRS2~Forest.cov+Mean_Elevation+Previous_Use, data=data3, Gamma("identity"))
m93<-glm(RRS2~Forest.cov+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m94<-glm(RRS2~Forest.cov+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m95<-glm(RRS2~Forest.cov+Matrix+Previous_Use, data=data3, Gamma("identity"))

m96<-glm(RRS2~Continent_Island+Biogeographic_realms+Nbiome, data=data3, Gamma("identity"))
m97<-glm(RRS2~Continent_Island+Biogeographic_realms+Mean_Elevation, data=data3, Gamma("identity"))
m98<-glm(RRS2~Continent_Island+Biogeographic_realms+Mean_Rainfall, data=data3, Gamma("identity"))
m99<-glm(RRS2~Continent_Island+Biogeographic_realms+Matrix, data=data3, Gamma("identity"))
m100<-glm(RRS2~Continent_Island+Biogeographic_realms+Previous_Use, data=data3, Gamma("identity"))
m101<-glm(RRS2~Continent_Island+Nbiome+Mean_Elevation, data=data3, Gamma("identity"))
m102<-glm(RRS2~Continent_Island+Nbiome+Mean_Rainfall, data=data3, Gamma("identity"))
m103<-glm(RRS2~Continent_Island+Nbiome+Matrix, data=data3, Gamma("identity"))
m104<-glm(RRS2~Continent_Island+Nbiome+Previous_Use, data=data3, Gamma("identity"))
m105<-glm(RRS2~Continent_Island+Mean_Elevation+Mean_Rainfall, data=data3, Gamma("identity"))
m106<-glm(RRS2~Continent_Island+Mean_Elevation+Matrix, data=data3, Gamma("identity"))
m107<-glm(RRS2~Continent_Island+Mean_Elevation+Previous_Use, data=data3, Gamma("identity"))
m108<-glm(RRS2~Continent_Island+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m109<-glm(RRS2~Continent_Island+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m110<-glm(RRS2~Continent_Island+Matrix+Previous_Use, data=data3, Gamma("identity"))

m111<-glm(RRS2~Biogeographic_realms+Nbiome+Mean_Elevation, data=data3, Gamma("identity"))
m112<-glm(RRS2~Biogeographic_realms+Nbiome+Mean_Rainfall, data=data3, Gamma("identity"))
m113<-glm(RRS2~Biogeographic_realms+Nbiome+Matrix, data=data3, Gamma("identity"))
m114<-glm(RRS2~Biogeographic_realms+Nbiome+Previous_Use, data=data3, Gamma("identity"))
m115<-glm(RRS2~Biogeographic_realms+Mean_Elevation+Mean_Rainfall, data=data3, Gamma("identity"))
m116<-glm(RRS2~Biogeographic_realms+Mean_Elevation+Matrix, data=data3, Gamma("identity"))
m117<-glm(RRS2~Biogeographic_realms+Mean_Elevation+Previous_Use, data=data3, Gamma("identity"))
m118<-glm(RRS2~Biogeographic_realms+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m119<-glm(RRS2~Biogeographic_realms+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m120<-glm(RRS2~Biogeographic_realms+Matrix+Previous_Use, data=data3, Gamma("identity"))
m121<-glm(RRS2~Nbiome+Mean_Elevation+Mean_Rainfall, data=data3, Gamma("identity"))
m122<-glm(RRS2~Nbiome+Mean_Elevation+Matrix, data=data3, Gamma("identity"))
m123<-glm(RRS2~Nbiome+Mean_Elevation+Previous_Use, data=data3, Gamma("identity"))
m124<-glm(RRS2~Nbiome+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m125<-glm(RRS2~Nbiome+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m126<-glm(RRS2~Nbiome+Matrix+Previous_Use, data=data3, Gamma("identity"))
m127<-glm(RRS2~Mean_Elevation+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m128<-glm(RRS2~Mean_Elevation+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m129<-glm(RRS2~Mean_Rainfall+Matrix+Previous_Use, data=data3, Gamma("identity"))

m130<-glm(RRS2~Age+Forest.cov+Continent_Island+Biogeographic_realms, data=data3, Gamma("identity"))
m131<-glm(RRS2~Age+Forest.cov+Continent_Island+Nbiome, data=data3, Gamma("identity"))
m132<-glm(RRS2~Age+Forest.cov+Continent_Island+Mean_Elevation, data=data3, Gamma("identity"))
m133<-glm(RRS2~Age+Forest.cov+Continent_Island+Mean_Rainfall, data=data3, Gamma("identity"))
m134<-glm(RRS2~Age+Forest.cov+Continent_Island+Matrix, data=data3, Gamma("identity"))
m135<-glm(RRS2~Age+Forest.cov+Continent_Island+Previous_Use, data=data3, Gamma("identity"))

m136<-glm(RRS2~Age+Forest.cov+Biogeographic_realms+Nbiome, data=data3, Gamma("identity"))
m137<-glm(RRS2~Age+Forest.cov+Biogeographic_realms+Mean_Elevation, data=data3, Gamma("identity"))
m138<-glm(RRS2~Age+Forest.cov+Biogeographic_realms+Mean_Rainfall, data=data3, Gamma("identity"))
m139<-glm(RRS2~Age+Forest.cov+Biogeographic_realms+Matrix, data=data3, Gamma("identity"))
m140<-glm(RRS2~Age+Forest.cov+Biogeographic_realms+Previous_Use, data=data3, Gamma("identity"))
m141<-glm(RRS2~Age+Forest.cov+Nbiome+Mean_Elevation, data=data3, Gamma("identity"))
m142<-glm(RRS2~Age+Forest.cov+Nbiome+Mean_Rainfall, data=data3, Gamma("identity"))
m143<-glm(RRS2~Age+Forest.cov+Nbiome+Matrix, data=data3, Gamma("identity"))
m144<-glm(RRS2~Age+Forest.cov+Nbiome+Previous_Use, data=data3, Gamma("identity"))
m145<-glm(RRS2~Age+Forest.cov+Mean_Elevation+Mean_Rainfall, data=data3, Gamma("identity"))
m146<-glm(RRS2~Age+Forest.cov+Mean_Elevation+Matrix, data=data3, Gamma("identity"))
m147<-glm(RRS2~Age+Forest.cov+Mean_Elevation+Previous_Use, data=data3, Gamma("identity"))
m148<-glm(RRS2~Age+Forest.cov+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m149<-glm(RRS2~Age+Forest.cov+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m150<-glm(RRS2~Age+Forest.cov+Matrix+Previous_Use, data=data3, Gamma("identity"))

m151<-glm(RRS2~Age+Continent_Island+Biogeographic_realms+Nbiome, data=data3, Gamma("identity"))
m152<-glm(RRS2~Age+Continent_Island+Biogeographic_realms+Mean_Elevation, data=data3, Gamma("identity"))
m153<-glm(RRS2~Age+Continent_Island+Biogeographic_realms+Mean_Rainfall, data=data3, Gamma("identity"))
m154<-glm(RRS2~Age+Continent_Island+Biogeographic_realms+Matrix, data=data3, Gamma("identity"))
m155<-glm(RRS2~Age+Continent_Island+Biogeographic_realms+Previous_Use, data=data3, Gamma("identity"))
m156<-glm(RRS2~Age+Continent_Island+Nbiome+Mean_Elevation, data=data3, Gamma("identity"))
m157<-glm(RRS2~Age+Continent_Island+Nbiome+Mean_Rainfall, data=data3, Gamma("identity"))
m158<-glm(RRS2~Age+Continent_Island+Nbiome+Matrix, data=data3, Gamma("identity"))
m159<-glm(RRS2~Age+Continent_Island+Nbiome+Previous_Use, data=data3, Gamma("identity"))
m160<-glm(RRS2~Age+Continent_Island+Mean_Elevation+Mean_Rainfall, data=data3, Gamma("identity"))
m161<-glm(RRS2~Age+Continent_Island+Mean_Elevation+Matrix, data=data3, Gamma("identity"))
m162<-glm(RRS2~Age+Continent_Island+Mean_Elevation+Previous_Use, data=data3, Gamma("identity"))
m163<-glm(RRS2~Age+Continent_Island+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m164<-glm(RRS2~Age+Continent_Island+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m165<-glm(RRS2~Age+Continent_Island+Matrix+Previous_Use, data=data3, Gamma("identity"))

m166<-glm(RRS2~Age+Biogeographic_realms+Nbiome+Mean_Elevation, data=data3, Gamma("identity"))
m167<-glm(RRS2~Age+Biogeographic_realms+Nbiome+Mean_Rainfall, data=data3, Gamma("identity"))
m168<-glm(RRS2~Age+Biogeographic_realms+Nbiome+Matrix, data=data3, Gamma("identity"))
m169<-glm(RRS2~Age+Biogeographic_realms+Nbiome+Previous_Use, data=data3, Gamma("identity"))
m170<-glm(RRS2~Age+Biogeographic_realms+Mean_Elevation+Mean_Rainfall, data=data3, Gamma("identity"))
m171<-glm(RRS2~Age+Biogeographic_realms+Mean_Elevation+Matrix, data=data3, Gamma("identity"))
m172<-glm(RRS2~Age+Biogeographic_realms+Mean_Elevation+Previous_Use, data=data3, Gamma("identity"))
m173<-glm(RRS2~Age+Biogeographic_realms+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m174<-glm(RRS2~Age+Biogeographic_realms+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m175<-glm(RRS2~Age+Biogeographic_realms+Matrix+Previous_Use, data=data3, Gamma("identity"))
m176<-glm(RRS2~Age+Nbiome+Mean_Elevation+Mean_Rainfall, data=data3, Gamma("identity"))
m177<-glm(RRS2~Age+Nbiome+Mean_Elevation+Matrix, data=data3, Gamma("identity"))
m178<-glm(RRS2~Age+Nbiome+Mean_Elevation+Previous_Use, data=data3, Gamma("identity"))
m179<-glm(RRS2~Age+Nbiome+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m180<-glm(RRS2~Age+Nbiome+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m181<-glm(RRS2~Age+Nbiome+Matrix+Previous_Use, data=data3, Gamma("identity"))
m182<-glm(RRS2~Age+Mean_Elevation+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m183<-glm(RRS2~Age+Mean_Elevation+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m184<-glm(RRS2~Age+Mean_Elevation+Matrix+Previous_Use, data=data3, Gamma("identity"))
m185<-glm(RRS2~Age+Mean_Rainfall+Matrix+Previous_Use, data=data3, Gamma("identity"))

m186<-glm(RRS2~Forest.cov+Continent_Island+Biogeographic_realms+Nbiome, data=data3, Gamma("identity"))
m187<-glm(RRS2~Forest.cov+Continent_Island+Biogeographic_realms+Mean_Elevation, data=data3, Gamma("identity"))
m188<-glm(RRS2~Forest.cov+Continent_Island+Biogeographic_realms+Mean_Rainfall, data=data3, Gamma("identity"))
m189<-glm(RRS2~Forest.cov+Continent_Island+Biogeographic_realms+Matrix, data=data3, Gamma("identity"))
m190<-glm(RRS2~Forest.cov+Continent_Island+Biogeographic_realms+Previous_Use, data=data3, Gamma("identity"))
m191<-glm(RRS2~Forest.cov+Continent_Island+Nbiome+Mean_Elevation, data=data3, Gamma("identity"))
m192<-glm(RRS2~Forest.cov+Continent_Island+Nbiome+Mean_Rainfall, data=data3, Gamma("identity"))
m193<-glm(RRS2~Forest.cov+Continent_Island+Nbiome+Matrix, data=data3, Gamma("identity"))
m194<-glm(RRS2~Forest.cov+Continent_Island+Nbiome+Previous_Use, data=data3, Gamma("identity"))
m195<-glm(RRS2~Forest.cov+Continent_Island+Mean_Elevation+Mean_Rainfall, data=data3, Gamma("identity"))
m196<-glm(RRS2~Forest.cov+Continent_Island+Mean_Elevation+Matrix, data=data3, Gamma("identity"))
m197<-glm(RRS2~Forest.cov+Continent_Island+Mean_Elevation+Previous_Use, data=data3, Gamma("identity"))
m198<-glm(RRS2~Forest.cov+Continent_Island+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m199<-glm(RRS2~Forest.cov+Continent_Island+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m200<-glm(RRS2~Forest.cov+Continent_Island+Matrix+Previous_Use, data=data3, Gamma("identity"))

m201<-glm(RRS2~Forest.cov+Biogeographic_realms+Nbiome+Mean_Elevation, data=data3, Gamma("identity"))
m202<-glm(RRS2~Forest.cov+Biogeographic_realms+Nbiome+Mean_Rainfall, data=data3, Gamma("identity"))
m203<-glm(RRS2~Forest.cov+Biogeographic_realms+Nbiome+Matrix, data=data3, Gamma("identity"))
m204<-glm(RRS2~Forest.cov+Biogeographic_realms+Nbiome+Previous_Use, data=data3, Gamma("identity"))
m205<-glm(RRS2~Forest.cov+Biogeographic_realms+Mean_Elevation+Mean_Rainfall, data=data3, Gamma("identity"))
m206<-glm(RRS2~Forest.cov+Biogeographic_realms+Mean_Elevation+Matrix, data=data3, Gamma("identity"))
m207<-glm(RRS2~Forest.cov+Biogeographic_realms+Mean_Elevation+Previous_Use, data=data3, Gamma("identity"))
m208<-glm(RRS2~Forest.cov+Biogeographic_realms+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m209<-glm(RRS2~Forest.cov+Biogeographic_realms+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m210<-glm(RRS2~Forest.cov+Biogeographic_realms+Matrix+Previous_Use, data=data3, Gamma("identity"))
m211<-glm(RRS2~Forest.cov+Nbiome+Mean_Elevation+Mean_Rainfall, data=data3, Gamma("identity"))
m212<-glm(RRS2~Forest.cov+Nbiome+Mean_Elevation+Matrix, data=data3, Gamma("identity"))
m213<-glm(RRS2~Forest.cov+Nbiome+Mean_Elevation+Previous_Use, data=data3, Gamma("identity"))
m214<-glm(RRS2~Forest.cov+Nbiome+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m215<-glm(RRS2~Forest.cov+Nbiome+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m216<-glm(RRS2~Forest.cov+Nbiome+Matrix+Previous_Use, data=data3, Gamma("identity"))
m217<-glm(RRS2~Forest.cov+Mean_Elevation+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m218<-glm(RRS2~Forest.cov+Mean_Elevation+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m219<-glm(RRS2~Forest.cov+Mean_Elevation+Matrix+Previous_Use, data=data3, Gamma("identity"))
m220<-glm(RRS2~Forest.cov+Mean_Rainfall+Matrix+Previous_Use, data=data3, Gamma("identity"))

m221<-glm(RRS2~Continent_Island+Biogeographic_realms+Nbiome+Mean_Elevation, data=data3, Gamma("identity"))
m222<-glm(RRS2~Continent_Island+Biogeographic_realms+Nbiome+Mean_Rainfall, data=data3, Gamma("identity"))
m223<-glm(RRS2~Continent_Island+Biogeographic_realms+Nbiome+Matrix, data=data3, Gamma("identity"))
m224<-glm(RRS2~Continent_Island+Biogeographic_realms+Nbiome+Previous_Use, data=data3, Gamma("identity"))
m225<-glm(RRS2~Continent_Island+Biogeographic_realms+Mean_Elevation+Mean_Rainfall, data=data3, Gamma("identity"))
m226<-glm(RRS2~Continent_Island+Biogeographic_realms+Mean_Elevation+Matrix, data=data3, Gamma("identity"))
m227<-glm(RRS2~Continent_Island+Biogeographic_realms+Mean_Elevation+Previous_Use, data=data3, Gamma("identity"))
m228<-glm(RRS2~Continent_Island+Biogeographic_realms+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m229<-glm(RRS2~Continent_Island+Biogeographic_realms+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m230<-glm(RRS2~Continent_Island+Biogeographic_realms+Matrix+Previous_Use, data=data3, Gamma("identity"))
m231<-glm(RRS2~Continent_Island+Nbiome+Mean_Elevation+Mean_Rainfall, data=data3, Gamma("identity"))
m232<-glm(RRS2~Continent_Island+Nbiome+Mean_Elevation+Matrix, data=data3, Gamma("identity"))
m233<-glm(RRS2~Continent_Island+Nbiome+Mean_Elevation+Previous_Use, data=data3, Gamma("identity"))
m234<-glm(RRS2~Continent_Island+Nbiome+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m235<-glm(RRS2~Continent_Island+Nbiome+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m236<-glm(RRS2~Continent_Island+Nbiome+Matrix+Previous_Use, data=data3, Gamma("identity"))
m237<-glm(RRS2~Continent_Island+Mean_Elevation+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m238<-glm(RRS2~Continent_Island+Mean_Elevation+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m239<-glm(RRS2~Continent_Island+Mean_Elevation+Matrix+Previous_Use, data=data3, Gamma("identity"))
m240<-glm(RRS2~Continent_Island+Mean_Rainfall+Matrix+Previous_Use, data=data3, Gamma("identity"))

m241<-glm(RRS2~Biogeographic_realms+Nbiome+Mean_Elevation+Mean_Rainfall, data=data3, Gamma("identity"))
m242<-glm(RRS2~Biogeographic_realms+Nbiome+Mean_Elevation+Matrix, data=data3, Gamma("identity"))
m243<-glm(RRS2~Biogeographic_realms+Nbiome+Mean_Elevation+Previous_Use, data=data3, Gamma("identity"))
m244<-glm(RRS2~Biogeographic_realms+Nbiome+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m245<-glm(RRS2~Biogeographic_realms+Nbiome+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m246<-glm(RRS2~Biogeographic_realms+Nbiome+Matrix+Previous_Use, data=data3, Gamma("identity"))
m247<-glm(RRS2~Biogeographic_realms+Mean_Elevation+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m248<-glm(RRS2~Biogeographic_realms+Mean_Elevation+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m249<-glm(RRS2~Biogeographic_realms+Mean_Elevation+Matrix+Previous_Use, data=data3, Gamma("identity"))
m250<-glm(RRS2~Biogeographic_realms+Mean_Rainfall+Matrix+Previous_Use, data=data3, Gamma("identity"))
m251<-glm(RRS2~Nbiome+Mean_Elevation+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m252<-glm(RRS2~Nbiome+Mean_Elevation+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m253<-glm(RRS2~Nbiome+Mean_Elevation+Matrix+Previous_Use, data=data3, Gamma("identity"))
m254<-glm(RRS2~Nbiome+Mean_Rainfall+Matrix+Previous_Use, data=data3, Gamma("identity"))
m255<-glm(RRS2~Mean_Elevation+Mean_Rainfall+Matrix+Previous_Use, data=data3, Gamma("identity"))

m256<-glm(RRS2~Age+Biogeographic_realms+Nbiome+Mean_Elevation+Mean_Rainfall, data=data3, Gamma("identity"))
m257<-glm(RRS2~Age+Biogeographic_realms+Nbiome+Mean_Elevation+Matrix, data=data3, Gamma("identity"))
m258<-glm(RRS2~Age+Biogeographic_realms+Nbiome+Mean_Elevation+Previous_Use, data=data3, Gamma("identity"))
m259<-glm(RRS2~Age+Biogeographic_realms+Nbiome+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m260<-glm(RRS2~Age+Biogeographic_realms+Nbiome+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m261<-glm(RRS2~Age+Biogeographic_realms+Nbiome+Matrix+Previous_Use, data=data3, Gamma("identity"))
m262<-glm(RRS2~Age+Biogeographic_realms+Mean_Elevation+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m263<-glm(RRS2~Age+Biogeographic_realms+Mean_Elevation+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m264<-glm(RRS2~Age+Biogeographic_realms+Mean_Elevation+Matrix+Previous_Use, data=data3, Gamma("identity"))
m265<-glm(RRS2~Age+Biogeographic_realms+Mean_Rainfall+Matrix+Previous_Use, data=data3, Gamma("identity"))
m266<-glm(RRS2~Age+Nbiome+Mean_Elevation+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m267<-glm(RRS2~Age+Nbiome+Mean_Elevation+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m268<-glm(RRS2~Age+Nbiome+Mean_Elevation+Matrix+Previous_Use, data=data3, Gamma("identity"))
m269<-glm(RRS2~Age+Nbiome+Mean_Rainfall+Matrix+Previous_Use, data=data3, Gamma("identity"))
m270<-glm(RRS2~Age+Mean_Elevation+Mean_Rainfall+Matrix+Previous_Use, data=data3, Gamma("identity"))

m271<-glm(RRS2~Forest.cov+Continent_Island+Biogeographic_realms+Nbiome+Mean_Elevation, data=data3, Gamma("identity"))
m272<-glm(RRS2~Forest.cov+Continent_Island+Biogeographic_realms+Nbiome+Mean_Rainfall, data=data3, Gamma("identity"))
m273<-glm(RRS2~Forest.cov+Continent_Island+Biogeographic_realms+Nbiome+Matrix, data=data3, Gamma("identity"))
m274<-glm(RRS2~Forest.cov+Continent_Island+Biogeographic_realms+Nbiome+Previous_Use, data=data3, Gamma("identity"))
m275<-glm(RRS2~Forest.cov+Continent_Island+Biogeographic_realms+Mean_Elevation+Mean_Rainfall, data=data3, Gamma("identity"))
m276<-glm(RRS2~Forest.cov+Continent_Island+Biogeographic_realms+Mean_Elevation+Matrix, data=data3, Gamma("identity"))
m277<-glm(RRS2~Forest.cov+Continent_Island+Biogeographic_realms+Mean_Elevation+Previous_Use, data=data3, Gamma("identity"))
m278<-glm(RRS2~Forest.cov+Continent_Island+Biogeographic_realms+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m279<-glm(RRS2~Forest.cov+Continent_Island+Biogeographic_realms+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m280<-glm(RRS2~Forest.cov+Continent_Island+Biogeographic_realms+Matrix+Previous_Use, data=data3, Gamma("identity"))
m281<-glm(RRS2~Forest.cov+Continent_Island+Nbiome+Mean_Elevation+Mean_Rainfall, data=data3, Gamma("identity"))
m282<-glm(RRS2~Forest.cov+Continent_Island+Nbiome+Mean_Elevation+Matrix, data=data3, Gamma("identity"))
m283<-glm(RRS2~Forest.cov+Continent_Island+Nbiome+Mean_Elevation+Previous_Use, data=data3, Gamma("identity"))
m284<-glm(RRS2~Forest.cov+Continent_Island+Nbiome+Matrix+Previous_Use, data=data3, Gamma("identity"))
m285<-glm(RRS2~Forest.cov+Continent_Island+Mean_Elevation+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m286<-glm(RRS2~Forest.cov+Continent_Island+Mean_Elevation+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m287<-glm(RRS2~Forest.cov+Continent_Island+Mean_Elevation+Matrix+Previous_Use, data=data3, Gamma("identity"))
m288<-glm(RRS2~Forest.cov+Continent_Island+Mean_Rainfall+Matrix+Previous_Use, data=data3, Gamma("identity"))

m289<-glm(RRS2~Forest.cov+Biogeographic_realms+Nbiome+Mean_Elevation+Mean_Rainfall, data=data3, Gamma("identity"))
m290<-glm(RRS2~Forest.cov+Biogeographic_realms+Nbiome+Mean_Elevation+Matrix, data=data3, Gamma("identity"))
m291<-glm(RRS2~Forest.cov+Biogeographic_realms+Nbiome+Mean_Elevation+Previous_Use, data=data3, Gamma("identity"))
m292<-glm(RRS2~Forest.cov+Biogeographic_realms+Nbiome+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m293<-glm(RRS2~Forest.cov+Biogeographic_realms+Nbiome+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m294<-glm(RRS2~Forest.cov+Biogeographic_realms+Nbiome+Matrix+Previous_Use, data=data3, Gamma("identity"))
m295<-glm(RRS2~Forest.cov+Biogeographic_realms+Mean_Elevation+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m296<-glm(RRS2~Forest.cov+Biogeographic_realms+Mean_Elevation+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m297<-glm(RRS2~Forest.cov+Biogeographic_realms+Mean_Elevation+Previous_Use, data=data3, Gamma("identity"))
m298<-glm(RRS2~Forest.cov+Biogeographic_realms+Mean_Rainfall+Matrix+Previous_Use, data=data3, Gamma("identity"))
m299<-glm(RRS2~Forest.cov+Nbiome+Mean_Elevation+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m300<-glm(RRS2~Forest.cov+Nbiome+Mean_Elevation+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m301<-glm(RRS2~Forest.cov+Nbiome+Mean_Elevation+Matrix+Previous_Use, data=data3, Gamma("identity"))
m302<-glm(RRS2~Forest.cov+Nbiome+Mean_Rainfall+Matrix+Previous_Use, data=data3, Gamma("identity"))
m303<-glm(RRS2~Forest.cov+Mean_Elevation+Mean_Rainfall+Matrix+Previous_Use, data=data3, Gamma("identity"))

m304<-glm(RRS2~Continent_Island+Biogeographic_realms+Nbiome+Mean_Elevation+Mean_Rainfall, data=data3, Gamma("identity"))
m305<-glm(RRS2~Continent_Island+Biogeographic_realms+Nbiome+Mean_Elevation+Matrix, data=data3, Gamma("identity"))
m306<-glm(RRS2~Continent_Island+Biogeographic_realms+Nbiome+Mean_Elevation+Previous_Use, data=data3, Gamma("identity"))
m307<-glm(RRS2~Continent_Island+Biogeographic_realms+Nbiome+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m308<-glm(RRS2~Continent_Island+Biogeographic_realms+Nbiome+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m309<-glm(RRS2~Continent_Island+Biogeographic_realms+Nbiome+Matrix+Previous_Use, data=data3, Gamma("identity"))
m310<-glm(RRS2~Continent_Island+Biogeographic_realms+Mean_Elevation+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m311<-glm(RRS2~Continent_Island+Biogeographic_realms+Mean_Elevation+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m312<-glm(RRS2~Continent_Island+Biogeographic_realms+Mean_Elevation+Matrix+Previous_Use, data=data3, Gamma("identity"))
m313<-glm(RRS2~Continent_Island+Biogeographic_realms+Mean_Rainfall+Matrix+Previous_Use, data=data3, Gamma("identity"))
m314<-glm(RRS2~Continent_Island+Nbiome+Mean_Elevation+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m315<-glm(RRS2~Continent_Island+Nbiome+Mean_Elevation+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m316<-glm(RRS2~Continent_Island+Nbiome+Mean_Elevation+Matrix+Previous_Use, data=data3, Gamma("identity"))
m317<-glm(RRS2~Continent_Island+Nbiome+Mean_Rainfall+Matrix+Previous_Use, data=data3, Gamma("identity"))
m318<-glm(RRS2~Continent_Island+Mean_Elevation+Mean_Rainfall+Matrix+Previous_Use, data=data3, Gamma("identity"))

m319<-glm(RRS2~Biogeographic_realms+Nbiome+Mean_Elevation+Mean_Rainfall+Matrix, data=data3, Gamma("identity"))
m320<-glm(RRS2~Biogeographic_realms+Nbiome+Mean_Elevation+Mean_Rainfall+Previous_Use, data=data3, Gamma("identity"))
m321<-glm(RRS2~Biogeographic_realms+Nbiome+Mean_Elevation+Matrix+Previous_Use, data=data3, Gamma("identity"))
m322<-glm(RRS2~Biogeographic_realms+Nbiome+Mean_Rainfall+Matrix+Previous_Use, data=data3, Gamma("identity"))
m323<-glm(RRS2~Biogeographic_realms+Mean_Elevation+Mean_Rainfall+Matrix+Previous_Use, data=data3, Gamma("identity"))
m324<-glm(RRS2~Nbiome+Mean_Elevation+Mean_Rainfall+Matrix+Previous_Use, data=data3, Gamma("identity"))

		   
ranks<-m.AICc(list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,
				   m31,m32,m33,m34,m35,m36,m37,m38,m39,m40,m41,m42,m43,m44,m45,m46,m47,m48,m49,m50,m51,m52,m53,m54,m55,m56,m57,m58,m59,m60,
				   m61,m62,m63,m64,m65,m66,m67,m68,m69,m70,m71,m72,m73,m74,m75,m76,m77,m78,m79,m80,m81,m82,m83,m84,m85,m86,m87,m88,m89,m90,
				   m91,m92,m93,m94,m95,m96,m97,m98,m99,m100,
				   m101,m102,m103,m104,m105,m106,m107,m108,m109,m110,m111,m112,m113,m114,m115,m116,m117,m118,m119,m120,m121,m122,m123,m124,m125,m126,m127,m128,m129,m130,
				   m131,m132,m133,m134,m135,m136,m137,m138,m139,m140,m141,m142,m143,m144,m145,m146,m147,m148,m149,m150,m151,m152,m153,m154,m155,m156,m157,m158,m159,m160,
				   m161,m162,m163,m164,m165,m166,m167,m168,m169,m170,m171,m172,m173,m174,m175,m176,m177,m178,m179,m180,m181,m182,m183,m184,m185,m186,m187,m188,m189,m190,
				   m191,m192,m193,m194,m195,m196,m197,m198,m199,m200,
				   m201,m202,m203,m204,m205,m206,m207,m208,m209,m210,m211,m212,m213,m214,m215,m216,m217,m218,m219,m220,m221,m222,m223,m224,m225,m226,m227,m228,m229,m230,
				   m231,m232,m233,m234,m235,m236,m237,m238,m239,m240,m241,m242,m243,m244,m245,m246,m247,m248,m249,m250,m251,m252,m253,m254,m255,m256,m257,m258,m259,m260,
				   m261,m262,m263,m264,m265,m266,m267,m268,m269,m270,m271,m272,m273,m274,m275,m276,m277,m278,m279,m280,m281,m282,m283,m284,m285,m286,m287,m288,m289,m290,
				   m291,m292,m293,m294,m295,m296,m297,m298,m299,m300,
				   m301,m302,m303,m304,m305,m306,m307,m308,m309,m310,m311,m312,m313,m314,m315,m316,m317,m318,m319,m320,m321,m322,m323,m324),
				   as.integer(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25",
				   "26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50",
				   "51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75",
				   "76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100",
				   "101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120","121","122","123","124","125",
				   "126","127","128","129","130","131","132","133","134","135","136","137","138","139","140","141","142","143","144","145","146","147","148","149","150",
				   "151","152","153","154","155","156","157","158","159","160","161","162","163","164","165","166","167","168","169","170","171","172","173","174","175",
				   "176","177","178","179","180","181","182","183","184","185","186","187","188","189","190","191","192","193","194","195","196","197","198","199","200",
				   "201","202","203","204","205","206","207","208","209","210","211","212","213","214","215","216","217","218","219","220","221","222","223","224","225",
				   "226","227","228","229","230","231","232","233","234","235","236","237","238","239","240","241","242","243","244","245","246","247","248","249","250",
				   "251","252","253","254","255","256","257","258","259","260","261","262","263","264","265","266","267","268","269","270","271","272","273","274","275",
				   "276","277","278","279","280","281","282","283","284","285","286","287","288","289","290","291","292","293","294","295","296","297","298","299","300",
				   "301","302","303","304","305","306","307","308","309","310","311","312","313","314","315","316","317","318","319","320","321","322","323","324")),
				   length(data3$RRS2))

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
write.table(df_top_ranked, "Top_ranked.txt",quote=F)
write.table(df_null_model_AICc, "Null_model_AICc.txt",quote=F)
write.table(df_null_model_position, "Null_model_position.txt",quote=F)

## Top-model frequence
plausible_mod_frequency<-as.data.frame(ftable(df_top_ranked[,7]))
plausible_mod_table_AICc<-aggregate(df_top_ranked, list(df_top_ranked[,7]), mean)
plausible_Null_table_AICc<-aggregate(df_null_model_AICc, list(df_null_model_AICc[,7]), mean)
## Saving data.frames
write.table(plausible_mod_frequency, "Plausible_mod_frequency.txt",quote=F, sep="\t")
write.table(plausible_mod_table_AICc, "Plausible_mod_table_AICc_mean.txt",quote=F, sep="\t")
write.table(plausible_Null_table_AICc, "Plausible_Null_table_AICc_mean.txt",quote=F, sep="\t")

#% Variation explained
#Rop Ranked models for Amphibians
plausible_mod_frequency
#First 5:
#1 (5134) <- k = 2
#7 (2427) <- k = 3
#2 ( 612) <- k = 3
#3 ( 452) <- k = 4
#6 ( 440) <- k = 10
m1b<-lm(RRS2~1, data=data2)
summary(m1b)
Am1b<-anova(m1b)
Am1bss<-Am1b$"Sum Sq"
print(cbind(Am1b,PctExp=Am1bss/sum(Am1bss)*100)) #PctExp: % explain variation

m7b<-lm(RRS2~Mean_Elevation, data=data2)
summary(m7b)
Am7b<-anova(m7b)
Am7bss<-Am7b$"Sum Sq"
print(cbind(Am7b,PctExp=Am7bss/sum(Am7bss)*100)) #PctExp: % explain variation

m2b<-lm(RRS2~Age, data=data2)
summary(m2b)
Am2b<-anova(m2b)
Am2bss<-Am2b$"Sum Sq"
print(cbind(Am2b,PctExp=Am2bss/sum(Am2bss)*100)) #PctExp: % explain variation

m3b<-lm(RRS2~Forest.cov, data=data2)
summary(m3b)
Am3b<-anova(m3b)
Am3bss<-Am3b$"Sum Sq"
print(cbind(Am3b,PctExp=Am3bss/sum(Am3bss)*100)) #PctExp: % explain variation

m6b<-lm(RRS2~Nbiome, data=data2)
summary(m6b)
Am6b<-anova(m6b)
Am6bss<-Am6b$"Sum Sq"
print(cbind(Am6b,PctExp=Am6bss/sum(Am6bss)*100)) #PctExp: % explain variation
