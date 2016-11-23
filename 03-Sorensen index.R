#Calculating of Sorensen similarity index for 69 studies of 76
  #and Morisita-Horn similarity index for 31 studies of 76
install.packages("vegan")
library(vegan)

study <- read.csv("C:/Users/QUECA/Dropbox/Tesis MSc/Capitulo 1 - Review/Articulos base/R-beta data/FS52.csv", 
                  row.names=1)
MorHorn <- vegdist(study, method="horn", binary=FALSE, diag = TRUE, upper=FALSE)
MorHorn#dissimilarity of Morisita-Horn

study <- read.csv("C:/Users/QUECA/Dropbox/Tesis MSc/Capitulo 1 - Review/Articulos base/R-beta data/S4.csv", 
                  row.names=1)
Sor <- vegdist(study, binary=TRUE, diag = TRUE, upper=FALSE)
Sor#dissimilarity of Sorensen

#each value was compiled in the mydata.csv matrix (my data), plus the other 72 studies