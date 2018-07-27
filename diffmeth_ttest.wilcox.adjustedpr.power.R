###ttest_Stg4_Gast_OL_27###

## PACKAGES
## Package for adapting the format to perform ttest
library("reshape2")
## Package to calculate the power of analysis
library("pwr")

library("effsize")


## MAIN
#source <- read.csv(file="StgGast_vs_rest.diffmeth.ORG.csv", sep=";", dec=",", head=T, na.strings = "NA")
#if dataset where columns are separated by the character ";" - the decimals are indicated by "," - the names of the columns are specified (header = TRUE) - the lines do not all have the same number of variables -> "NA" are added
source <- read.csv(file="StgGast_vs_rest.diffmeth.ORG.csv", head=T, na.strings = "NA")
#if dataset where columns are separated by the character "," - the decimals are indicated by "." - the names of the columns are specified (header = TRUE) - the lines do not all have the same number of variables -> "NA" are added


source$ttest.4groups.p.adj<- p.adjust(source$Pr, method = "fdr")

ttest.4groups.pval<-c()
ttest.3groups.pval<-c()

wilcox.4groups.pval<-c()
wilcox.3groups.pval<-c()

#norm_R1<-c()
#norm_R2<-c()
#norm_R3<-c()
#norm_R4<-c()

norm_RS<-c()

ttest.d4<-c()
ttest.d3<-c()

for (x in source$line){
#for (x in 1:5){
  
  R<-c(source[x,13], source[x,14], source[x,15])
  S3<-c(source[x,19], source[x,20], source[x,21], source[x,22], source[x,23], source[x,24])
  S4<-c(source[x,26], source[x,17], source[x,18], source[x,19], source[x,20], source[x,21], source[x,22], source[x,23], source[x,24])
  all3<-c(R,S3)
  all4<-c(R,S4)
  #  R3<-c(source[x,23], source[x,24], source[x,25])
#  R4<-c(source[x,26], source[x,27], source[x,28])

  all3names<-c(rep("R",3),rep("S",6))
  Rs3<-data.frame(all3names, all3)
  all4names<-c(rep("R",3),rep("S",9))
  Rs4<-data.frame(all4names, all4)
  
## Capturing all four or just three [R2-R4] groups in a data frame per row ##  
#  Rs4<-data.frame(R,S4)
#  Rs3<-data.frame(R,S3)

## Test for normality
  
#  shapR1<-shapiro.test(Rs4$R1)
#  shapR2<-shapiro.test(Rs4$R2)
#  shapR3<-shapiro.test(Rs4$R3)
#  shapR4<-shapiro.test(Rs4$R4)
  
  
#  norm_R1[x]<-shapR1[2]
#  norm_R2[x]<-shapR2[2]
#  norm_R3[x]<-shapR3[2]
#  norm_R4[x]<-shapR4[2]
  
## Convert format for ttest with melt ##  
#  fRs4<-melt(Rs4, measure.vars=c("R","S4"))
#  fRs3<-melt(Rs3, measure.vars=c("R","S3"))
  
  
## test for normality all R groups together
  
  shapR<-shapiro.test(all3) 
#  norm_R[x] <-shapR[2]
  if (shapR[2]>0.05){
    norm_RS[x]<-c("Normal USE ttest")
  }else{
    norm_RS[x]<-c("Not normal USE wilcox")
  }
## We actually calculate ttest (ttest) including the 3 groups only, since for all four it is already in source$Pr ##  
  ttest_out3 <- t.test(x=R, y=S3, na.action=na.exclude) 
  ttest.3groups.pval[x]<- ttest_out3[3]

  
## Now, for Wilcoxon-Mann_Whitney, we calculate pval for 4 and 3 groups (Stg4 excluded) ##
  wilcox_out4 <- wilcox.test(x=R, y=S4, na.action=na.exclude, paired=FALSE) 
  wilcox.4groups.pval[x]<- wilcox_out4[3]
 
  wilcox_out3 <- wilcox.test(x=R, y=S3, na.action=na.exclude, paired=FALSE) 
  wilcox.3groups.pval[x]<- wilcox_out3[3] 
  
  
## Now, for the calculation of power for ttest
  d4<-cohen.d(Rs4$all4, Rs4$all4names, na.rm=TRUE)
  ttest.d4[x]<-d4[3]
  d3<-cohen.d(Rs3$all3, Rs3$all3names, na.rm=TRUE)
  ttest.d3[x]<-d3[3]
}

#source$normtestR1 <- unlist(norm_R1)
#source$normtestR2 <- unlist(norm_R2)
#source$normtestR3 <- unlist(norm_R3)
#source$normtestR4 <- unlist(norm_R4)

source$normtestRall<- unlist(norm_RS)

source$ttest.3groups.pval<- unlist(ttest.3groups.pval)

source$ttest.3groups.p.adj<- p.adjust(ttest.3groups.pval, method = "fdr")

source$wilcox.4groups.pval<-unlist(wilcox.4groups.pval)

source$wilcox.4groups.p.adj<- p.adjust(source$wilcox.4groups.pval, method = "fdr")

source$wilcox.3groups.pval<-unlist(wilcox.3groups.pval)

source$wilcox.3groups.p.adj<- p.adjust(source$wilcox.3groups.pval, method = "fdr")

source$cohens.d.4<-unlist(ttest.d4)

source$cohens.d.3<-unlist(ttest.d3)


## Power part
## Custom functions:

#create a new function = equivalent of standard.deviation.p in excel (.../n) - because population size is known (= "n" exact is known)
#sd.p=function(x,...){sd(x)*sqrt((length(x)-1)/length(x))}
#calculate and extract the specific value of power
powrt3 <- function(x){

  pw <- pwr.t2n.test(d=x, n1=3,n2=6,sig.level=0.05) 
  return(pw$power)
}
powrt4 <- function(x){
  
  pw <- pwr.t2n.test(d=x, n1=3,n2=9,sig.level=0.05) 
  return(pw$power)
}


source$ttest.pwr.4groups <- powrt4(source$cohens.d.4)
source$ttest.pwr.3groups <- powrt3(source$cohens.d.3)

### Finally, write results onto a file (both in US and european formats)
write.csv(source, file="StgGast_vs_rest.diffmeth.ORG.ttest.wilcox.adjpval.withNAs.csv", na="NA")
write.table(source, file="StgGast_vs_rest.diffmeth.ORG.ttest.wilcox.adjpval.withNAs.EUR.csv", sep=";", na="NA")

