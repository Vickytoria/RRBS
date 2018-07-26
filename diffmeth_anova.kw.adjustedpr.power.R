###anova_Stg4_Gast_OL_27###

## Package for adapting the format to perform anova
library("reshape2")
## Package to calculate the power of analysis
library("pwr")

### the bit that calculates the anova and extracts the pvalue 
### per line of the file, STARTS HERE

source <- read.csv(file="Stg4_Gast_OL_27_anova.diffmeth.ORG.EUR.csv", sep=";", dec=",", head=T, na.strings = "NA")
#dataset where columns are separated by the character ";" - the decimals are indicated by "," - the names of the columns are specified (header = TRUE) - the lines do not all have the same number of variables -> "NA" are added


source$anova.4groups.p.adj<- p.adjust(source$Pr, method = "fdr")

anova.4groups.pval<-c()
anova.3groups.pval<-c()

kw.4groups.pval<-c()
kw.3groups.pval<-c()

#norm_R1<-c()
#norm_R2<-c()
#norm_R3<-c()
#norm_R4<-c()

norm_R<-c()


for (x in source$line){
#for (x in 1:5){
  R1<-c(source[x,17], source[x,18], source[x,19])
  R2<-c(source[x,20], source[x,21], source[x,22])
  R3<-c(source[x,23], source[x,24], source[x,25])
  R4<-c(source[x,26], source[x,27], source[x,28])

  
## Capturing all four or just three [R2-R4] groups in a data frame per row ##  
  Rs4<-data.frame(R1,R2,R3,R4)
  Rs3<-data.frame(R2,R3,R4)

## Test for normality
  
#  shapR1<-shapiro.test(Rs4$R1)
#  shapR2<-shapiro.test(Rs4$R2)
#  shapR3<-shapiro.test(Rs4$R3)
#  shapR4<-shapiro.test(Rs4$R4)
  
  
#  norm_R1[x]<-shapR1[2]
#  norm_R2[x]<-shapR2[2]
#  norm_R3[x]<-shapR3[2]
#  norm_R4[x]<-shapR4[2]
  
## Convert format for aov with melt ##  
  fRs4<-melt(Rs4, measure.vars=c("R1","R2","R3","R4"))
  fRs3<-melt(Rs3, measure.vars=c("R2","R3","R4"))
  
  
## test for normality all R groups together
  
  shapR<-shapiro.test(fRs4$value) 
#  norm_R[x] <-shapR[2]
  if (shapR[2]>0.05){
    norm_R[x]<-c("Normal USE ANOVA")
  }else{
    norm_R[x]<-c("Not normal USE KW")
  }
## We actually calculate anova (aov) for the 3 groups only, since for all four it is already in source$Pr ##  
  aov_out3 <- aov(fRs3$value ~ fRs3$variable, fRs3, na.action=na.exclude) 
  saov_out3 <- summary(aov_out3)
  msaov_out3 <- matrix(unlist(saov_out3))
  anova.3groups.pval[x]<- msaov_out3[9,1]

  
## Now, for Kruskal-Wallis, we calculate pval for 4 and 3 groups ##
  kw_out4 <- kruskal.test(fRs4$value ~ fRs4$variable, fRs4, na.action=na.exclude) 
  kw.4groups.pval[x]<- kw_out4[3]
 
  kw_out3 <- kruskal.test(fRs3$value ~ fRs3$variable, fRs3, na.action=na.exclude) 
  kw.3groups.pval[x]<- kw_out3[3] 
  
}

#source$normtestR1 <- unlist(norm_R1)
#source$normtestR2 <- unlist(norm_R2)
#source$normtestR3 <- unlist(norm_R3)
#source$normtestR4 <- unlist(norm_R4)

source$normtestRall<- unlist(norm_R)

source<-cbind(source, anova.3groups.pval)

source$anova.3groups.p.adj<- p.adjust(source$anova.3groups.pval, method = "fdr")

source$kw.4groups.pval<-unlist(kw.4groups.pval)

source$kw.4groups.p.adj<- p.adjust(source$kw.4groups.pval, method = "fdr")

source$kw.3groups.pval<-unlist(kw.3groups.pval)

source$kw.3groups.p.adj<- p.adjust(source$kw.3groups.pval, method = "fdr")



## Power part
## Custom functions:

#create a new function = equivalent of standard.deviation.p in excel (.../n) - because population size is known (= "n" exact is known)
sd.p=function(x,...){sd(x)*sqrt((length(x)-1)/length(x))}
#calculate and extract the specific value of power
power <- function(x){
  pw <- pwr.anova.test(x, k=4,n=3,sig.level=0.05) 
  return(pw$power)
}

## Standard deviations are calculated per group and across groups, adding new columns to the input data
source <- data.frame(source, 
                   RSD4=apply(source[c(13:16)],1, sd.p, na.rm = TRUE),
                   RSD3=apply(source[c(14:16)],1, sd.p, na.rm = TRUE), 
                   R1SD=apply(source[c(17:19)],1, sd.p , na.rm = TRUE),
                   R2SD=apply(source[c(20:22)],1, sd.p , na.rm = TRUE), 
                   R3SD=apply(source[c(23:25)],1, sd.p , na.rm = TRUE), 
                   R4SD=apply(source[c(26:28)],1, sd.p , na.rm = TRUE)
)
#R1SD is the std. dev. of the average % of methylation of individuals samples of R1 (R1a, R1b and R1c) - etc...
#head(data, 10)

RxSD4 <- subset(source, select=c(R1SD, R2SD, R3SD, R4SD))
source$SDr4 <- apply(RxSD4,1,mean, na.rm = TRUE)

#SDr is the std. dev. of the replicats = the average of R1SD, R2SD, R3SD, R4SD
RxSD3 <- subset(source, select=c(R2SD, R3SD, R4SD))
source$SDr3 <- apply(RxSD3,1,mean, na.rm = TRUE)


source$f4 <- source$RSD4 / source$SDr4
source$f3 <- source$RSD3 / source$SDr3
#f = effect size use in pwr.anova.test


#"pwr.anova.test" function transformed - we take just the power data
source$ANOVApwr4groups <- power(source$f4)
source$ANOVApwr3groups <- power(source$f3)
#calcul power of each "line"


### Finally, write results onto a file (both in US and european formats)
write.csv(source, file="Stg4_Gast_OL_27_anova.diffmeth.ORG.ANOVA.KW.adjpval.PWRwithNAs.csv", na="NA")
write.table(source, file="Stg4_Gast_OL_27_anova.diffmeth.ORG.ANOVA.KW.adjpval.PWRwithNAs.EUR.csv", sep=";", na="NA")

### TILL HERE!!!!!!


#write.table(data, file="Stg4_Gast_OL_27_anova.diffmeth.ORG.EUR.adjustedpr.power.csv",row.names=F, na="NA",col.names=T, sep=";", dec=",")
#create a new csv with all the last adds

###########################################2
#data <- head(source, 1)
#pour n'utiliser que la 1er ligne de donn?es de "source"