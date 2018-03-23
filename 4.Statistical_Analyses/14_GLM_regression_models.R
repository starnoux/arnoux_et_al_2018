---
title : 'Glm on pi according to DEG or DV or BOTH'
author: "S. Arnoux"
date : "January 2018"
output : html_document
---
  
###########################
### LOADING R LIBRARIES ###
###########################

```r

rm(list=ls())

library(missMethyl)
library(edgeR)
library("DESeq2")
library(fields)
library(dplyr)
library(tidyr)
library(Rcpp)
library(data.table)
library(car)
library(lme4)
library(LambertW)

species = "Species_name"

### Data for Covergae 
PATH_STAT ="/path/to/Deseq_Folder/"
FILE = "Exp_Coverage.txt" 
Para_Filter = "/path/to/Filtered_liste.tab.txt"
POP = "Acc_Loc.txt"
setwd(PATH_STAT)

### Data for Egglib Stats
pathwayToPi = "/path/to/DNAsp_results/"
StatCrop = "DNAsp_results_PopA.out"
StatWild = "DNAsp_results_PopB.out"

#Here we filter the accessions that we will not use anyway. Pease is made to be used with SNPs only (not the same tissu RNAseqED.)

#Filter the outgroups 
FILTER = c("O") 

CountTable.cut = read.table(FILE, header=TRUE, row.names=1)
head(CountTable.cut)
Loc = read.table(POP, header=TRUE)
Loc$LocName = paste(Loc$Name,Loc$POP,sep="_") 
names(CountTable.cut)=Loc$LocName
head(CountTable.cut); dim(CountTable.cut)

condition.pop = Loc$POP

#Filter according to the SNPs and not paralogs
parafilter = read.table(Para_Filter,header=FALSE)
para.CountTable.cut =CountTable.cut[rownames(CountTable.cut) %in% parafilter$V1,]

sub.CountTable.cut = subset(para.CountTable.cut, select =  grep(paste(FILTER,collapse="|"),names(para.CountTable.cut),invert=T))
head(sub.CountTable.cut); dim(sub.CountTable.cut)

condition.pop.cut = droplevels(Loc$POP[grep(paste(FILTER,collapse="|"),Loc$LocName,invert=T)])
condition.pop.cut 

level=combn(levels(factor(condition.pop.cut)), 2)

x=level[,1]  

dble.CountTable.cut = subset(sub.CountTable.cut, select = grep(paste("_",x,collapse="|",sep = ""), names(sub.CountTable.cut)))
dble.condition.pop.cut = droplevels(condition.pop.cut[grep(paste(x,collapse="|"),condition.pop.cut)])

cds.cut = DESeqDataSetFromMatrix(countData = dble.CountTable.cut, colData = data.frame(dble.condition.pop.cut), design = ~dble.condition.pop.cut)
cds.cut = DESeq(cds.cut, test = "Wald", fitType = c("parametric"), parallel = TRUE)

res <- results(cds.cut, contrast = c(paste("dble.condition.pop.cut"),as.vector(x)))
dim <- dim(res)
sum <- summary(res)

############################################# 
###############   stats   ###################
############################################# 

Stat_Crop = read.table(paste(pathwayToPi,StatCrop,sep=""),header = TRUE)
colnames(Stat_Crop)[1] <- "Gene"
Stat_Wild = read.table(paste(pathwayToPi,StatWild,sep=""),header = TRUE)
colnames(Stat_Wild)[1] <- "Gene"

## dt frame creation
Stat_Crop$Gene = as.character(Stat_Crop$Gene)
Stat_Wild$Gene = as.character(Stat_Wild$Gene)
dtCrop=data.frame(Gene=Stat_Crop$Gene,PiC=as.numeric(Stat_Crop$Pi),TajimaDC=as.numeric(as.character(Stat_Crop$TajimaD)))
dtWild=data.frame(Gene=Stat_Wild$Gene,PiW=as.numeric(Stat_Wild$Pi),TajimaDW=as.numeric(as.character(Stat_Wild$TajimaD)))
Both_Pi = merge(dtCrop,dtWild, by="Gene")
DEG_data = data.frame(Gene=rownames(res),LFC=res$log2FoldChange)
Delta_st = merge(Both_Pi,DEG_data,by="Gene")
Delta_st$GauLFC = Gaussianize(Delta_st$LFC)[,1]
Delta_st$GauLFC = Delta_st$LFC
Delta_st$chr = paste("Chr.",as.character(substr(Stat_Crop$Gene,7,8)),sep="")
Delta_st = Delta_st[!(Delta_st$chr =='Chr.00'),]

#################
### Pi Filter ###
#################
maxPiWild = quantile(Delta_st$PiW, 1 - 1e-5,na.rm=TRUE)[[1]]
maxPiCrop = quantile(Delta_st$PiC, 1 - 1e-5,na.rm=TRUE)[[1]]

##Filter
Delta_st$PiC[(Delta_st$PiC) >= maxPiCrop] <- NA
Delta_st$PiW[(Delta_st$PiW) >= maxPiWild] <- NA

Pi_only_Delta = na.omit(Delta_st[,c("Gene","PiC","PiW","chr","GauLFC")])
Pi_only_Delta = Pi_only_Delta[!(Pi_only_Delta$PiC==0 & Pi_only_Delta$PiW==0),]
#16376
formula_pi = Gaussianize(Pi_only_Delta$PiW-Pi_only_Delta$PiC)
qqPlot(formula_pi)
formula_pi = Pi_only_Delta$PiW-Pi_only_Delta$PiC

##############  GLM  ############## 

model1a = glm(as.formula(formula_pi ~ GauLFC + chr), data=Pi_only_Delta)
model1c = glm(as.formula(GauLFC ~ formula_pi + chr), data=Pi_only_Delta) # NOT SO GOOD FOR THE DISTRIBUTION GauLFC

summary(model1a)
summary(model1c)

plot(resid(model1a))
plot(resid(model1c))

##############  Creation of the function for the prediction  ############## 
mypredictdf <- function (model, newdata, level=0.95){
  pred <- stats::predict(model, newdata = newdata, se =TRUE, type = "link")
  std <- qnorm(level/2 + 0.5)
  data.frame(newdata,
             y = model$family$linkinv(as.vector(pred$fit)),
             ymin = model$family$linkinv(as.vector(pred$fit - std * pred$se)),
             ymax = model$family$linkinv(as.vector(pred$fit + std * pred$se)), 
             se = as.vector(pred$se))
}

##############  PLOT PREDICTION VS DATA  ############## 
px <- with(Pi_only_Delta, seq(from=min(GauLFC), to=max(GauLFC), length=100))
pdf <- expand.grid(GauLFC=px, chr=unique(Pi_only_Delta$chr))
pdf <- mypredictdf(model1a, newdata=pdf)

g <- ggplot(data=pdf, aes(group=chr))
g <- g + geom_point(data=Pi_only_Delta, aes(x=GauLFC, y=formula_pi),col="grey69")
g <- g + geom_ribbon(aes(x=GauLFC, ymin=ymin, ymax=ymax, color=as.character(chr)),
                     alpha=0.2)
g <- g + geom_line(aes(x=GauLFC, y=y, color=as.character(chr)))
g <- g + facet_wrap(~chr)
g <- g + ggtitle(paste("GLM Predicted probabilities in ",species,sep="")) +
  xlab("LFC Gene Expression") + ylab(expression(paste("Delta ",Pi,sep=""))) 
g <- g + theme(panel.border = element_blank()) + theme_bw()
g <- g + theme(legend.position="none",text=element_text(family="Helvetica",size = 14))

imageoutput=file.path("~/Documents/Solution/vcf/DNAsp/",paste(species,"Pi_GLM.tiff", sep = ""))
tiff(file=imageoutput,height = 17, width = 17, units = 'cm',res = 300)
print(g)
dev.off()


#################
### D Filter ###
#################

D_only_Delta = na.omit(Delta_st[,c("Gene","TajimaDC","TajimaDW","chr","GauLFC")])

formula_D = D_only_Delta$TajimaDW-D_only_Delta$TajimaDC
qqPlot(formula_D)

##############  GLM  ############## 

model2a = glm(as.formula(formula_D ~ GauLFC + chr), data=D_only_Delta)
model2c = glm(as.formula(GauLFC ~ formula_D + chr), data=D_only_Delta) # NOT SO GOOD FOR THE DISTRIBUTION GauLFC

summary(model2a)
summary(model2c)

plot(resid(model2a))
plot(resid(model2c))


##############  PLOT PREDICTION PREDICTION VS DATA  ############## 
px <- with(D_only_Delta, seq(from=min(GauLFC), to=max(GauLFC), length=100))
pdf <- expand.grid(GauLFC=px, chr=unique(D_only_Delta$chr))
pdf <- mypredictdf(model2a, newdata=pdf)

g <- ggplot(data=pdf, aes(group=chr))
g <- g + geom_point(data=D_only_Delta, aes(x=GauLFC, y=formula_D),col="grey69")
g <- g + geom_ribbon(aes(x=GauLFC, ymin=ymin, ymax=ymax, color=as.character(chr)),
                     alpha=0.2)
g <- g + geom_line(aes(x=GauLFC, y=y, color=as.character(chr)))
g <- g + facet_wrap(~chr)
g <- g + ggtitle(paste("GLM Predicted probabilities in ",species,sep="")) +
  xlab("LFC Gene Expression") + ylab(paste("Delta Tajima'D",sep="")) 
g <- g + theme(panel.border = element_blank()) + theme_bw()
g <- g + theme(legend.position="none",text=element_text(family="Helvetica",size = 14))

imageoutput=file.path("~/Documents/Solution/vcf/DNAsp/",paste(species,"TajD_GLM.tiff", sep = ""))
tiff(file=imageoutput,height = 17, width = 17, units = 'cm',res = 300)
print(g)
dev.off()
```
