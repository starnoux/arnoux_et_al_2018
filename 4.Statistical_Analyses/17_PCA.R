
---
title : 'Production of PCA plot (on Genetic diversity)'
author: "S. Arnoux"
date : "January 2018"
output : html_document
---

```r

rm(list = ls())

setwd("/Users/stephaniearnoux/Documents/Solution/vcf/VCF_CLEAN_RMD_NoMiss/")
library("SNPRelate")

MMvcf.fn<-"eggplant.vcf"
PMvcf.fn<-"pepper.vcf"
LAvcf.fn<-"Tomato.vcf"
MM_C = "S. melongena" 
MM_W = "S. insanum"
PM_C = "C. annuum" 
PM_W = "C. frustescens \n & chinense"
LA_C = "S. lycopersicum" 
LA_W = "S. peruvianum"
MM_col.list <- c("#A164D0","dimgrey") # Eggplant Colours
PM_col.list <- c("#FF8C00","dimgrey") # Pepper Colours
LA_col.list <- c("firebrick2","dimgrey") # Tomato Colours

#### MM ####
snpgdsVCF2GDS(MMvcf.fn, "MMccm.gds",  method="biallelic.only")
MMgenofile <- snpgdsOpen("MMccm.gds")
MMsample.id <- read.gdsn(index.gdsn(MMgenofile, "sample.id"))
MMsample.id
MMpop_code <- scan("MM_Group.txt", what=character())
MMbind=data.frame( MMsample.id, MMpop_code)

MMccm_pca<-snpgdsPCA(MMgenofile, sample.id=levels(droplevels(MMbind[MMbind$MMpop_code != "X",]$MMsample.id)), autosome.only=FALSE)
MMtab <- data.frame(sample.id = MMccm_pca$sample.id, 
                       pop = factor(MMpop_code)[match(MMccm_pca$sample.id, MMsample.id)],
                       EV1 = MMccm_pca$eigenvect[,1],    # the first eigenvector
                       EV2 = MMccm_pca$eigenvect[,2],    # the second eigenvector
                       stringsAsFactors = FALSE)
MMpc.percent <- MMccm_pca$varprop*100
head(round(MMpc.percent, 2))

snpgdsClose(MMgenofile)


#### PM ####
snpgdsVCF2GDS(PMvcf.fn, "PMccm.gds",  method="biallelic.only")
PMgenofile <- snpgdsOpen("PMccm.gds")
PMsample.id <- read.gdsn(index.gdsn(PMgenofile, "sample.id"))
PMsample.id

PMpop_code <- scan("PM_Group.txt", what=character())
PMbind =data.frame( PMsample.id, PMpop_code)
PMccm_pca<-snpgdsPCA(PMgenofile,sample.id=levels(droplevels(PMbind[PMbind$PMpop_code != "X",]$PMsample.id)), autosome.only=FALSE, num.thread=4)
PMtab <- data.frame(sample.id = PMccm_pca$sample.id, 
pop = factor(PMpop_code)[match(PMccm_pca$sample.id, PMsample.id)],
EV1 = PMccm_pca$eigenvect[,1],    # the first eigenvector
EV2 = PMccm_pca$eigenvect[,2],    # the second eigenvector
stringsAsFactors = FALSE)
PMpc.percent <- PMccm_pca$varprop*100
head(round(PMpc.percent, 2))

snpgdsClose(MMgenofile)

#### LA ####
snpgdsVCF2GDS(LAvcf.fn, "LAccm.gds",  method="biallelic.only")
LAgenofile <- snpgdsOpen("LAccm.gds")
LAsample.id <- read.gdsn(index.gdsn(LAgenofile, "sample.id"))
LAsample.id

LApop_code <- scan("LA_Group.txt", what=character())
LAbind= data.frame( LAsample.id, LApop_code)

LAccm_pca<-snpgdsPCA(LAgenofile,sample.id=levels(droplevels(LAbind[LAbind$LApop_code != "X",]$LAsample.id)), autosome.only=FALSE, num.thread=4)
LAtab <- data.frame(sample.id = LAccm_pca$sample.id, 
                       pop = factor(LApop_code)[match(LAccm_pca$sample.id, LAsample.id)],
                       EV1 = LAccm_pca$eigenvect[,1],    # the first eigenvector
                       EV2 = LAccm_pca$eigenvect[,2],    # the second eigenvector
                       stringsAsFactors = FALSE)
LApc.percent <- LAccm_pca$varprop*100
head(round(LApc.percent, 2))

snpgdsClose(LAgenofile)


##########################################
##########################################
###########     Plot    ##################
##########################################
##########################################

#### WARNING ####     #####     #####     #####     #####     #####     #####     #####     #####      #####
############################################################################################################
#### WARNING #### HERE BELOW YOU NEED TO CHECK AND MAKE THE FIGURE WITH THE SPECIFIC AXIS. #################
############################################################################################################
#### WARNING ####     #####     #####     #####     #####     #####     #####     #####     #####      #####

pch.list =c(18,20)
imageoutput=file.path("~/Documents/Solution/vcf/VCF_CLEAN_RMD_NoMiss/",paste("PCA_Pi_woTitle.tiff", sep = ""))
tiff(file=imageoutput,height = 10, width = 23, units = 'cm',res = 300)
layout(matrix(c(1,2,3), nrow = 1, ncol = 3, byrow = TRUE))
par(oma = c(1.5,1.5,2,1),mar=c(4.6,5.1,2,1))
plot(MMtab$EV2, MMtab$EV1, col=MM_col.list[as.integer(MMtab$pop)], cex.lab=1.4, xlab="", ylab="eigenvector 1",pch=pch.list[as.integer(MMtab$pop)],bty="n",xaxt="n",yaxt="n",cex=2,cex.lab=2)
axis(1, at=c(-0.4,-0.2,0,0.2,0.4,0.6) ,cex.axis=1.2) #0.9
axis(2, at=c(-0.2,0,0.2,0.4), cex.axis=1.2) #0.8

legend("right", legend=c(MM_C,MM_W), pch=c(21,23), col=MM_col.list[1:2],bty="n",pt.cex=1.5,cex=1.2)
title(main= "Eggplant",lwd=2,lty=5,cex.main=2)

par(mar=c(4.6,2,2,1))
plot(PMtab$EV2, PMtab$EV1, col=PM_col.list[as.integer(PMtab$pop)], cex.lab=1.4, xlab="eigenvector 2", ylab="",ylim=c(-0.2,0.8),pch=pch.list[as.integer(PMtab$pop)],bty="n",xaxt="n",yaxt="n",cex=2,cex.lab=2)
axis(1, at=c(-0.2,0,0.2,0.4,0.6,0.8) ,cex.axis=1.2) #0.9
axis(2, at=c(-0.2,0,0.2,0.4,0.6,0.8), cex.axis=1.2) #0.8

legend("topright", legend=c(PM_C,PM_W), pch=c(23,21), col=PM_col.list[1:2],bty="n",pt.cex=1.5,cex=1.2)
title(main= "Pepper",lwd=2,lty=5,cex.main=2)

par(mar=c(4.6,2,2,1))
plot(LAtab$EV2, LAtab$EV1, col=LA_col.list[as.integer(LAtab$pop)], cex.lab=1.4, xlab="", ylab="",pch=pch.list[as.integer(LAtab$pop)],bty="n",xaxt="n",yaxt="n",cex=2,cex.lab=2)
#text(LAtab$EV2, LAtab$EV1, labels = LAtab$sample.id, pos = 4)
axis(1, at=c(-0.4,-0.2,0,0.2,0.4,0.6) ,cex.axis=1.2) #0.9
axis(2, at=c(-0.3,-0.2,-0.1,0,0.1,0.2), cex.axis=1.2) #0.8

legend("right", legend=c(LA_C,LA_W), pch=c(23,21), col=LA_col.list[1:2],bty="n",pt.cex=1.5,cex=1.2)
title(main= "Tomato",lwd=2,lty=5,cex.main=2)
#title(main= "PCA on SNPs",lwd=2,lty=5,font.main=4,cex.main=2, outer=TRUE)
dev.off()

```
