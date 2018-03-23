---
title : 'Genome-Wide figures'
author: "S. Arnoux"
date : "January 2018"
output : html_document
---


###############################################
### STARTING by giving the path to the data ###
###############################################

```r

rm(list=ls())

library(fields)
library(plyr) 
library(data.table)
library("zoo")
library("grDevices")

## Give the working folder
PATH_STAT = "/path/to/DNAsp_Results/"
setwd(PATH_STAT)

## Give the chromosomes locations 
ChrLOC=read.table("/path/to/reference/ref.Gene_Loc.tab")

## Give sample details 
pop2 = "Pop_A_name" 
pop3 = "Pop_B_name"

Species= "species_name"
title = "Prefix_to_the_output_files"

## Settle the file names
Name_File2 = "DNAsp_results_PopA.out" 
Name_File3 = "DNAsp_results_PopB.out"

## Colours for the domesticated population // The wild will be in grey
CropCol_Light = "#A164D0" 
CropCol_Dark = "#5D04A1"

################### NOW DONT TOUCH BUT RUN IT THROUGH
## Reads the files
stat2 =  read.table(Name_File2, header = TRUE)
colnames(stat2)[1] <- "Gene"
stat3 =  read.table(Name_File3, header = TRUE)
colnames(stat3)[1] <- "Gene"

#################
maxPiWild = quantile(stat2$Pi, 1 - 1e-5,na.rm=TRUE)[[1]]
maxPiCrop = quantile(stat3$Pi, 1 - 1e-5,na.rm=TRUE)[[1]]

which((stat2$Pi) >= maxPiCrop)
which((stat3$Pi) >= maxPiWild)

stat2$Pi[(stat2$Pi) >= maxPiCrop] <- NA
which((stat2$Pi) >= maxPiCrop)
stat3$Pi[(stat3$Pi) >= maxPiWild] <- NA
which((stat3$Pi) >= maxPiWild)


## Get the locations 
dtChrLOC = data.table(ChrLOC, key = "V1")
pdtCrop <- data.table(stat2, key = "Gene") #lyc
dtCrop <- pdtCrop[dtChrLOC]
chr <- substr(dtCrop$Gene,7,8)
dtCrop$V2 = as.numeric(as.character(dtCrop$V2))
pFdtCrop <- cbind(dtCrop,chr=as.numeric(as.character(chr)))
#FdtCrop <- na.omit(pFdtCrop)

pdtWild = data.table(stat3, key = "Gene") #pim
dtWild <- pdtWild[dtChrLOC]
chr <- substr(dtWild$Gene,7,8)
dtWild$V2 = as.numeric(as.character(dtWild$V2))
pFdtWild <- cbind(dtWild,chr=as.numeric(as.character(chr)))
#FdtWild <- na.omit(pFdtWild)
pFdtCrop = pFdtCrop[!(pFdtCrop$Pi=='NA'),]
pFdtWild = pFdtWild[!(pFdtWild$Pi=='NA'),]

### Sanety check
#listofx = c(FdtCrop, FdtWild)
numchroms=length(unique(pFdtWild$chr))
chrlist=unique(pFdtWild$chr)[-1]
chrlist2=chrlist
x = pFdtWild
BadGENES=NULL
for (i in 1:length(chrlist)){
  Carp = x[chr==i,]
  for (j in 1:(nrow(Carp)-2)){
    if((Carp[j,]$V2 > Carp[j+1,]$V2) & (Carp[j+1,]$V2 < Carp[j+2,]$V2)){
      BadGENES=rbind(BadGENES, Carp[j+1,])
    }else if((Carp[j,]$V2 < Carp[j+1,]$V2) & (Carp[j+1,]$V2 > Carp[j+2,]$V2)){
      BadGENES=rbind(BadGENES, Carp[j+1,])
    }else{
      BadGENES=BadGENES
    }
  }
  for (j in nrow(Carp)){
    if(Carp[j,]$V2 < Carp[j-1,]$V2){
      BadGENES=rbind(BadGENES, Carp[j,])
    }else{
      BadGENES=BadGENES
    }
  }
}

head(BadGENES) ; nrow(BadGENES)
#113
FdtCrop <- pFdtCrop[!Gene%in%BadGENES$Gene,]
FdtCrop = FdtCrop[order(FdtCrop$V2),]
FdtCrop = FdtCrop[order(FdtCrop$chr),] 
FdtWild <- pFdtWild[!Gene%in%BadGENES$Gene,]
FdtWild = FdtWild[order(FdtWild$V2),]
FdtWild = FdtWild[order(FdtWild$chr),] 


###############################
###############################
###  Stats Pi & Tajima's D  ###
###############################
###############################
## Pi per chr
aggregate(FdtWild$Pi~FdtWild$chr, FUN=mean)[,2]
aggregate(FdtCrop$Pi~FdtCrop$chr, FUN=mean)[,2]
#Global Wild Pi
mean(FdtWild$Pi, na.rm=TRUE) # 0.0005008926 
max(FdtWild$Pi, na.rm=TRUE) #  0.04935065
min(FdtWild$Pi, na.rm=TRUE) # 0
median(FdtWild$Pi, na.rm=TRUE)# 0
# Global Crop Pi
mean(FdtCrop$Pi, na.rm=TRUE) # 0.0003082977
max(FdtCrop$Pi, na.rm=TRUE) # 0.04761905
min(FdtCrop$Pi, na.rm=TRUE) # 0
median(FdtCrop$Pi, na.rm=TRUE) # 0

## D per chr
aggregate(as.numeric(as.character(FdtWild$TajimaD))~FdtWild$chr, FUN=mean)[,2]
aggregate(as.numeric(as.character(FdtCrop$TajimaD))~FdtCrop$chr, FUN=mean)[,2]
#Global Wild D
mean(FdtWild$TajimaD, na.rm=TRUE) #  0.2172857
max(FdtWild$TajimaD, na.rm=TRUE) #  2.779694
min(FdtWild$TajimaD, na.rm=TRUE) #  -2.002529
median(FdtWild$TajimaD, na.rm=TRUE) #  0.236818
# Global Crop D
mean(FdtCrop$TajimaD, na.rm=TRUE) # 0.1425619
max(FdtCrop$TajimaD, na.rm=TRUE) # 2.804832
min(FdtCrop$TajimaD, na.rm=TRUE) # -2.173252
median(FdtCrop$TajimaD, na.rm=TRUE) # -0.3414383


###############################
###############################
### Pi NUCLEOTIDE DIVERSITY ###
###############################
###############################

Rolling.crop = rollmean(FdtCrop$Pi,50,align="left")
data.pi.crop=data.frame(cbind(as.numeric(FdtCrop$chr), as.numeric(FdtCrop$V2), Rolling.crop))
colnames(data.pi.crop) = c("CHROM","POS","PI"); dim(data.pi.crop); head(data.pi.crop)

Rolling.wild = rollmean(FdtWild$Pi,50,align="left")
data.pi.wild=data.frame(cbind(as.numeric(FdtWild$chr), as.numeric(FdtWild$V2), Rolling.wild))
colnames(data.pi.wild) = c("CHROM","POS","PI"); dim(data.pi.wild); head(data.pi.wild)

###############################

CROP=data.pi.crop
WILD=data.pi.wild

###############################
################
## TICKS CROP ##

#attach(CROP)
CROP$pos = NA

# Ticks Creation
ticks=NULL
lastbase=0

for (i in 1:length(chrlist))
{
  if (i==1) 
  {CROP[CROP$CHROM==chrlist[i], ]$pos=CROP[CROP$CHROM==chrlist[i], ]$POS}
  else 
  {
    lastbase=lastbase+tail(subset(CROP,CHROM==chrlist[i-1])$POS, 1)
    CROP[CROP$CHROM==chrlist[i], ]$pos=CROP[CROP$CHROM==chrlist[i], ]$POS+lastbase
  }
  ticks=c(ticks, ((((max(CROP[CROP$CHROM==chrlist[i], ]$pos))-lastbase)/2)+lastbase))
}

################
## TICKS WILD ##

#attach(WILD)
WILD$pos = NA

lastbase=0

for (i in 1:length(chrlist))
{
  if (i==1) 
  {WILD[WILD$CHROM==chrlist[i], ]$pos=WILD[WILD$CHROM==chrlist[i], ]$POS}
  else 
  {
    lastbase=lastbase+tail(subset(WILD,CHROM==chrlist[i-1])$POS, 1)
    WILD[WILD$CHROM==chrlist[i], ]$pos=WILD[WILD$CHROM==chrlist[i], ]$POS+lastbase
  }
}

###############################
#### WARNING ####     #####     #####     #####     #####     #####     #####     #####     #####      #####
############################################################################################################
#### WARNING #### HERE BELOW YOU NEED TO CHECK AND MAKE THE FIGURE WITH THE SPECIFIC AXIS. #################
############################################################################################################
#### WARNING ####     #####     #####     #####     #####     #####     #####     #####     #####      #####

imageoutput=file.path("~/Documents/Solution/vcf/DNAsp/",paste(title,"Pi_GenomeWide.tiff", sep = ""))
tiff(file=imageoutput,height = 10, width = 23, units = 'cm',res = 300)
layout(matrix(c(1,1,1,1,1,2), nrow = 1, ncol = 6, byrow = TRUE))
par(oma = c(2,2,2,0),mar=c(4.1,4.1,1,0))
with(CROP, plot(pos,rep(0.1,nrow(CROP)),xlim=c(0,max(CROP[CROP$CHROM==chrlist[12], ]$pos)),ylim=c(0,0.0075),main="",bty="n",ylab= expression(paste("Mean nucleotide diversity (", pi,")",sep="")),xlab="Chromosomes",xaxt="n",yaxt="n",col="white",cex.lab=1.2))
axis(1, at=ticks, lab=chrlist2, cex.axis=1) #0.9
axis(2, at=c(0,0.0025,0.005,0.0075), cex.axis=1,line = -1.5,las=1) 

#############################
##    SIGNIFICANCE LINES   ##
#i in even chromosomes and j in uneven chromosmess... You might need to check the ylim and to adjust the figures.
#for (i in c(8))
#{
#  rect(min(WILD[WILD$CHROM==chrlist[i],]$pos)-min(WILD[WILD$CHROM==chrlist[i],]$POS),0.2785,max(WILD[WILD$CHROM==chrlist[i],]$pos),0.28,col="dimgrey",border = NA)
#}
#for (j in c(9))
#{
#  rect(min(WILD[WILD$CHROM==chrlist[j],]$pos)-min(WILD[WILD$CHROM==chrlist[j],]$POS),0.2785,max(WILD[WILD$CHROM==chrlist[j],]$pos),0.28,col="darkgrey",border = NA)
#}

################
## Lines WILD ##

for (i in (1:12))
{
  points(WILD[WILD$CHROM==chrlist[((2*i)-1)],]$pos,WILD[WILD$CHROM==chrlist[((2*i)-1)],]$PI,type="l",col="darkgrey",lwd=1.3)
}
for (i in (1:12))
{
  points(WILD[WILD$CHROM==chrlist[2*i],]$pos,WILD[WILD$CHROM==chrlist[2*i],]$PI,type="l",col="dimgrey",lwd=1.3)
}

################
## Lines Crop ##

for (i in (1:12))
{
  i = (2*i)-1
  points(CROP[CROP$CHROM==chrlist[i],]$pos,CROP[CROP$CHROM==chrlist[i],]$PI,type="l",col=CropCol_Light,lwd=1.3)
} #pch=20
for (i in (1:12))
{
  points(CROP[CROP$CHROM==chrlist[2*i],]$pos,CROP[CROP$CHROM==chrlist[2*i],]$PI,type="l",col=CropCol_Dark,lwd=1.3)
}

#################################
### Pi NUC. DIV. DISTRIBUTION ###
#################################

par(mar=c(4.2,0,1,2.1))
plot(x=density(CROP[!CROP$CHROM==0,]$PI)$y,y=density(CROP[!CROP$CHROM==0,]$PI)$x, col="white",ylim=c(0,0.0075),ylab="",main = "",xlab="Density", xaxt="n",yaxt="n",cex.lab=1,col.axis="white",bty="n")
polygon(x=density(WILD[!WILD$CHROM==0,]$PI)$y,y=density(WILD[!WILD$CHROM==0,]$PI)$x, col=adjustcolor( "darkgrey", alpha.f = 0.8), border="dimgrey",lwd=1.3)
polygon(x=density(CROP[!CROP$CHROM==0,]$PI)$y,y=density(CROP[!CROP$CHROM==0,]$PI)$x, col=adjustcolor( CropCol_Light, alpha.f = 0.8), border=CropCol_Dark,lwd=1.3)

legend("topright", inset=0.05, c(pop2,pop3), fill=c(adjustcolor( CropCol_Light, alpha.f = 0.8),adjustcolor( "darkgrey", alpha.f = 0.8)), border=c(CropCol_Dark,"dimgrey"),bty="n" ,cex=1)

title(paste(Species,sep=""),cex.main=1.8, outer=TRUE)

dev.off()


###############################
###############################
###       D Tajimas         ###
###############################
###############################
FdtWild_D = FdtWild[,c("Gene","TajimaD","V2","chr")]
FdtCrop_D = FdtCrop[,c("Gene","TajimaD","V2","chr")]

FdtCrop_D =na.omit(FdtCrop_D)
FdtWild_D =na.omit(FdtWild_D)

Rolling.crop = rollmean(FdtCrop_D$TajimaD,50,align="left")
data.D.crop=data.frame(cbind(as.numeric(FdtCrop_D$chr), as.numeric(FdtCrop_D$V2), Rolling.crop ))
colnames(data.D.crop) = c("CHROM","POS","D"); dim(data.D.crop); head(data.D.crop)

Rolling.wild = rollmean(FdtWild_D$TajimaD,50,align="left")
data.D.wild=data.frame(cbind(as.numeric(FdtWild_D$chr), as.numeric(FdtWild_D$V2), Rolling.wild))
colnames(data.D.wild) = c("CHROM","POS","D"); dim(data.D.wild); head(data.D.wild)

###############################

CROP=data.D.crop
WILD=data.D.wild

###############################
##################
## TICKS CROP ##

#attach(CROP)
CROP$pos = NA

# Ticks Creation
ticks=NULL
lastbase=0

for (i in 1:length(chrlist))
{
  if (i==1) 
  {CROP[CROP$CHROM==chrlist[i], ]$pos=CROP[CROP$CHROM==chrlist[i], ]$POS}
  else 
  {
    lastbase=lastbase+tail(subset(CROP,CHROM==chrlist[i-1])$POS, 1)
    CROP[CROP$CHROM==chrlist[i], ]$pos=CROP[CROP$CHROM==chrlist[i], ]$POS+lastbase
  }
  ticks=c(ticks, ((((max(CROP[CROP$CHROM==chrlist[i], ]$pos))-lastbase)/2)+lastbase))
}


#################
## TICKS WILD ##

#attach(WILD)
WILD$pos = NA

lastbase=0

for (i in 1:length(chrlist))
{
  if (i==1) 
  {WILD[WILD$CHROM==chrlist[i], ]$pos=WILD[WILD$CHROM==chrlist[i], ]$POS}
  else 
  {
    lastbase=lastbase+tail(subset(WILD,CHROM==chrlist[i-1])$POS, 1)
    WILD[WILD$CHROM==chrlist[i], ]$pos=WILD[WILD$CHROM==chrlist[i], ]$POS+lastbase
  }
}

#########
#### WARNING ####     #####     #####     #####     #####     #####     #####     #####     #####      #####
############################################################################################################
#### WARNING #### HERE BELOW YOU NEED TO CHECK AND MAKE THE FIGURE WITH THE SPECIFIC AXIS. #################
############################################################################################################
#### WARNING ####     #####     #####     #####     #####     #####     #####     #####     #####      #####

imageoutput=file.path("~/Documents/Solution/vcf/DNAsp/",paste(title,"D_GenomeWide.tiff", sep = ""))
tiff(file=imageoutput,height = 10, width = 23, units = 'cm',res = 300)
layout(matrix(c(1,1,1,1,1,2), nrow = 1, ncol = 6, byrow = TRUE))
par(oma = c(2,2,2,0),mar=c(4.1,4.1,1,0))
with(CROP, plot(pos,rep(0.1,nrow(CROP)),xlim=c(0,max(CROP[CROP$CHROM==chrlist[12], ]$pos)),ylim=c(-1.5,1.5),main="",bty="n",ylab= paste("Tajima's D mean",sep=""),xlab="Chromosomes",xaxt="n",yaxt="n",col="white",cex.lab=1.2))
axis(1, at=ticks, lab=chrlist2, cex.axis=1) #0.9
axis(2, at=c(-1.5,-1,-0.5,0,0.5,1,1.5), cex.axis=1,line = -1.5,las=1) #0.8

#############################
##    SIGNIFICANCE LINES   ##
#i in even chromosomes and j in uneven chromosmes...
 for (i in c(8))
 {
   rect(min(WILD[WILD$CHROM==chrlist[i],]$pos)-min(WILD[WILD$CHROM==chrlist[i],]$POS),1.4752,max(WILD[WILD$CHROM==chrlist[i],]$pos-1),1.49,col="dimgrey",border = NA)
 }
 for (j in c(9))
 {
   rect(min(WILD[WILD$CHROM==chrlist[j],]$pos)-min(WILD[WILD$CHROM==chrlist[j],]$POS),1.4752,max(WILD[WILD$CHROM==chrlist[j],]$pos-1),1.49,col="darkgrey",border = NA)
 }


################
## LINES WILD ##
for (i in (1:12))
{
  points(WILD[WILD$CHROM==chrlist[((2*i)-1)],]$pos,WILD[WILD$CHROM==chrlist[((2*i)-1)],]$D,col="darkgrey",lwd=1.3,pch=16,cex=.85)
} #pch=20
}
for (i in (1:12))
{
  points(WILD[WILD$CHROM==chrlist[2*i],]$pos,WILD[WILD$CHROM==chrlist[2*i],]$D,col="dimgrey",lwd=1.3,pch=16,cex=.85)
} #pch=20
}


################
## LINES CROP ##
for (i in (1:12))
{
  i = (2*i)-1
  points(CROP[CROP$CHROM==chrlist[i],]$pos,CROP[CROP$CHROM==chrlist[i],]$D,col=CropCol_Light,lwd=1.3,pch=16,cex=.85)
} #pch=20
for (i in (1:12))
{
  points(CROP[CROP$CHROM==chrlist[2*i],]$pos,CROP[CROP$CHROM==chrlist[2*i],]$D,col=CropCol_Dark,lwd=1.3,pch=16,cex=.85)
}

#######################
### D. DISTRIBUTION ###
#######################
par(mar=c(4.2,0,1,2.1))
plot(x=density(Rolling.wild)$y,y=density(Rolling.wild)$x, col="white",ylim=c(-1.5,1.5),ylab="",main = "",xlab="Density", xaxt="n",yaxt="n",cex.lab=1,col.axis="white",bty="n")
polygon(x=density(Rolling.wild)$y,y=density(Rolling.wild)$x, col=adjustcolor( "darkgrey", alpha.f = 0.8), border="dimgrey",lwd=1.3)
polygon(x=density(Rolling.crop)$y,y=density(Rolling.crop)$x, col=adjustcolor( CropCol_Light, alpha.f = 0.8), border=CropCol_Dark,lwd=1.3)

legend("topright", inset=0.05, c(pop2,pop3), fill=c(adjustcolor( CropCol_Light, alpha.f = 0.8),adjustcolor( "darkgrey", alpha.f = 0.8)), border=c(CropCol_Dark,"dimgrey"),bty="n" ,cex=1)

title(paste(Species,sep=""),cex.main=1.8, outer=TRUE)

dev.off()  

```






