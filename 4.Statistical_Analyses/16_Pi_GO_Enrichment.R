---
title : 'Gene ontology enrichement analysis on Pi shifted genes'
author: "S. Arnoux"
date : "January 2018"
output : html_document
---

###############################################
### STARTING by giving the path to the data ###
###############################################

```r

# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
library("zoo")
library("grDevices")
library(fields)
library(plyr) 
library(data.table)
library(missMethyl)
library(edgeR)
library("DESeq2")
library("Rcpp")
library("vsn")
library("gplots")
library("airway")
library("pasilla")
library("Biobase")
library("BiocParallel")
library("ggplot2")
library("reshape")
register(MulticoreParam(3))
library(goseq)
library(GO.db)
library(ape)
library(Biostrings)


```

#####################
### LOAD THE DATA ###
#####################

```
rm(list=ls())
## Give the working folder
PATH_STAT = "/path/to/DNAsp_results/"
REF = "/path/to/reference/reference.fa"
GO_SLIM = "/path/to/reference/ref_interproSc_GO_pfam.txt" 
species = "Species_Name"
Colspecies = "#A164D0"
ChrLOC=read.table("/path/to/reference/ref.Gene_Loc.tab")
setwd(PATH_STAT)

## Give sample details 
pop2 = "Pop_A_name" 
pop3 = "Pop_B_name"

Species= "species_name"
title = "Prefix_to_the_output_files"

## Colours for the domesticated population // The wild will be in grey
CropCol_Light = "#A164D0"
CropCol_Dark = "#5D04A1"

## Settle the file names
Name_File2 = "DNAsp_results_PopA.out" 
Name_File3 = "DNAsp_results_PopB.out"


## Reads the files
stat2 =  read.table(Name_File2, header = TRUE)
colnames(stat2)[1] <- "Gene"
stat3 =  read.table(Name_File3, header = TRUE)
colnames(stat3)[1] <- "Gene"


### Data for Coverage 
Crop_Wild_Up = "/path/to/Deseq_Folder/Pop_A_B_DownReg.txt"
Crop_Wild_Down = "/path/to/Deseq_Folder/Pop_A_B_UpReg.txt"
CropUp_reg = read.table(Crop_Wild_Down, header = TRUE)
CropDown_reg = read.table(Crop_Wild_Up, header = TRUE)

```
#####################################
### Sanity Check and data reading ###
#####################################

```
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
#113 badly annotated ...
FdtCrop <- pFdtCrop[!Gene%in%BadGENES$Gene,]
FdtCrop = FdtCrop[order(FdtCrop$V2),]
FdtCrop = FdtCrop[order(FdtCrop$chr),] 
FdtWild <- pFdtWild[!Gene%in%BadGENES$Gene,]
FdtWild = FdtWild[order(FdtWild$V2),]
FdtWild = FdtWild[order(FdtWild$chr),] 


FdtWild$DEG = 0
FdtWild$DEG <- ifelse(FdtWild$Gene %in% rownames(CropUp_reg), 1, 0 )
FdtWild$DEG <- ifelse(FdtWild$Gene %in% rownames(CropDown_reg), -1, FdtWild$DEG)
FdtCrop$DEG = 0
FdtCrop$DEG <- ifelse(FdtCrop$Gene %in% rownames(CropUp_reg), 1, 0 )
FdtCrop$DEG <- ifelse(FdtCrop$Gene %in% rownames(CropDown_reg), -1, FdtCrop$DEG)


```
########################
### Delta Pi / Taj D ###
########################

```

delta.pi.dist=data.frame(as.character(FdtWild$Gene),as.numeric(FdtCrop$chr), as.numeric(FdtCrop$V2), as.numeric(FdtCrop$DEG), as.numeric(FdtWild$Pi),as.numeric(FdtCrop$Pi))
colnames(delta.pi.dist) = c("GENE","CHROM","POS", "DEG", "PI_Wild","PI_Crop"); dim(delta.pi.dist); head(delta.pi.dist)
delta.pi.dist = na.omit(delta.pi.dist)


```
#######################
### GET THE EVAL GO ###
#######################

```

d = delta.pi.dist
d$ITAG <- as.character(d$GENE) #just in case...sometime factors act wierd in  

# Need length of each ITAG, because goseq adjusts for this use Biostrings to calculate this

itagSeqs <- readDNAStringSet(file = REF) 
itagLength <- nchar(itagSeqs) #length of each ITAG
names(itagLength) <- names(itagSeqs)
names(itagLength) <- substr(names(itagLength),1,20) # Adjust if your annotation is too long
d$ITAG[!d$ITAG %in% names(itagLength)] #Looks OK

# This file is on smart site at Maloof/Sinha Tomato Group Resources / Data / GOannotation 
GOinterpro_annex_slim <-  read.delim(GO_SLIM,row.names=NULL) #,as.is=T)
summary(GOinterpro_annex_slim)


# Check to see if ITAG names are a problem.

sum(GOinterpro_annex_slim$ITAG %in% d$ITAG) #32

#head(GOinterpro_annex_slim$ITAG)
head(d$ITAG)

# The number from the next command should match the number from the previous command.

sum(substr(GOinterpro_annex_slim$ITAG,1,20) %in% substr(d$ITAG,1,20)) #2407, good.

# Wrapper function to actually do the GO evalutation. Change the default arguments to match your file. PiC = fold change info. PiW = adjusted DE PiWues. Check you d file to make sure the names match.

colnames(d)

##Filter
maxPiWild = quantile(d$PI_Wild, 1 - 1e-5,na.rm=TRUE)[[1]]
maxPiCrop = quantile(d$PI_Crop, 1 - 1e-5,na.rm=TRUE)[[1]]

which((d$PI_Crop) >= maxPiCrop)
# 11220 14263 15632
which((d$PI_Wild) >= maxPiWild)
# [1]  2523  9125  9966 11944 12120 12515 17304

##Filter
d$PI_Crop[(d$PI_Crop) >= maxPiCrop] <- NA
d$PI_Wild[(d$PI_Wild) >= maxPiWild] <- NA
d=na.omit(d)
##### WE CREATE AN EVALGO FUNCTION THAT PERFORM THE ENRICHMENT TEST #####

GO.sets <- ls(pattern="GO[[:alnum:]]")
gene.names=d$GENE
PiC=d$PI_Crop
PiW=d$PI_Wild
PiC.thresh=quantile(d$PI_Crop, 0.95,na.rm=TRUE)[[1]]
PiC.threshmin=0.001*max(PiC)
PiW.thresh=quantile(d$PI_Wild, 0.95,na.rm=TRUE)[[1]]
PiW.threshmin=0.001*max(PiC)/max(PiW)*max(PiW)
ilength=itagLength
go.cutoff=.1
keep.GO="BP"
type="GO"
go.terms=get(GO.sets[1])

#add GO: header if needed

head(go.terms)
if (type=="GO" & length(grep("GO",go.terms$GO[1]))==0) {
  go.terms$GO <- gsub("([0-9]{7})","GO:\\1",go.terms$GO)
}

#remove extra spaces
go.terms$GO <- gsub(" +","",go.terms$GO)

#get length list to match gene names
ilength <- ilength[names(ilength) %in% gene.names]

#filter go terms to match gene list
go.terms <- go.terms[go.terms$ITAG %in% gene.names,]
#head(go.terms)

#convert go terms to list
go.list <- strsplit(as.character(go.terms$GO),split=",")
head(go.list)
names(go.list) <- go.terms$ITAG

#filter genes based on criterion
up <- as.integer(PiC > PiC.thresh & PiW < PiW.threshmin) 
# up <- as.integer(PiC > PiC.thresh & PiW < PiW.thresh & PiW > PiW.threshmin)
names(up) <- gene.names #upSelected genes

down <- as.integer(PiW > PiW.thresh & PiC < PiC.threshmin) 
# down <- as.integer(PiW > PiW.thresh & PiC < PiC.thresh & PiC > PiC.threshmin) 
names(down) <- gene.names #downregulated genes

#calculate bias function
up.pwf <- nullp(up,bias.data=ilength,plot.fit=F)
down.pwf <- nullp(down,bias.data=ilength,plot.fit=F)

#calculate p-values for over-representation
up.go <- goseq(up.pwf, gene2cat=go.list, use_genes_without_cat=TRUE, method="Wallenius")
down.go <- goseq(down.pwf, gene2cat=go.list, use_genes_without_cat=TRUE, method="Wallenius")

if (type=="GO") { #add GO term description
  up.go$description <- Term(up.go$category)
  up.go$ontology <- Ontology(up.go$category)
  down.go$description <- Term(down.go$category)
  down.go$ontology <- Ontology(down.go$category)
  
  #filter for GO categories of interest
  up.go <- up.go[up.go$ontology==keep.GO,]
  down.go <- down.go[down.go$ontology==keep.GO,]
  
  #remove NAs
  up.go <- up.go[!is.na(up.go$ontology),]
  down.go <- down.go[!is.na(down.go$ontology),]
}   

if (type=="mapman") {#add mapman description
  up.go <- merge(up.go,bincodes,by.x="category",by.y="BINCODE",sort=F)
  down.go <- merge(down.go,bincodes,by.x="category",by.y="BINCODE",sort=F)
}

#truncate to go.cutoff threshold
up.go.res <- up.go[up.go$over<go.cutoff,]
down.go.res <- down.go[down.go$over<go.cutoff,]

results.up = up.go[up.go$over<go.cutoff, c(1,2,4,5,6,7)]
dim(results.up); head(results.up)
results.down = down.go[down.go$over<go.cutoff, c(1,2,4,5,6,7)]
dim(results.down); head(results.down)
length(up[up==1])  #A 196 188
length(down[down==1]) #B 374 370
########
##### LIST OF A and B Genes !!!
########
OUTGO.A<-paste(title,"_Pi_A_CROP_list.txt",sep="")
write.table(names(up[up==1]), OUTGO.A, row.names = FALSE, col.names = "GENE")
OUTGO.B<-paste(title,"_Pi_B_CROP_list.txt",sep="")
write.table(names(down[down==1]), OUTGO.B, row.names = FALSE, col.names = "GENE")

########
##### LIST OF GOTerms !!!
########
OUTGO.UP<-paste(title,"_Pi_CROP_GO.filt.wall.A.summary.txt",sep="")
write.table(results.up, OUTGO.UP, row.names = FALSE, col.names = TRUE)
OUTGO.DOWN<-paste(title,"_Pi_CROP_GO.filt.wall.B.summary.txt",sep="")
write.table(results.down, OUTGO.DOWN, row.names = FALSE, col.names = TRUE)

#### WARNING ####     #####     #####     #####     #####     #####     #####     #####     #####      #####
############################################################################################################
#### WARNING #### HERE BELOW YOU NEED TO CHECK AND MAKE THE FIGURE WITH THE SPECIFIC AXIS. #################
############################################################################################################
#### WARNING ####     #####     #####     #####     #####     #####     #####     #####     #####      #####

imageoutput=file.path("~/Documents/Solution/vcf/DNAsp/",paste(title,"Pi_GO_filt.wall_zoom.tiff", sep = ""))
tiff(file=imageoutput,height = 10, width = 10, units = 'cm',res = 300)
par(mar=c(4.2,5.1,3.1,2.1))
plot(d$PI_Crop~d$PI_Wild, xlim=c(0,0.03),xaxt="n", ylim=c(0,0.03),yaxt="n", cex.axis=1,pch=20, cex.lab=1.4,xlab =expression(paste(pi," Wild",sep="")),ylab = expression(paste(pi," Crop",sep="")), col = ifelse((d$PI_Crop > PiC.thresh & d$PI_Wild < PiW.threshmin | d$PI_Wild > PiW.thresh & d$PI_Crop < PiC.threshmin ),Colspecies,'dimgray'),las=1)
#xlim=c(0,0.06),  ylim=c(0,0.06), yaxt="n", xaxt="n",
axis(1, at=c(0,0.02,0.04,0.05), cex.axis=0.8,las=1) #xaxis
axis(2, at=c(0,0.02,0.04,0.05), cex.axis=0.8,las=1) #0.8
axis(1, at=c(0,0.05,0.1,0.15), cex.axis=0.8,las=1) #xaxis
axis(2, at=c(0,0.05,0.1,0.15), cex.axis=0.8,las=1) #0.8

abline(h=PiC.threshmin,col="gray20",lwd=1.2,lty=5)
abline(h=PiC.thresh,col="gray20",lwd=1.2,lty=5)
abline(v=PiW.threshmin,col="gray20",lwd=1.2,lty=5)
abline(v=PiW.thresh,col="gray20",lwd=1.2,lty=5)

title(main= paste(Species),lwd=1.8,lty=5,cex.main=1.4)

dev.off()

```
