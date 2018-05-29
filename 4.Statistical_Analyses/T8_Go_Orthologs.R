###############################################
### STARTING by giving the path to the data ###
###############################################

library("zoo")
library("grDevices")
library(fields)
library(plyr) 
library(data.table)
library(missMethyl)
library(edgeR)
library("DESeq2")
library("Rcpp")
# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
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
PATH_STAT = "/Users/stephaniearnoux/Documents/Solution/vcf/DNAsp/"
REF = "~/Documents/Solution/Transcriptomics/Coverage/ITAG3.2_PerGene/ITAG3.2_CDS.fasta"
GO_SLIM = "~/Documents/Solution/Transcriptomics/Coverage/ITAG3.2_PerGene/ITAG3.2_InterproSc_GO_pfam.txt" #SMART

ChrLOC=read.table("/Users/stephaniearnoux/Documents/Solution/references/ITAG3.2_gene_Loc.tab")
setwd(PATH_STAT)
## Give sample details 

Species= "Orthologs"
title = "Orthologs_"

## Settle the file names
Name_File2 = "LA_DNAsp_Crop.recode.VCF.out" 
Name_File3 = "LA_DNAsp_Per.recode.VCF.out"

################### NOW DONT TOUCH BUT RUN IT THROUGH
## Reads the files
ListeTotalOrtho = read.table("/Users/stephaniearnoux/Documents/Solution/references/Ortholog_Analizes/2018_Tot_ortho.proteinortho")
ListeSolycOrtho = ListeTotalOrtho[2]
colnames(ListeSolycOrtho) <- "GENE"

Wild_Crop_Down = "~/Documents/Solution/Orthologs/List_Ortho_DOWN_TOT_Solyc.txt"
Wild_Crop_Up = "~/Documents/Solution/Orthologs/List_Ortho_Up_TOT_Solyc.txt"
Wild_Crop_A = "~/Documents/Solution/Orthologs/List_Ortho_A_TOT_Solyc.txt"
Wild_Crop_B = "~/Documents/Solution/Orthologs/List_Ortho_B_TOT_Solyc.txt"

CropUp_reg = read.table(Wild_Crop_Up, header = TRUE)
colnames(CropUp_reg) <- "GENE"
CropDown_reg = read.table(Wild_Crop_Down, header = TRUE)
colnames(CropDown_reg) <- "GENE"
CropA_reg = read.table(Wild_Crop_A, header = TRUE)
colnames(CropA_reg) <- "GENE"
CropB_reg = read.table(Wild_Crop_B, header = TRUE)
colnames(CropB_reg) <- "GENE"

```
#######################
### GET THE EVAL GO ###
#######################

```
  
d = ListeSolycOrtho
d$ITAG <- as.character(d$GENE) #just in case...sometime factors act wierd in  

# Need length of each ITAG, because goseq adjusts for this use Biostrings to calculate this

itagSeqs <- readDNAStringSet(file = REF) 
itagLength <- nchar(itagSeqs) #length of each ITAG
names(itagLength) <- names(itagSeqs)

# head(names(itagLength))
names(itagLength) <- substr(names(itagLength),1,18) # not needed if you use seqinR
# head(d$ITAG)
d$ITAG[!d$ITAG %in% names(itagLength)] #Looks OK

# This file is on smart site at Maloof/Sinha Tomato Group Resources / Data / GOannotation 
GOinterpro_annex_slim <-  read.delim(GO_SLIM,row.names=NULL) #,as.is=T)
# head(GOinterpro_annex_slim) 
summary(GOinterpro_annex_slim)


# Check to see if ITAG names are a problem.

sum(GOinterpro_annex_slim$ITAG %in% d$ITAG) #32

#head(GOinterpro_annex_slim$ITAG)
head(d$ITAG)

# The number from the next command should match the number from the previous command.

sum(substr(GOinterpro_annex_slim$ITAG,1,20) %in% substr(d$ITAG,1,20)) # 8222
colnames(d)

##### WE CREATE AN EVALGO FUNCTION THAT PERFORM THE ENRICHMENT TEST #####
GO.sets <- ls(pattern="GO[[:alnum:]]")
gene.names=d$ITAG
#p.thresh=.01
ilength=itagLength
go.cutoff=.05
keep.GO="BP"
type="GO"
go.terms=get(GO.sets[1])

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
  
up <- as.integer(d$ITAG %in% CropUp_reg$GENE)
names(up) <- gene.names
down <- as.integer(d$ITAG %in% CropDown_reg$GENE)
names(down) <- gene.names
AA <- as.integer(d$ITAG %in% CropA_reg$GENE)
names(AA) <- gene.names
BB <- as.integer(d$ITAG %in% CropB_reg$GENE)
names(BB) <- gene.names

up.pwf <- nullp(up,bias.data=ilength,plot.fit=F)
down.pwf <- nullp(down,bias.data=ilength,plot.fit=F)
AA.pwf <- nullp(AA,bias.data=ilength,plot.fit=F)
BB.pwf <- nullp(BB,bias.data=ilength,plot.fit=F)

#calculate p-values for over-representation
up.go <- goseq(up.pwf, gene2cat=go.list, use_genes_without_cat=TRUE, method="Wallenius")
down.go <- goseq(down.pwf, gene2cat=go.list, use_genes_without_cat=TRUE, method="Wallenius")
AA.go <- goseq(AA.pwf, gene2cat=go.list, use_genes_without_cat=TRUE, method="Wallenius")
BB.go <- goseq(BB.pwf, gene2cat=go.list, use_genes_without_cat=TRUE, method="Wallenius")

if (type=="GO") { #add GO term description
    up.go$description <- Term(up.go$category)
    up.go$ontology <- Ontology(up.go$category)
    down.go$description <- Term(down.go$category)
    down.go$ontology <- Ontology(down.go$category)
    AA.go$description <- Term(AA.go$category)
    AA.go$ontology <- Ontology(AA.go$category)
    BB.go$description <- Term(BB.go$category)
    BB.go$ontology <- Ontology(BB.go$category)
    #filter for GO categories of interest
    up.go <- up.go[up.go$ontology==keep.GO,]
    down.go <- down.go[down.go$ontology==keep.GO,]
    AA.go <- AA.go[AA.go$ontology==keep.GO,]
    BB.go <- BB.go[BB.go$ontology==keep.GO,]
    #remove NAs
    up.go <- up.go[!is.na(up.go$ontology),]
    down.go <- down.go[!is.na(down.go$ontology),]
    AA.go <- AA.go[!is.na(AA.go$ontology),]
    BB.go <- BB.go[!is.na(BB.go$ontology),]
}  

if (type=="mapman") {#add mapman description
      up.go <- merge(up.go,bincodes,by.x="category",by.y="BINCODE",sort=F)
      down.go <- merge(down.go,bincodes,by.x="category",by.y="BINCODE",sort=F)
      AA.go <- merge(AA.go,bincodes,by.x="category",by.y="BINCODE",sort=F)
      BB.go <- merge(BB.go,bincodes,by.x="category",by.y="BINCODE",sort=F)
      }

results.up = up.go[up.go$over<go.cutoff, c(1,2,4,5,6,7)]
dim(results.up); head(results.up)
results.down = down.go[down.go$over<go.cutoff, c(1,2,4,5,6,7)]
dim(results.down); head(results.down)
results.AA = AA.go[AA.go$over<go.cutoff, c(1,2,4,5,6,7)]
dim(results.AA); head(results.AA)
results.BB = BB.go[BB.go$over<go.cutoff, c(1,2,4,5,6,7)]
dim(results.BB); head(results.BB)

length(up[up==1])  #90 #406 with GO
length(down[down==1])  #48 #470 with GO
length(AA[AA==1])  #1 #ALL:41
length(BB[BB==1])  #17 #ALL:137

########
#### SUM1MARY GO TERMS !!!
########
OUTGO.UP<-paste("/Users/stephaniearnoux/Documents/Solution/Orthologs/",title,"_CROP_GO.filt.wall.UP.TOT.summary.txt",sep="")
write.table(results.up, OUTGO.UP, row.names = FALSE, col.names = TRUE)
OUTGO.DOWN<-paste("/Users/stephaniearnoux/Documents/Solution/Orthologs/",title,"_CROP_GO.filt.wall.DOWN.TOT.summary.txt",sep="")
write.table(results.down, OUTGO.DOWN, row.names = FALSE, col.names = TRUE)
OUTGO.A<-paste("/Users/stephaniearnoux/Documents/Solution/Orthologs/",title,"_CROP_GO.filt.wall.A.TOT.summary.txt",sep="")
write.table(results.AA, OUTGO.A, row.names = FALSE, col.names = TRUE)
OUTGO.B<-paste("/Users/stephaniearnoux/Documents/Solution/Orthologs/",title,"_CROP_GO.filt.wall.B.TOT.summary.txt",sep="")
write.table(results.BB, OUTGO.B, row.names = FALSE, col.names = TRUE)


