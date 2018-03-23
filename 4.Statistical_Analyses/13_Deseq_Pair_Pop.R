---
title : 'Differencial Expression analyses in RNAseq'
author: "S. Arnoux using R script from Sauvage et al. 2016 \n C. Sauvage, A. Rau and J. Chadoeuf"
date : "January 2017"
output : html_document
---

###########################
### LOADING R LIBRARIES ###
###########################

```r  

#source("http://bioconductor.org/biocLite.R")
#install.packages("VennDiagram")
library("Rcpp")
library("DESeq2")
library("vsn")
library("RColorBrewer")
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
library(VennDiagram)
library(vsn)
library(hexbin)
library("pheatmap")
library("RColorBrewer")
library("genefilter")
library("PoiClaClu")
library(missMethyl) 
library(edgeR)

#####################
### LOAD THE DATA ###
#####################

rm(list=ls())

## Give the working folder  
PATH_STAT = "/path/to/Deseq_Folder/"  
FILE = "Exp_Coverage.txt"  
POP = "Acc_Loc.txt" 
#The follodwing texte file is from the step 3.g.Paralogous countig extraction
Para_Filter = "/path/to/Filtered_liste.tab.txt" 
setwd(PATH_STAT)

REF = "/path/to/reference/reference.fa"
GO_SLIM = "/path/to/reference/ref_interproSc_GO_pfam.txt" 
species = "Species_name"

## It will remove the columns with the Outgroups and keep only two populations
FILTER = c("Outgroup_Name")
LA_col.list <- c("firebrick2","dimgrey")
pch.list =c(18,20)
LA_C = "Population_nameA" 
LA_W = "Population_nameB"

TITLE_VENNDIAGRAM = "Differentially expressed genes"
subTITLE_VENNDIAGRAM ="Genes differentially regulated compare to Crop expressions"

## Number of gene [Top Variant] when it comes to pheatmap
NUMGENE = 600
pal=RColorBrewer::brewer.pal(3,"Dark2")

#####################
### RUN THE DESEQ ###
#####################

### First of all we start by reading the files and filtering the data
CountTable.cut = read.table(FILE, header=TRUE, row.names=1)
head(CountTable.cut)
Loc = read.table(POP, header=TRUE)
Loc$LocName = paste(Loc$Name,Loc$POP,sep="_") 
names(CountTable.cut)=Loc$LocName
head(CountTable.cut); dim(CountTable.cut)

condition.pop = Loc$POP

#Filter according to the SNPs and not-paralogs
parafilter = read.table(Para_Filter,header=FALSE)
para.CountTable.cut = CountTable.cut[rownames(CountTable.cut) %in% parafilter$V1,]
#sub.CountTable.cut = CountTable.cut
sub.CountTable.cut = subset(para.CountTable.cut, select =  grep(paste(FILTER,collapse="|"),names(para.CountTable.cut),invert=T))
head(sub.CountTable.cut); dim(sub.CountTable.cut)

condition.pop.cut = droplevels(Loc$POP[grep(paste(FILTER,collapse="|"),Loc$LocName,invert=T)])

condition.pop.cut #= Loc$POP

#para.CountTable.cut = subset(CountTable.cut, select =  grep(paste(FILTER,collapse="|"),names(CountTable.cut),invert=T))
RC.CountTable.cut = subset(CountTable.cut, select =  grep(paste(FILTER,collapse="|"),names(para.CountTable.cut),invert=T))
RC.CountTable.cut$tot = rowSums(RC.CountTable.cut)

## Gives you the number of expressed genes (RC) on the total genes possible (from the CDS)
RC = nrow(RC.CountTable.cut[!RC.CountTable.cut$tot==0,])
RC / nrow(CountTable.cut)

### Just to check the distribution of the raw read counts

par(mfrow=c(3,4))
for (i in (1:12))
{
    hist(log10(sub.CountTable.cut[,i]), col="grey", main = condition.pop.cut[i], xlab="raw read count")
}

### Then we plot the individuals using a pca based on gene expression data

cds.cut = DESeqDataSetFromMatrix(countData = sub.CountTable.cut, colData = data.frame(condition.pop.cut), design = ~condition.pop.cut)
cds.cut
cds.cut = DESeq(cds.cut, test = "Wald", fitType = c("parametric"), parallel = TRUE)

rld.cut= rlog(cds.cut, blind = TRUE)
PCA = plotPCA(rld.cut, intgroup=c("condition.pop.cut"), returnData = TRUE)

#### WARNING ####     #####     #####     #####     #####     #####     #####     #####     #####      #####
############################################################################################################
#### WARNING #### HERE BELOW YOU NEED TO CHECK AND MAKE THE FIGURE WITH THE SPECIFIC AXIS. #################
############################################################################################################
#### WARNING ####     #####     #####     #####     #####     #####     #####     #####     #####      #####

tiff(paste("~/Desktop/PCA_rld_all_",species,".tiff", sep = ""), width = 9, height = 9, units = "cm", res = 300)
par(mar=c(4.6,2,2,1))
plot(PCA$PC2, PCA$PC1, col=LA_col.list[as.integer(PCA$group)], xlab="", ylab="",pch=pch.list[as.integer(PCA$group)],xaxt ="n", yaxt="n",bty="n",cex=1.2,cex.lab=1.4)
axis(1, at=c(-50,-30,-10,10) ,cex.axis=0.9) #0.9
axis(2, at=c(-20,0,20,40), cex.axis=0.9) #0.8

legend("left", legend=c(LA_C,LA_W), pch=c(23,21), col=LA_col.list[1:2],bty="n",pt.cex=0.8,cex=0.8)
title(main= paste(species),lwd=2,lty=5,cex.main=1.2)
dev.off()

### HERE WE PERFORMED THE DIFFERENTIAL GENES EXPRESSION ANALYSIS 

vsd.cut = varianceStabilizingTransformation(cds.cut, blind=TRUE)
rlogMat <- assay(rld.cut)
vstMat <- assay(vsd.cut)


par(mai=ifelse(1:4 <= 2, par("mai"), 0))
px     <- counts(cds.cut)[,1] / sizeFactors(cds.cut)[1]
ord    <- order(px)
ord    <- ord[px[ord] < 150]
ord    <- ord[seq(1, length(ord), length=50)]
last   <- ord[length(ord)]
vstcol <- c("blue", "black")
tiff('Matplot_vsd_all.tiff', width = 2700, height = 1800, units = "px", res = 400)
matplot(px[ord], cbind(assay(vsd.cut)[, 1], log2(px))[ord, ], type="l", lty=1, col=vstcol, xlab="n", ylab="f(n)")
legend("bottomright", legend = c(expression("variance stabilizing transformation"), expression(log[2](n/s[1]))), fill=vstcol)
dev.off()

# Standard deviation over mean plot 
par(mfrow=c(1,3))
notAllZero <- (rowSums(counts(cds.cut))>0)
tiff(paste("SdPLot_log2_cds_all_",species,".tiff",sep = ""), width = 2700, height = 1800, units = "px", res = 400)
meanSdPlot(log2(counts(cds.cut,normalized=TRUE)[notAllZero,] + 1))
dev.off()
tiff(paste("SdPLot_rld_all_",species,".tiff", sep = ""), width = 2700, height = 1800, units = "px", res = 400)
meanSdPlot(assay(rld.cut[notAllZero,]))
dev.off()
tiff(paste("SdPLot_dsv_all_",species,".tiff", sep = ""), width = 2700, height = 1800, units = "px", res = 400)
meanSdPlot(assay(vsd.cut[notAllZero,]))
dev.off()

# Plot of data distribution and correlation between the two dds and rld
dds <- estimateSizeFactors(cds.cut)
tiff(paste("plot_dds_vs_rld_all_",species,".tiff",sep = ""), width = 2700, height = 1800, units = "px", res = 400)
par( mfrow = c( 1, 2 ) )
plot(log2( 1 + counts(dds, normalized=TRUE)[ , c(1,22)] ),
     pch=16, cex=0.3, main="estimateSizeFactors cds")
plot(assay(rld.cut)[ , c(1,22)],
     pch=16, cex=0.3, main="rlog cds")
dev.off()
sampleDists <- dist( t(rlogMat))
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )


# Heatmap and Matrix per gene
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
colnames(sampleDistMatrix) <- NULL
tiff(paste("pheatmap_rld_all_",species,".tiff",sep = ""), width = 2700, height = 1800, units = "px", res = 400)
pheatmap(sampleDistMatrix,
          clustering_distance_rows=sampleDists,
          clustering_distance_cols=sampleDists,
          col=colors)
dev.off()


# Heatmap and Matrix per gene with a poisson distribution adj.
poisd <- PoissonDistance(t(counts(cds.cut)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
colnames(samplePoisDistMatrix) <- NULL
row.names(samplePoisDistMatrix) <- row.names(sampleDistMatrix)
tiff(paste("pheatmap_poisson_cds_all_",species,".tiff",sep = ""), width = 2700, height = 1800, units = "px", res = 400)
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows=poisd$dd,
         clustering_distance_cols=poisd$dd,
         col=colors)
dev.off()

# Etimation of the genes with the greatest variance in expression
topVarGenes <- head(order(-rowVars(assay(rld.cut))),NUMGENE)
mat <- assay(rld.cut)[topVarGenes, ]
mat <-mat - rowMeans(mat)
df <- as.data.frame(colData(rld.cut)[,"condition.pop.cut"])
names(df) <- c("Group")
rownames(df) <- rownames(t(mat))
tiff(paste0("pheatmap_annotations_",NUMGENE,"_rld_all",species,".tiff",sep = ""), width = 2700, height = 1800, units = "px", res = 400)
pheatmap(mat, annotation_col=df, show_rownames = F,main = paste0("Diff.Var. of",NUMGENE,"genes in Tomato",sep = " "))
dev.off()

# Plot Dispersion # 
tiff(paste0("plotDispEsts_cds_all_",species,".tiff",sep = ""), width = 2700, height = 1800, units = "px", res = 400)
par(mfrow=c(1,1))
plotDispEsts(cds.cut)
dev.off()
tiff(paste0("plotDispEsts_rld_",species,".tiff",sep = ""), width = 2700, height = 1800, units = "px", res = 400)
par(mfrow=c(1,1))
plotDispEsts(rld.cut)
dev.off()

### HERE WE PERFORMED THE ANALYSES OF ALL OUR POPULATION 2 PER 2.

level=combn(levels(factor(condition.pop.cut)), 2)
x = level[,1]

res <- results(cds.cut, contrast = c(paste("condition.pop.cut"),as.vector(x)))
NameFile <- paste(paste("ResDESeq",x[1],x[2], sep="_"),".df", sep="")
write.table(res,NameFile)

dim <- dim(res)
resOrdered <- res[order(res$padj),]
Cond <- paste( "Pop", x[1], x[2], sep="_")
OUT<-paste( "Total.summary.", Cond, ".txt", sep="")
#head(resOrdered)
sum <- summary(res)
  
cond = condition.pop.cut
filt.data = sub.CountTable.cut
#filt.data = subset(CountTable.cut, rowSums(CountTable.cut)>5000)  
  
## Voom transformation
y <- DGEList(counts=filt.data)
y <- calcNormFactors(y)
des <- model.matrix(~0+cond)
v <- voom(y, des, plot=FALSE)$E

## Differential variability model
fitvar.contr <- varFit(y, design=des)
contr <- makeContrasts(condLyc. - condPim., levels=colnames(des))
fitvar.contr <- contrasts.varFit(fitvar.contr, contrasts=contr)
summary(decideTests(fitvar.contr))
topDV <- topVar(fitvar.contr)
topDV <- topDV[topDV$Adj.P.Value < 0.05,]
topDV

#We can also plot the voom-transformed counts corresponding to the 7 differentially variable genes.

par(mfrow=c(3,3), mar=c(4,4,2,2)) 
for(i in 1:9) {
  stripchart(v[rownames(v) == rownames(topDV)[i],] ~ as.factor(cond), method="jitter",
             group.names=c("crop", "Pim"), pch=16, cex=1.5, col=pal, ylab="Voom-counts",
             vertical=TRUE, cex.axis=1.5, cex.lab=1.5)
  title(rownames(topDV)[i], cex.main=1.5)
}


#We also examine whether these differentially variable genes were identified by the differential expression analysis. It is interesting to note that there are four novel genes identified with the DiffVar analysis!

diffvar_genes <- rownames(topDV)
diff_genes <- rownames(res[which(res$padj < 0.01),])
diffvar_genes %in% diff_genes
diffvar_genes[diffvar_genes %in% diff_genes == TRUE]

  
cat(paste(Cond), file = OUT, append = TRUE, sep = "\n")
capture.output(summary(res), file = OUT)

DISP <- (dispersions(cds.cut))
normalizationFactors(cds.cut)
estimateSizeFactors(cds.cut, locfunc=median)
norm.count <- counts(cds.cut, normalized=T)
#dim(norm.count); head(norm.count); class(norm.count)

norm.factors <- sizeFactors(cds.cut)

mcols(res)$description
resSig <- subset(resOrdered, padj < 0.01)
resSig

up <- subset(resSig, log2FoldChange > 0)
up
write.table(up[order(up$log2FoldChange, decreasing=TRUE),], paste(Cond,"UpReg.txt",sep=""))

down <- subset(resSig, log2FoldChange < 0)
down
write.table(down[order(down$log2FoldChange, decreasing=TRUE),], paste(Cond,"DownReg.txt",sep=""))

list.down = sort((rownames(down)), decreasing = FALSE)
write(list.down,
paste(Cond,"Up.txt",sep=""))
list.up = sort((rownames(up)), decreasing = FALSE)
write(list.up,
paste(Cond,"Down.txt",sep=""))

### MA PLOT ###
tiff(filename = paste("MAPlot_",Cond,"tiff",sep=""), width = 2700, height = 1800, units = "px", res = 400) 
par(mar=c(4.5, 4.5, 2, 0.5))
plotMA(res, main=paste("MA-plot of log2 fold changes",Cond,sep=" "), ylim=c(-8,8))
dev.off() 

### PLOT COUNTS ###
tiff(filename = paste("CountsPlot_",Cond,"tiff",sep=""), width = 2700, height = 1800, units = "px", res = 400)
pc <- plotCounts(cds.cut, gene=which.min(res$padj), intgroup="condition.pop.cut", returnData=TRUE)
ggplot(pc, aes(x=condition.pop.cut, y=count)) +
geom_point(position=position_jitter(w=0.1,h=0), aes(colour = condition.pop.cut),size=3) +     scale_y_log10(breaks=c(500,1000,2000,4000,8000))
dev.off() 

### THEN WE REPRESENT THE FOLD CHANGE VS PVALUES RELATIONSHIP USING VOLCANO PLOT ###
pvaladj = data.frame(res[6])
log2FC = data.frame(res[2])

tiff(filename = paste("Volcano_",Cond,"tiff",sep=""), width = 2700, height = 1800, units = "px", res = 400)
plot(log2FC[,1], -log10(pvaladj[,1]), xlab="log2(Fold Change)", ylab="-log10(adjusted pvalue)", main = paste("Volcano plot (Fold chang vs pvalues) for", Cond , sep = " "), pch = 16, col=ifelse(pvaladj[,1]<0.001, "#FF000050", "#00000050"))
abline(h=-log10(0.001), col="red")
dev.off() 

d = data.frame(resOrdered)
dim(d); head(d); summary(d)

## head(d, n=1)
d$ITAG = rownames(d)

d$ITAG <- as.character(d$ITAG) #just in case...sometime factors act wierd in matching statments
# head(d$ITAG)

# Need length of each ITAG, because goseq adjusts for this use Biostrings to calculate this

itagSeqs <- readDNAStringSet(file = REF) 
itagLength <- nchar(itagSeqs) #length of each ITAG
names(itagLength) <- names(itagSeqs)

# head(names(itagLength))
names(itagLength) <- substr(names(itagLength),1,18)

# Create GO term list in format needed for goseq. This file is on smart site at Maloof/Sinha Tomato Group Resources / Data / GOannotation

# This file is on pfam sites
GOinterpro_annex_slim <-  read.delim(GO_SLIM,row.names=NULL) #,as.is=T)
# head(GOinterpro_annex_slim)
summary(GOinterpro_annex_slim)


# Check to see if ITAG names are a problem.
sum(GOinterpro_annex_slim$ITAG %in% d$ITAG) #5568

#head(GOinterpro_annex_slim$ITAG)
head(d$ITAG)

# The number from the next command should match the number from the previous command.

sum(substr(GOinterpro_annex_slim$ITAG,1,18) %in% substr(d$ITAG,1,18)) #2407, good.

colnames(d)
d=na.omit(d) #19332 #19628


 ##### #####  ##### ##### ##### ##### ##### ##### ##### ##### ##### 
 ##### ##### ##### WE PERFORM THE ENRICHMENT TEST ##### ##### ##### 
 ##### #####  ##### ##### ##### ##### ##### ##### ##### ##### ##### 
#### You might want to change the threshold below

GO.sets <- ls(pattern="GO[[:alnum:]]")
gene.names=d$ITAG
FC=d$log2FoldChange
pval=d$padj
FC.thresh=0
p.thresh=.01
ilength=itagLength
go.cutoff=.05
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
up <- as.integer(FC > FC.thresh & pval < p.thresh) #upregulated genes
names(up) <- gene.names
down <- as.integer(FC < - FC.thresh & pval < p.thresh) #downregulated genes
names(down) <- gene.names

print(summary(up))
print(summary(down))

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

results.up = up.go[up.go$over<go.cutoff, c(1,2,4,5,6,7)]
dim(results.up); head(results.up)
results.down = down.go[down.go$over<go.cutoff, c(1,2,4,5,6,7)]
dim(results.down); head(results.down)
length(up[up==1])  #3116
length(down[down==1]) #4034

OUTGO.DOWN<-paste(Cond,"GO.fdr.wall.DOWN.summary.txt",sep="")
OUTGO.UP<-paste(Cond,"GO.fdr.wall.UP.summary.txt",sep="")

write.table(results.up, OUTGO.UP, row.names = FALSE, col.names = TRUE) 

write.table(results.down, OUTGO.DOWN, row.names = FALSE, col.names = TRUE) 


fDown <- list.files(path = PATH_STAT, pattern = "Down.txt$", full.names=TRUE)
fUp <- list.files(path = PATH_STAT , pattern = "Up.txt$", full.names=TRUE)

myfilelistD <- lapply(fDown, read.table)
myfilelistU <- lapply(fUp, read.table)

names(myfilelistD) <- list.files(path = PATH_STAT, pattern = "Down.txt$", full.names=FALSE)
names(myfilelistU) <- list.files(path = PATH_STAT, pattern = "Up.txt$", full.names=FALSE)

Number_Down=NULL
Number_Up=NULL
for(x in 1:length(myfilelistD)){
  Number_Down = rbind(Number_Down,as.vector(c(gsub(".Down.txt","",names(myfilelistD)[x]), nrow(myfilelistD[[x]]))))
}
for(x in 1:length(myfilelistU)){
  Number_Up = rbind(Number_Up,as.vector(c(gsub(".Up.txt","",names(myfilelistU)[x]), nrow(myfilelistU[[x]]))))
}
DEG_Num = as.data.frame(cbind(Accesions = Number_Up[,1],Up_Reg = Number_Up[,2],Down_Reg = Number_Down[,2]))

write.table(DEG_Num, "Summary_DEG_Num.txt", row.names = FALSE, col.names = TRUE)

```

