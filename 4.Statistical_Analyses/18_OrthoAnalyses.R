---
title: "Othology_Analyses"
author: "Stéphanie Arnoux"
date: "5/27/2018"
output: html_document
---

##################################################
### LOADING Results from proteinortho software ###
##################################################

```r
rm(list=ls())

PATH_to_LIST = "/path/to/folder_with_Orthologs/"
setwd(PATH_to_LIST)

ListeTotalOrtho = read.table("/path/to/folder_with_Orthologs/Common_Ortho.proteinortho")

#####################################################
#################### DEG ANALYSES ################### 
#####################################################
PATH_to_Cov = "/path/to/Deseq_Folder/"
####### GET THE DATA -> Change Tomato, eggplant and pepper with your species 1, 2 and 3
Tom = "ITAG3.2_PerGene/DESeq_FilteredPara_vbeta/Lyc_Pim/Pop_Lyc._Pim.Down.txt"
Egg = "SMEL_PerGene/DESeq_FilteredPara_vbeta/mel_ins/Pop_ins_melUp.txt"
Pep = "CaZl1_PerGene/DESeq_FilteredPara_vbeta/Crop_Wild_dadi/Pop_Crop_WildDown.txt"
TomD = "ITAG3.2_PerGene/DESeq_FilteredPara_vbeta/Lyc_Pim/Pop_Lyc._Pim.Up.txt"
EggD = "SMEL_PerGene/DESeq_FilteredPara_vbeta/mel_ins/Pop_ins_melDown.txt"
PepD = "CaZl1_PerGene/DESeq_FilteredPara_vbeta/Crop_Wild_dadi/Pop_Crop_WildUp.txt"
#PepDV = "pepper/Ortho_Pop_Crop_WildBOTH_list.df"
Tom.DUp = read.table(paste(PATH_to_Cov,Tom, sep=""), header = TRUE)
Egg.DUp = read.table(paste(PATH_to_Cov,Egg, sep=""), header = TRUE)
Pep.DUp = read.table(paste(PATH_to_Cov,Pep, sep=""), header = TRUE)
Tom.Down = read.table(paste(PATH_to_Cov,TomD, sep=""), header = TRUE)
Egg.Down = read.table(paste(PATH_to_Cov,EggD, sep=""), header = TRUE)
Pep.Down = read.table(paste(PATH_to_Cov,PepD, sep=""), header = TRUE)

####### MAKE COMPARISONS
listSoly.Tom.DUp = ListeTotalOrtho[which(as.character(Tom.DUp[,1]) %in% as.character(ListeTotalOrtho$V2) == TRUE),]
listSoly.Egg.DUp = ListeTotalOrtho[which(as.character(Egg.DUp[,1]) %in% as.character(ListeTotalOrtho$V3) == TRUE),]
listSoly.Pep.DUp = ListeTotalOrtho[which(as.character(Pep.DUp[,1]) %in% as.character(ListeTotalOrtho$V1) == TRUE),]
DUpPepEgg = merge(listSoly.Egg.DUp,listSoly.Pep.DUp, by = "V3")
nrow(DUpPepEgg)
write.table(DUpPepEgg[,3],file = "List_Ortho_Up_MM_PM_Solyc.txt",row.names = FALSE, col.names = "Genes_Name")

DUpPepTom = merge(listSoly.Tom.DUp,listSoly.Pep.DUp, by = "V3")
nrow(DUpPepTom)
write.table(DUpPepTom[,3],file = "List_Ortho_Up_LA_PM_Solyc.txt",row.names = FALSE, col.names = "Genes_Name")

DUpEggTom = merge(listSoly.Tom.DUp,listSoly.Egg.DUp, by = "V3")
nrow(DUpEggTom)
write.table(DUpEggTom[,3],file = "List_Ortho_Up_LA_MM_Solyc.txt",row.names = FALSE, col.names = "Genes_Name")

DUpShared = merge(DUpPepEgg, listSoly.Tom.DUp, by = "V3")
nrow(DUpShared)
write.table(DUpShared[,3],file = "List_Ortho_Up_LA_MM_PM_Solyc.txt",row.names = FALSE, col.names = "Genes_Name")

UpTOT = c(as.character(DUpShared[,3]),as.character(DUpPepEgg[,3]), as.character(DUpEggTom[,3]), as.character(DUpPepTom[,3]))
write.table(UpTOT,file = "List_Ortho_UP_TOT_Solyc.txt",row.names = FALSE, col.names = "Genes_Name")

listSoly.Egg.Down = ListeTotalOrtho[which(as.character(Egg.Down[,1]) %in% as.character(ListeTotalOrtho$V3) == TRUE),]
listSoly.Pep.Down = ListeTotalOrtho[which(as.character(Pep.Down[,1]) %in% as.character(ListeTotalOrtho$V1) == TRUE),]
listSoly.Tom.Down = ListeTotalOrtho[which(as.character(Tom.Down[,1]) %in% as.character(ListeTotalOrtho$V2) == TRUE),]
DownPepEgg = merge(listSoly.Egg.Down,listSoly.Pep.Down, by = "V1")
nrow(DownPepEgg)
write.table(DownPepEgg[,2],file = "List_Ortho_DOWN_MM_PM_Solyc.txt",row.names = FALSE, col.names = "Genes_Name")

DownTomPep = merge(listSoly.Tom.Down,listSoly.Pep.Down, by = "V1")
nrow(DownTomPep)
write.table(DownTomPep[,2],file = "List_Ortho_DOWN_LA_PM_Solyc.txt",row.names = FALSE, col.names = "Genes_Name")

DownTomEgg = merge(listSoly.Tom.Down,listSoly.Egg.Down, by = "V1")
nrow(DownTomEgg)
write.table(DownTomEgg[,2],file = "List_Ortho_DOWN_LA_MM_Solyc.txt",row.names = FALSE, col.names = "Genes_Name")

DownShared = merge(DownPepEgg, listSoly.Tom.Down, by = "V1")
nrow(DownShared)
write.table(DownShared[,2],file = "List_Ortho_DOWN_LA_MM_PM_Solyc.txt",row.names = FALSE, col.names = "Genes_Name")

TOT = c(as.character(DownShared[,2]),as.character(DownPepEgg[,2]), as.character(DownTomPep[,2]), as.character(DownTomEgg[,2]))
dim(TOT)
write.table(TOT,file = "List_Ortho_DOWN_TOT_Solyc.txt",row.names = FALSE, col.names = "Genes_Name")

###############################################
#################### A & B #################### 
###############################################
PATH_to_Pi = "/path/to/DNAsp/"
####### GET THE DATE
TomA = "LA_ALL_DNAsp_Pi_A_CROP_list.txt"
EggA = "MM_DNAsp_Pi_A_CROP_list.txt"
PepA = "PM_DNAsp_Pi_A_CROP_list.txt"
TomB = "LA_ALL_DNAsp_Pi_B_CROP_list.txt"
EggB = "MM_DNAsp_Pi_B_CROP_list.txt"
PepB = "PM_DNAsp_Pi_B_CROP_list.txt"

Tom.A = read.table(paste(PATH_to_Pi,TomA, sep=""), header = TRUE)
Egg.A = read.table(paste(PATH_to_Pi,EggA, sep=""), header = TRUE)
Pep.A = read.table(paste(PATH_to_Pi,PepA, sep=""), header = TRUE)
Tom.B = read.table(paste(PATH_to_Pi,TomB, sep=""), header = TRUE)
Egg.B = read.table(paste(PATH_to_Pi,EggB, sep=""), header = TRUE)
Pep.B = read.table(paste(PATH_to_Pi,PepB, sep=""), header = TRUE)

####### MAKE COMPARISONS
listSoly.Tom.A = ListeTotalOrtho[which(as.character(Tom.A[,1]) %in% as.character(ListeTotalOrtho$V2) == TRUE),]
listSoly.Egg.A = ListeTotalOrtho[which(as.character(Egg.A[,1]) %in% as.character(ListeTotalOrtho$V3) == TRUE),]
listSoly.Pep.A = ListeTotalOrtho[which(as.character(Pep.A[,1]) %in% as.character(ListeTotalOrtho$V1) == TRUE),]
APepEgg = merge(listSoly.Egg.A,listSoly.Pep.A, by = "V3")
nrow(APepEgg)
write.table(APepEgg[,3],file = "List_Ortho_A_MM_PM_Solyc.txt",row.names = FALSE, col.names = "Genes_Name")

APepTom = merge(listSoly.Tom.A,listSoly.Pep.A, by = "V3")
nrow(APepTom)
write.table(APepTom[,3],file = "List_Ortho_A_LA_PM_Solyc.txt",row.names = FALSE, col.names = "Genes_Name")

AEggTom = merge(listSoly.Tom.A,listSoly.Egg.A, by = "V3")
nrow(AEggTom)
write.table(AEggTom[,3],file = "List_Ortho_A_LA_MM_Solyc.txt",row.names = FALSE, col.names = "Genes_Name")

AShared = merge(APepEgg, listSoly.Tom.A, by = "V3")
nrow(AShared)
write.table(AShared[,3],file = "List_Ortho_A_LA_MM_PM_Solyc.txt",row.names = FALSE, col.names = "Genes_Name")

ATOT = c(as.character(AShared[,3]),as.character(AEggTom[,3]), as.character(APepTom[,3]), as.character(APepEgg[,3]))
write.table(ATOT,file = "List_Ortho_A_TOT_Solyc.txt",row.names = FALSE, col.names = "Genes_Name")


listSoly.Egg.B = ListeTotalOrtho[which(as.character(Egg.B[,1]) %in% as.character(ListeTotalOrtho$V3) == TRUE),]
listSoly.Pep.B = ListeTotalOrtho[which(as.character(Pep.B[,1]) %in% as.character(ListeTotalOrtho$V1) == TRUE),]
listSoly.Tom.B = ListeTotalOrtho[which(as.character(Tom.B[,1]) %in% as.character(ListeTotalOrtho$V2) == TRUE),]
BPepEgg = merge(listSoly.Egg.B,listSoly.Pep.B, by = "V1")
nrow(BPepEgg)
write.table(BPepEgg[,2],file = "List_Ortho_B_MM_PM_Solyc.txt",row.names = FALSE, col.names = "Genes_Name")

BTomPep = merge(listSoly.Tom.B,listSoly.Pep.B, by = "V1")
nrow(BTomPep)
write.table(BTomPep[,2],file = "List_Ortho_B_LA_PM_Solyc.txt",row.names = FALSE, col.names = "Genes_Name")

BTomEgg = merge(listSoly.Tom.B,listSoly.Egg.B, by = "V1")
nrow(BTomEgg)
write.table(BTomEgg[,2],file = "List_Ortho_B_LA_MM_Solyc.txt",row.names = FALSE, col.names = "Genes_Name")

BShared = merge(BPepEgg, listSoly.Tom.B, by = "V1")
nrow(BShared)
write.table(BShared[,2],file = "List_Ortho_B_LA_MM_PM_Solyc.txt",row.names = FALSE, col.names = "Genes_Name")

TOT = c(as.character(BShared[,2]),as.character(BTomEgg[,2]), as.character(BTomPep[,2]), as.character(BPepEgg[,2]))
length(TOT)
write.table(TOT,file = "List_Ortho_B_TOT_Solyc.txt",row.names = FALSE, col.names = "Genes_Name")


nrow(Tom.A);nrow(Tom.B);nrow(Egg.A);nrow(Egg.B);nrow(Pep.A);nrow(Pep.B)

```

