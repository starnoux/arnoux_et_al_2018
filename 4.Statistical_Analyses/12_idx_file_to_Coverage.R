---
title : 'Preliminary data formating'
author: "S. Arnoux"
date : "January 2018"
output : html_document
---

```r
rm(list=ls())

## Give the working folder
PATH_STAT = "/path/to/idx_files/"  

setwd(PATH_STAT)  

## This first file will contain all the third columns of your idx files (meaning the coverage).
Name_File = "/path/to/Deseq_Folder/Exp_Coverage.txt" 

## This second file will contain all the fourth columns of your idx files (meaning the coverage).
Name_File_UnMapped = "/path/to/Deseq_Folder/Unmapped_Coverage.txt"

########## Tomato ##############
# LA8001=idx_LACMVSel.txt
# LA8002=idx_LAHirB.txt
# LA8003=idx_LAPI247087.txt

## Get the list of files starting by 'idx_' in the folder specified 
idx_List <- list.files(path = PATH_STAT, pattern="^idx_",  full.names=TRUE)

## Read with read.table the list and create a grouped list
my_file_list <- lapply(idx_List, read.table)

## Name the list with the accession names and rempves the 'idx_' and other '_' or '.txt' from the name to get a clean name. example 'LA1028'
Name_list <- list.files(path = PATH_STAT, pattern="^idx_",  full.names=FALSE)
names(my_file_list) <- gsub("_","",(gsub(".txt","",(gsub("idx_","",Name_list)))))

## Extraction of the Gene Names
Start_df <- as.data.frame(my_file_list[[1]]$V1)
names(Start_df) <- "Genes" 

## Function to get the mapped reads
TOT <- NULL
IdX_2_depth_File <- function(Acc){
  Col = Acc$V3
  names(Col) = names(Acc) 
  TOT = c(TOT,Col)
}
Data_Cov = as.data.frame(lapply(my_file_list, IdX_2_depth_File))
head(Data_Cov)

write.table(cbind(Start_df,Data_Cov),file=Name_File, row.names = FALSE)

## Function to get the unmapped reads
TOT <- NULL
IdX_2_Uncov_File <- function(Acc){
  Col = Acc$V4
  names(Col) = names(Acc) 
  TOT = c(TOT,Col)
}
Data_UnCov = as.data.frame(lapply(my_file_list, IdX_2_Uncov_File))
head(Data_UnCov)

write.table(cbind(Start_df,Data_UnCov),file=Name_File_UnMapped, row.names = FALSE)

## Accesions file
write.table(cbind(Name=names(my_file_list),POP=c(rep("x",length(names(my_file_list))))),file = "Acc_Loc.txt", row.names = FALSE)

```
