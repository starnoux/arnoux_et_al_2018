#!/usr/bin/Rscript

#title : 'Skin_Coverage_vbeta'
#author: "S. Arnoux"
#date : "January 2018"


library("optparse")
rm(list=ls())
option_list = list(
make_option(c("-d", "--directory"), type="character", default=NULL, 
              help="directory location", metavar="character"),
make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input file type epluche.txt", metavar="character"),
make_option(c("-o", "--out"), type="character", default="out.list", 
              help="output file name [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$directory)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

InputEpluche  = opt$input
OutputList = opt$out

#"liste_des_contig_tomato_couverts_a_6.liste"


# This script is here to filter the minimum coverage


setwd(opt$directory) 

don <-read.table(InputEpluche, header=TRUE, sep="\t", dec=".")

dim(don)
# 27170 ## 31847 

names(don)
#  1 "contig"   
#  2 "nb_snp" : SNP number per contig  
#  3 "nb_htz"  : number of heterozigous situations observed on SNP number (nb_snp)
#  4 "nbobs"  : number of total relevant situations observed  on SNP number (nb_snp)
#  5 "nb_all_ref"  : number of occurance when the reference allele appeared  
#  6 "rapport"  : heterozygosity  ## NOT PRESENT
#  7 "InbreedCoef" mean of InbredCoef on contig
tiff(paste('HistogramInfo_', InputEpluche, '.tiff', sep = ""), width = 2700, height = 1800, units = "px", res = 400)
par(mfrow = c(3,1))
hist(don[,6] , main = "InbreedCoef", xlab = "Coefficient", col = "gray")
hist(log10(don[,2]) , main = "Log 10 Nbr de SNP", xlab = "Nbr de SNP", col = "gray")
dev.off()

# Here we add an element to the epluche file with the Individual mean of 
don$Ind_moy<-don[,4]/don[,2]
hist(don$Ind_moy, main = "Ind_moy", xlab = "")


# We keep the given threshold , 6 (because there is a lot of wild accessions)

subdon<-don[which(don[,7]>6),]

dim(subdon)
# 13186 (47%) ## 27221 (85%) # 29592

write.table(subdon, OutputList,row.names=T, quote=F)
