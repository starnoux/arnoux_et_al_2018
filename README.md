# Pipeline guideline {arnoux_et_al_2018}  

### In this github, you will find all the scripts that were used for the data analyses in our paper.  

The github is devided in two sections, the one treating the raw data to call the variants and the other treating all the statistical abalyses with R.

### A. From a RNA-seq to a clean vcf: Sequences alignement, variant calling and paralogous site cleaning.  
**This first guideline is giving details on the steps to align the sequences and obtain clean and reliable variant call format files (vcf).**    

*WARNING:* Use the guidline **_A_RNAseq_to_vcf.md_** that is the detailed pipeline script.
  
All the scripts are named from *01_fastqc_quality.py* to *10_subVCFfromFilter.vbeta.py* and they are split in between 3 folders, representing the different alignment and cleaning steps:  
1.  Cleaning_raw_data
2.  Mapping_And_SNPs_Calling
3.  Paralog_Cleaning


### B. Statistical analyses on vcf.  

This second guideline is giving details on the statistical analyses:  
*WARNING:* Use the guidline **_B_Statistical_Analyses.md_** that is the detailed pipeline script.
  
4.  Statistical analyses  




 
 




