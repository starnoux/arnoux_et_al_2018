# A. From a RNA-seq to a clean vcf:
## Sequences alignement, variant calling and paralogous site cleaning.

#### Example of input file name:  
  
The raw data were named as following example:  
> 'PM0910_AGTCAA_L007_R1.fastq.gz'  
> 'PM0910_AGTCAA_L007_R2.fastq.gz'  
- The species prefix is 'PM' and the 4 following number '0910' are the individual ID.  
- The 'AGTCAA' is the adaptator sequence that will be removed in the sequences.  
- The 'L007' is an indication on the line where it was sequenced, this is not used in our pipeline.  
- The '_R1.' and '_R2.' are the informations about the sequences respectively Forward and reverse sequencing.  

### 1. Cleaning the raw data. 
All necessary files or scripts are in the folder **_1.Cleaning_raw_data/_**
  
  **a. Unzip**  
Your input is a fastq.gz and therefore before everything we'll unzip these files.  
You need zcat (gzip) installed on your system.  
```bash
./Unzip_fastq_files.py -i /path/to/input/folder/ -pr PM
```   
> ###### *python v2.7; zcat v1.6*  

  **b. Quality control**  
You want to check the quality of your data and therefore you will produce fastqc data files.  
For this step you can find the file 'contaminants_illumina_enriched.fa' but I would advise to get this details depending on your sequencing methods and therefore on their specific contaminants.
```bash
./01_fastqc_quality.py --help
```
> ###### *FastQC v0.11.5*  

  **c. Generation of adapter**    
To generate the specific adapter files you need to know the sequencing indexes of your data set. And to report it inside the python script.
```
./02_generate_adapter.py --help
```
  
  **d. Trimming of possible residual adapter**  
Once the adapters are created, your files will be trimmed and renamed 'paired_*File-Name*'.
```
./03_trimmomatic.py --help
```
> ###### *Trimmomatic v0.33*  

  **e. Unzip the 'paired_File'**  
Your input is a paired_*File-Name*.fastq.gz and therefore before to map it you will unzip it.
```
./Unzip_fastq_files.py -i /path/to/input/folder/ -pr paired_
```   

### 2. Mapping data and SNPs Calling  
All necessary files or scripts are in the folder **_2.Mapping_And_SNPs_Calling/_**
  
- **Python pipeline**  
This python script had been written to obtain the SNP calling following the GATK best practices. Here are the hard filtering parameters: '*QD<2.0||FS>60.0||MQ<40.0||ReadPosRankSum<-8.0*'. It might eventually have to be regulated according to your folders and system. But don't give up and if all goes wrong just take it step by step, it will work out I promise.
```  
./04_mapping_refCDS.py --help 
```   
> ###### *picard v2.0.1; GenomeAnalysisTK v3.8; bwa v0.7.12; samtools v0.1.19; vcftools v0.1.12b*

### 3. VCF Filtering and Paralogous sites Cleaning    

All necessary files or scripts are in the folder **_3.Paralog_Cleaning/_**
  
  **a. Preliminary steps**  
  + Indels Removals
Now  that you obtained you vcf, you will need to filter the indels and the SNPs that don't behave in a binary way.
We will thereafter use vcftools to filter all indels '-remove-indels', to remove the quality '-minQ' and depth 'min-meanDP'. For more details on possible options take look at the manual.
*[Vcftools manual](http://vcftools.sourceforge.net/man_latest.html)*
```bash
./vcftools --vcf file.vcf --out filtered_file  --remove-indels --remove-filtered-all --min-meanDP 10 --minQ 20 --remove-filtered-geno-all --recode  
``` 
   + Stars '\*' Removals
Don't ask me why but it seems that some the deletion are represented by stars, and these stars are no welcome in the next step, so if vcftools didn't get rid of them use the next step.
```bash
sed "/\*/d" filtered_file.recode.vcf > filtered_file.recode.nost.vcf
```   
  **b. vcf file dissection**  
The given threshold value will be excluding all heterozygous sites that have less than the ‘threshold value’ reads over all the individuals. In our case if 'threshold = 4' then the site must be present in at least 4 alleles of the total number of the population allele number. (Diploid -> total allele over population = individuals \* 2)  
* Script from Sylvain Glémin (sylvain.glemin@univ-montp2.fr) / C++ version included in the read2snp software 
[READS2SNP](http://kimura.univ-montp2.fr/PopPhyl/resources/tools/)*.
```bash
perl 05_willy_waller_2006.vbeta.pl /path/to/filtered_file.recode.nost.vcf dissect_file.output.txt 4_or_otherThreshold
```  
  **c. Skin the bad coverage**  
In this step the coverage thershold 6 will be used to remove to low covered contigs, and a list will be created.      
```bash
./06_Skin_Coverage_vbeta.r -d /path/to/Folder/input/ -i dissect_file.output.txt -o out.list
```
  **d. Creation of a filtered vcf**  
In this step the filtered contigs *'out.list'* from the previous step will be used to filter the vcf file.  
```bash
./07_subVCFfromList.vbeta.py -i /path/to/Folder/input/ -vcf filtered_file.recode.nost.vcf -list out.list -new NewVcf_Filtered.vcf
```
  **e. formatting of the vcf.out file**   
In this step there is mostly formatting the previous file to make the next step possible.
* Script from Sylvain Glémin (sylvain.glemin@univ-montp2.fr) / C++ version included in the read2snp software [READS2SNP](http://kimura.univ-montp2.fr/PopPhyl/resources/tools/)*.
```bash
perl 08_formatting_parafilter.vbeta.pl /path/to/NewVcf_Filtered.vcf /path/to/NewVcf_Filtered.vcf.out  
```
  **f. Estimating the paralagous sites and creating the to-filter.list**   
In this script, we create the files to manage all the filtering in parallele subfiles. This step was used on the server. We used the following paramters for the resampling and the paralogous testing.   
>err<-0.001 # error threshold <br/>fis<-0.5 # FIS threshold <br/>PREC<-0.001 # optimisation threshold <br/>POLYPARA<-TRUE # Booleen to know if we want the polymorphism details of the paralogous sites <br/>Threshold<-0.01 # p.value threshold for the difference between the model with or without paralogs above which we do not test the paralog polymorphisms  

* Script from Sylvain Glémin (sylvain.glemin@univ-montp2.fr) / C++ version included in the read2snp software [READS2SNP](http://kimura.univ-montp2.fr/PopPhyl/resources/tools/)*.
```bash
./09_SplitFilterParalog_vbeta.R -d /path/to/Folder/input/ -i NewVcf_Filtered.vcf.out
```
###### WARNINGS: I had to do manually the merging of the .res files, that happened most probably because of the resources allocated to the server. Just check that the number of lines in your .vcf.out.Merged.out is corresponding to the sum of the .vcf.out.  As this could happen depending on the server you are working on, just be careful and do the mergin steps on R manually and not with the automatized script.    
  **g. Paralogous countig extraction**  
The pralalogous countig have been statistically tested with the previous steps and the list of paralagous can be now extracted according to the threshold we define (Here we use 0.05 as default but it can be change with the option -thr).  
```bash  
./10_subVCFfromFilter.vbeta.py -i /path/to/Folder/input/ -vcf NewVcf_Filtered.vcf -list NewVcf_Filtered.vcf.out.Merged.res -new Filtered_liste.tab.txt  
```  
 **h. Last filtering**  
To have a rapid and efficient filtering of our big vcf files from the paralogous countig sites, we will use bcftools. *[bcftools manual](https://samtools.github.io/bcftools/bcftools.html)*  
```bash  
bgzip NewVcf_Filtered.vcf
bcftools index NewVcf_Filtered.vcf.gz
bcftools view -R Filtered_liste.tab.txt NewVcf_Filtered.vcf.gz -o VCF_ParalogFreeAndFiltered.vcf
```
> ###### *bcftools v1.3.1*  

