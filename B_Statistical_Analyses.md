# B. From vcf to statistical analyses
## Producing figures and statistical analyses for each of the species.
  
### 4.  Statistical_Annalyses.  
  
All necessary files or scripts are in the folder **_4.Statistical_Annalyses/_**
Keep in mind if you use these scripts that you might need to check the figure axis limits to make sure that you are in the good range.
  
### a.  Reference annotation  
  
#### -> Translate the reference  
  
Your input is a reference genome fasta.  
```bash
export PATH=/usr/bin/java_1.8:$PATH]  
/path/to/interproscan-5.26-65.0/bin/nucleotide/translate -i /path/to/reference.fa -o /path/to/reference.prot.fa
```
> ###### *interproscan v5.26-65.0; java 1.8*  
  
#### -> Chunk the file in smaller subfiles  
  
In this step we used a script developped by Don Gilbert and named **_split_multifasta.pl _** that is accessible in the [Indiana University webpage](http://iubio.bio.indiana.edu/gmod/genogrid/scripts/split_multifasta.pl). And detailed in the [project webpage](http://iubio.bio.indiana.edu:7122/gmod/genogrid/). This script split the file in chunk of 30000 sequences as this is the limit for interproscan.  
```bash
./split_multifasta.pl --in /path/to/reference.prot.fa --output_dir=/path/to/reference_folder/ --seqs_per_file=30000
```
  
#### -> Bash loop to run the inteproscan on the subfiles  
  
WARNING: THE FOLLOWING IS DANGEROUS BECAUSE RUN EVERYTHING IN PARALLELE ... If you have 10 files then it is ok otherwise be CAREFUL !! if you have 60... That would be 60 jobs * cpu 8
This loop has to be adjusted according to the server you use and its memory capacities.
I renamed the 'subfiles.fa' -> 'subfiles.fsa' to be able to only select them with the extension name.
```bash
for i in {1..9}; do
  echo $i
  for j in `ls -rt SME_r2.5.1_Test/${i}*.fsa`; do
    echo ${j}
    num=`echo $(echo $j | sed 's/SME_r2.5.1_Test\///' | sed 's/.fsa//')`
    nohup /path/to/interproscan-5.26-65.0/interproscan.sh -t p -cpu 8 -i ${j} -dp -iprlookup -goterms --pathways -b /path/to/Annotation/ref_interproscan_${num} & > nohup${j}.out
  done
done
```
  
#### -> Merging the interpsocan '.tsv' result files  
  
This file will help merge the files and then sort the genes per name as in the reference genome.
```bash
cat ref_interproscan_*.tsv > ref_interproscan_Glob.tsv
sort -nk1 ref_interproscan_Glob.tsv  > ref_interproscan_Ord_Glob.tsv
```  
#### -> Extract only pfam annotations    
  
The files that will be produced will be used in the following statistical analyses.
  
The following is managed for the SMEL reference genome, the TAB needs to be replaced with a real 'tab':
HELP FOR MAC USERS ctrl+v+tab writes down a tab ;)
```bash
grep 'Pfam' ref_interproscan_Ord_Glob.tsv | grep 'GO\:' | awk -F'\t' '{print $1, $14}' | uniq | sed 's/\|GO:/\,/g' | sed 's/GO\:/TAB/g'| awk -F"_| " '{print $1,"_",$2,$4}' | sed 's/ //g' > ref_interproSc_GO_pfam.txt
``` 
The following is managed for the CaZl1 reference genome:
```bash
grep 'Pfam' ref_interproscan_Ord_Glob.tsv | grep 'GO\:' | awk -F'\t' '{print $1, $14}' | uniq | sed 's/\|GO:/\,/g' | sed 's/GO\://g'| awk -F"_| " '{print $1,$3}' |  sed 's/ /	/g'> ref_interproSc_GO_pfam.txt
```  
  
#### -> Obtain the Gene Locations from the reference '.gff' file 
  
It will produce a ref.Gene_Loc.tab that will be used for the future graphs.
```bash
./11_get_Loc_Chrom_tab_gff3.py -i /path/to/reference/ -tabin reference.gff3 -tabout ref.Gene_Loc.tab
```  

### b. Deseq analyses
  
#### -> Preliminary data formating  
  
At first, we want to use the '.idx' files produces during the alignment to extract and concatenate the informations on expression. 
```bash
cd /path/to/.idx_Files/
```
Then run **MANUALLY** (read it and follow line per line) the script **_12_idx_file_to_Coverage.R_**  
> ###### *R v3.3.2*   
It will help you understand what the script does, and you can check that your files have the same format as expected.  
At the end of this step modify the created file 'Acc_Loc.txt' to define the family, group or origin location of your accessions.   

#### -> Deseq Analyses  
  
These analyses will give you information on the genes that are up- and down-regulated. You can run the scripts entirely as far as you fill up the preliminary informations. Just open and check out the instructions of the script **_13_Deseq_Pair_Pop.R_**.    
  
This script will statistically estimate the level of gene expression differences and define which genes are Differentially Expressed Genes (DEG).  
 -  It gives lists of down- and up-regulated genes *'PopA_PopB.UP.txt'* where the PopB is up-regulated compare to the PopA reference.  
 -  It gives list of Gene Ontology of these previous DEG.  
 -  And it produces figures:  
    -  Principal Component Analysis (PCA) of gene expression level per individuals
    -  Heatmap of expression level per gene of DEG 
    -  Standiard deviation mean plot
    -  Plot of data distribution between the raw data and the rlog transformed data
    - ... *(cf.details in script)*       
  
#### -> GLM regression models     
In the following script on R, **_14_GLM_regression_models.R_**, we used the DNAsp results and tested the effect of Tajima's D on the expression and the effect of Pi on the expression (we tested in both cases the reciprocal hypothesis).   
  
#### -> Genome-Wide figures  
  
In the following script on R, **_15_DNAsp_Stat_Fig.R_** , we used the statistics per gene that were produced with DNAsp per population of each species. If you use this script, you have to know that we 'smoothed' the data over 50 genes with the 'rollmean' function and therefore that you might change this function parametes.    
As well, it is pretty important to check the graph/figure limits in order to have the plot in the range of the plot.  
  
**WARNING:** One element HAS to be changed according to the reference annotation. line 59 & 66, you'll need to indicate the chromosome location. Using the following 'substr' function you need to give a precise indication as in the following example. 

>ex.1: chr <- substr(dtWild$Gene,7,8)  #dtWild is the dataframe of the wild (PopB) DNAsp results  
>SMEL_0*03*g194330.1.01  
>chr <- substr(dtWild$Gene,6,7)   
>Solyc*12*g100353.1.1  
  
As well, this script will give you all the mean per chromosomes and the data informations (nucleotid diversity(*Pi*) and Tajima's D). Once you know the chromosomes that have statistics relevant changes (cf. *14_GLM_regression_models.R* ) that you want to put in front, then just write it up in the 'SIGNIFICANCE LINES' section.
  
#### -> Gene ontology enrichement analysis on Pi 'shifted' genes  

The Pi shifted genes, were defined as genes that would have a high nucleotid diversity in one pop and an extremly (or null) nucleotid diversity in the other population. The shifted group are split in rwo. The group A that are genes highly diverse in the crop population but not in wild population, and the group B that are genes with no nucleotid diversity in crop population but really high nucleotid diversity in wild population.   

In the following script on R, **_16_Pi_GO_Enrichment.R_** , we will statistically test the gene ontology enrichment of our shifted gene groups A and B. A plot will be produced with the abline following the threshold and the coloured dot being the shifted genes.

#### -> Production of PCA plot (on Genetic diversity)   
  
Using the 'SNPRelate' library on R, we produced PCA plots, the script is **_17_PCA.R_**.

 




