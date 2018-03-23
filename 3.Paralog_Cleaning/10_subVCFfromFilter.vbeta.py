#!/usr/bin/python2.7

import os
import glob
import sys
import argparse
import re

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


courriel = "ste.arnoux@gmail.com"
 
parser=argparse.ArgumentParser(                                                           #creation of a parser which gather the arguments of the program
    description=""" This program will generate sequences for all the fasta files in sequences with a name like this LD_LL0000.
    contact: """+courriel+"""\n\n""")

#parser.add_argument("-l",type=str,default="P",                                              help="prefix letter of the LD or LW or LO (default: %(default)s)")
parser.add_argument("-i",           type=str,       default="/path/to/Folder/input/",                         help="input directory (default: %(default)s)!")
parser.add_argument("-vcf",           type=str,       default="NewVcf_Filtered.vcf",                         help="input directory (default: %(default)s)!")
parser.add_argument("-list",           type=str,       default="NewVcf_Filtered.out.Merged.res",                         help="list of selected elements (default: %(default)s)!")
parser.add_argument("-new",           type=str,       default="Filtered_liste.tab.txt",                         help="New vcf file name (default: %(default)s)!")
parser.add_argument("-thr",           type=int,       default=0.05,                         help="theshold (default: %(default)s)!")

args=parser.parse_args()


with open (args.i+args.list, "r+b" ) as source:
	next(source)
	with open(args.i+args.new, 'w') as targetG , open(args.i+'Para_liste.txt', 'w') as targetPara:
		for line in source:
			col=line.split(" ")
			if col[0] == 'NA':
				pass
			if col[2] == 'NA':
				pass
			elif (int(float(col[5])) > args.thr):
				targetG.write(col[0].strip('"')+ '\t' +col[1]+ '\n' )
			elif (int(float(col[5])) <= args.thr):
				#targetPara.write(col[0]+ ' ' +col[1]+ ' \n' )
				targetPara.write(line)
				#print line #(Only if it is necessary to check if the script works fine)


source.close()
targetG.close()
targetPara.close()

### FOLLOWING DO 
# bgzip LA_BCalib_rmindels_mP10_mQ20_nost_CSFILTER_Py_6.vcf
# bcftools index LA_BCalib_rmindels_mP10_mQ20_nost_CSFILTER_Py_6.vcf.gz
# bcftools view -R Filtered_liste.tab.txt LA_BCalib_rmindels_mP10_mQ20_nost_CSFILTER_Py_6.vcf.gz -o LA_BCalib_BOBOTEST.vcf










							

