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
    description=""" This program will clean the paralogous from the vcf inital. """+courriel+"""\n\n""")

#parser.add_argument("-l",type=str,default="P",                                              help="prefix letter of the LD or LW or LO (default: %(default)s)")
parser.add_argument("-i",           type=str,       default="/path/to/input/folder/",                         help="input directory (default: %(default)s)!")
parser.add_argument("-vcf",           type=str,       default="filtered_file.recode.nost.vcf",                         help="input directory (default: %(default)s)!")
parser.add_argument("-list",           type=str,       default="out.list",                         help="list of selected elements (default: %(default)s)!")
parser.add_argument("-new",           type=str,       default="NewVcf_Filtered.vcf",                         help="New vcf file name (default: %(default)s)!")
parser.add_argument("-thr",           type=str,       default='6',                         help="threshold used in the list of contig (default: %(default)s)!")

args=parser.parse_args()


with open (args.i+args.list, "r+b" ) as source:
	with open(args.i+'liste'+args.thr+'.txt', 'w') as target:
		for line in source:
			col=line.split(" ")
			target.write(col[0]+ ' \n')


source.close()
target.close()
			
with open(args.i+'liste'+args.thr+'.txt', "r+b" ) as liste, open(args.i+args.vcf, 'r+b') as vcf, open(args.i+args.new, 'w') as vcfnew:
	next(liste)
	toto = list(word.split(' \n')[0] for word in liste)
	#print (list(word.split(' \n')[0] for word in liste))
	for line in vcf:
		if line[0]=='#':
			vcfnew.write(line)
		else :
			for t in toto:
				if t in line:
					vcfnew.write(line)
