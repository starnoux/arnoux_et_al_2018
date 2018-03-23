#!/usr/bin/python2.7

import os
import sys
import argparse
import re

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

### TRY TO GET THE FORMAT TAB WITH THE GENE LOCATIONS ##
courriel = "ste.arnoux@gmail.com"
 
parser=argparse.ArgumentParser(                                                           #creation of a parser which gather the arguments of the program
    description=""" This program will generate a tab-delimited file with the genes names and their physical location on the chromosomes. 
        """+courriel+"""\n\n""")

#parser.add_argument("-l",type=str,default="P",                                              help="prefix letter of the LD or LW or LO (default: %(default)s)")
parser.add_argument("-i",           type=str,       default="/path/to/reference/",                         help="input directory (default: %(default)s)!")
parser.add_argument("-tabin",           type=str,       default="reference.gff3",                         help="this is the input name (GFF FILE) (default: %(default)s)!")
parser.add_argument("-tabout",           type=str,       default="ref.Gene_loc.tab",                         help="output name (if not yet created it will be) (default: %(default)s)!")

args=parser.parse_args()

with open (args.i+args.tabin, "r+b" ) as source, open(args.i+args.tabout, 'w+') as target:
	ID       = ""
	# This is the ID that needs to be CDS.
	GName    = ""
	lstcol   = ""
	line_out = ""
	# Before to start the loop the best way is to define the variable as previously.
	for line in source:
		if line[0] != "#" :  
			lstcol = line.split("\t")
			if len(lstcol) > 1 :
			#Out of range because of the comment lines that are not spliteed in columns!!! I could have test at first if the line started with # see previous comment.
				ID = lstcol[2]
				if ID == "mRNA" : 
					if len(lstcol[8].split(";")) > 1 :
						GName = lstcol[8].split(";")[0].split("=")[1]
					line_out = GName + '\t' + lstcol[3] 
				# Here we defined what is the firt alternate allele in the 4th column
		#AA = lstcol[4][:-1]
		# the last columns of interest is the AA. In order to get we need to observe what the column 4 indicate but PAY ATTENTION at the end of the file there are some \n that return the line to the next one SO we indicate [:-1] in order to remove the last character. 
		# We want to have as return, all the lines in the source (for line in source:), AND on each line we want to have the column 0 and 1, seperated per tabs. 
					target.write(line_out + '\n')
		# In this last line after demanding to get the ref,alt1 or alt2 added at the end of the line ( line += AA ) we add as well a \n that will create a seperation for each line with each others. 


source.close()
target.close()
# time to close the files, it is done baby!
