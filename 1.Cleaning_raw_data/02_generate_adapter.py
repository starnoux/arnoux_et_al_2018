#!/usr/bin/python2.7

import re
import glob
import os
import sys
import argparse

################
#   Adapters   #
################

truseq_index_1 = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" 											#part 1 of the indexed adaptor
truseq_index_2 = "ATCTCGTATGCCGTCTTCTGCTTG"														#part 2 of the indexed adaptor

####################
#   Informations   #
####################
   
courriel = "francois-xavier.babin@outlook.fr"

#creation of a parser which gather the arguments of the program 
parser=argparse.ArgumentParser(																						
	description=""" program to generate the adapter file""",
	epilog="""
	contact: """+courriel+"""\n\n""")

parser.add_argument("-i", 				type=str, 		default="/path/to/input/folder/", 					help="input directory (default: %(default)s)")				#input_directory settings
parser.add_argument("-o", 				type=str, 		default="/path/to/input/folder/adapter/", 								help="output directory (default: %(default)s)")				#output_directory settings	 														help="min size to keep (default: %(default)s)!")			#minimum size of reads to keep settings
parser.add_argument("-pr", 				type=str, 		default="PM", 										help="prefix (default: %(default)s)!")						#input file settings (the default value is an used file in the project)
parser.add_argument("-contaminants", 	type=str, 		default="/path/to/contaminants_illumina_enriched.fa", 		help="contaminants file (default: %(default)s)")

args=parser.parse_args()

#################
#   Functions   #
#################

#dictionnary of complementary bases
compl={	
	"A":"T",
	"C":"G",
	"G":"C",
	"T":"A",
	"a":"t",
	"t":"a",
	"g":"c",
	"c":"g"
	}

#function to get the reverse complement of a sequence
def reverse_complement(seq) :

	cseq=""
	i= len(seq[1:])
	while i >= 0:
		cseq += compl[seq[i]] 
		i -= 1
	return cseq

#Function try except
def try_except(message,command):
	try:
		print(message+"\n")
		command = os.system(command)

		if command == 0:
			print(message+" done ...\n")
		else:
			print(message+" failed ...\n")
	except OSError as e:
		print("Execution failed ", e)

#####################
#   Verifications   #
#####################

#checking if the contaminant file exists
if not os.path.isfile(args.o+args.contaminants) :
	print('the contaminant file ('+args.contaminants+') doest not exist !')

###################
#   Instructions  #
###################

paired_fastq_files = glob.glob(args.i+args.pr+"*1.fastq")
index_adapter=""

for files in paired_fastq_files :
	
	file_name = os.path.splitext(os.path.basename(files))[0]
	index = re.search("[A-Z0-9]{6}", file_name).group()
	index_adapter += ">"+file_name+"_forward\n"+truseq_index_1+index+truseq_index_2+"\n"+">"+file_name+"_reverse\n"+reverse_complement(truseq_index_1+index+truseq_index_2)+"\n"

	#creation of the adapter file
	try_except("Creating adapter file ...","cat "+args.o+args.contaminants+" > "+args.o+"adapters.txt")

adapter_file=open(args.o+"adapters.txt", "a")
adapter_file.write(index_adapter)

