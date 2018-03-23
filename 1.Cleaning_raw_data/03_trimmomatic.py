#!/usr/bin/python2.7

import glob
import os
import sys
import argparse
import re

################
#   Adapters   #
################

truseq_index_1 = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" 											#part 1 of the indexed adaptor
truseq_index_2 = "ATCTCGTATGCCGTCTTCTGCTTG"														#part 2 of the indexed adaptor

####################
#   Informations   #
####################
   
courriel = 'francois-xavier.babin@outlook.fr'

#creation of a parser which gather the arguments of the program 
parser=argparse.ArgumentParser(																						
	description=''' program to trim adapters and clean fastq files in directory''',
	epilog="""The input files must be in 'fastq' format.
	contact: """+courriel+"""\n\n""")

parser.add_argument('-i', 			type=str, 		default="/path/to/input/folder/", 										help='input directory (default: %(default)s)')				#input_directory settings
parser.add_argument('-input_type', 	type=str, 		default='phred33', 													help='phred input type (phred33/64) (default: %(default)s)!')
parser.add_argument('-f', 			type=str, 		default="/path/to/input/folder/adapter/adapters.txt",									help='adapters file (default: %(default)s)!')
parser.add_argument('-q', 			type=int, 		default=20, 														help='quality level (default: %(default)s)!')				#quality level settings
parser.add_argument('-v', 			type=int, 		default=7, 															help='overlap (default: %(default)s)!')						#overlap settings
parser.add_argument('-m', 			type=int, 		default=20, 														help='min size to keep (default: %(default)s)!')			#minimum size of reads to keep settings
parser.add_argument('-pr', 			type=str, 		default="LA", 														help='prefix (default: %(default)s)!')						#input file settings (the default value is an used file in the project)
parser.add_argument('-t', 			type=int, 		default=8, 															help='threads (default: %(default)s)!')
parser.add_argument('-trimmomatic', type=str, 		default="/path/to/trimmomatic-0.33.jar", 			help='trimmomatic location (default: %(default)s)!')

args=parser.parse_args()

###################
#   Instructions  #
###################

paired_fastq_files = glob.glob(args.i+"/"+args.pr+"*1.fastq")

for files in paired_fastq_files :
	
	file_name = os.path.splitext(os.path.basename(files))[0]

	input_forward = args.i+file_name[:-1]+"1.fastq"
	input_reverse = args.i+file_name[:-1]+"2.fastq"
	paired_output_forward = args.i+"paired_"+file_name[:-1]+"1.fastq.gz"
	paired_output_reverse = args.i+"paired_"+file_name[:-1]+"2.fastq.gz"
	unpaired_output_forward = args.i+"unpaired_"+file_name[:-1]+"1.fastq.gz"
	unpaired_output_reverse = args.i+"unpaired_"+file_name[:-1]+"2.fastq.gz"

	#run trimmomatic
	try : 
		print('filtering  ...')  											
		trimmomatic_command = os.system("java_1.8 -jar "+args.trimmomatic+" PE -threads "+str(args.t)+" -"+args.input_type+" "+input_forward+" "+input_reverse+" "+paired_output_forward+" "+unpaired_output_forward+" "+paired_output_reverse+" "+unpaired_output_reverse+" ILLUMINACLIP:"+args.f+":2:30:10 MINLEN:"+str(args.m)+" SLIDINGWINDOW:"+str(args.v)+":"+str(args.q))		
		if trimmomatic_command == 0 :
			print('done')
		else :
			print('trimmomatic command failed !')
	except OSError as e:
		print("Execution failed ", e)
