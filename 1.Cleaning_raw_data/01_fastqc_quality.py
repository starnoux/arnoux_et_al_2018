#!/usr/bin/python2.7

import os
import argparse
import glob

####################
#   Informations   #
####################
   
courriel = 'francois-xavier.babin@outlook.fr'

#creation of a parser which gather the arguments of the program 
parser=argparse.ArgumentParser(																						
	description=''' program to generate fastqc qualities corresponding to fastq files in the specified directory''',
	epilog="""The input files must be in 'fastq' format.
	contact: """+courriel+"""\n\n""")

parser.add_argument('-i', 	type=str, 		default="/path/to/input/folder/", 						help='input directory (default: %(default)s)')				
parser.add_argument('-o', 	type=str, 		default="/path/to/input/folder/qualities/", 				help='output directory (default: %(default)s)')			
parser.add_argument('-t', 	type=int, 		default=1, 												help='threads (beware do not use too many threads if your computer is not able to !!!) (default: %(default)s)')					
parser.add_argument('-pr', 	type=str, 		default="MM", 											help='prefix (default: %(default)s)')
parser.add_argument('-f', 	type=str, 		default="/path/to/fastqc", 										help='fastqc location (default: %(default)s)')
parser.add_argument('-c', 	type=str, 		default="/path/to/contaminants_illumina_enriched.fa", 	help='contaminants file (default: %(default)s)')

args=parser.parse_args()

################
#   Function   #
################

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

#if the input directory doesn't exists it will create it
if not os.path.exists(args.i):
	try_except("the "+args.i+" doesn\'t exists ... creating an empty directory", "mkdir -p "+args.i)

#if the qualities directory doesn't exists it will create it
if not os.path.exists(args.o):
	try_except("the "+args.o+" doesn\'t exists ... creating an empty directory","mkdir -p "+args.o)

if not os.path.isfile(args.c) :
	print('the contaminant file ('+args.c+') doest not exist !')
		
###################
#   Instructions  #
###################

#fastqc qualities generation
try_except("Generating fastqc qualities ...","parallel --gnu -j"+str(args.t)+" '"+args.f+" -o "+args.o+" --nogroup -c "+args.c+"' ::: "+args.i+args.pr+"*.fastq")
