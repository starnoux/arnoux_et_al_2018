#!/usr/bin/python2.7

import os
import argparse
import glob
from os.path import basename

####################
#   Informations   #
####################

courriel = 'ste.arnoux@gmail.com'
 
parser=argparse.ArgumentParser(																	#creation of a parser which gather the arguments of the program
	description=''' program to unzip forward and reverse compressed files ''',
	epilog='''The compressed files must be in 'fastq.gz' format.\n\ncontact: '''+courriel+'''\n\n''')

parser.add_argument('-i', 	type=str, 		default="/path/to/input/folder/", 	help='input_directory (default: %(default)s)')
parser.add_argument('-pr', 	type=str, 		default="PM", 				help='prefix (default: %(default)s)')						

#parsed arguments are gathered in args
args=parser.parse_args()																		

#####################
#   Verifications   #
#####################

#if the input directory doesn't exists it will create it
if not os.path.exists(args.i):
	try :
		print('the '+args.i+' doesn\'t exists ... creating an empty directory')
		input_folder_command = os.system('mkdir '+args.i)
		
		if input_folder_command == 0 :
			print('done')
		else :
			print(args.i+' creation failed !')
	except OSError as e:
		print("Execution failed ", e)
	
###################
#   Instructions  #
###################

fastqgz_files = glob.glob(args.i+args.pr+"*.fastq.gz")

for files in fastqgz_files:																
		
	file_name = basename(files).split('.')[0]									#remove the path from the file

	try : 
		print("Unzipping "+files+" ...")
		unzip_command = os.system("zcat "+files+" > "+args.i+"/"+file_name+".fastq")	
		
		if unzip_command == 0 :
			print('done')
		else :
			print("Unzipping "+files+" failed !")
			break
	except OSError as e:
		print("Execution failed ", e)
