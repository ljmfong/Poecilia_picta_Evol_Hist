#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
#==============================================================================
import argparse
import sys
import csv
import os
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infolder", type=str,
					help="Infolder of Trinity files")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def list_folder(infolder):
	'''Returns a list of all files in a folder with the full path'''
	return [os.path.join(infolder, f) for f in os.listdir(infolder) if f.endswith(".fasta")]

def read_transcripts(infile):
	SG={}
	try:
		with open(infile, "r") as trinity:
			for line in trinity.readlines():
				if line[0] == ">":
					name = line[1:].rstrip().split()[0]
					header = line[1:].rstrip()
					SG[name] = ["", header]
				else:
					SG[name][0] += line.rstrip()
			return SG
	except IOError:
		print "!----ERROR----!"
		print "File %s does not exit!" % infile
		sys.exit(1)
	except KeyboardInterrupt:
		sys.exit(1)

def get_isoform_gene(isoform_dict):
	dict = defaultdict(list)
	for isoform in isoform_dict:
		gene = isoform.split("_")[0]+"_"+isoform.split("_")[1]+"_"+isoform.split("_")[2]+"_"+isoform.split("_")[3]
		dict[gene].append(isoform)
	return dict

def average(gene_dict):
	isoform_number=[]
	for gene in gene_dict:
		isoform_number.append(len(gene_dict[gene]))
	average = float(sum(isoform_number))/len(gene_dict)
	print "Average no of isoforms per gene = ", average
	max = sorted(isoform_number, reverse=True)
	print "Maximum no of isoforms per gene = ", max[0]
	min = sorted(isoform_number)
	print "Minimum no of isoforms per gene = ", min[0]

#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
	
	infiles = list_folder(args.infolder)
	print "Number of infiles (.fasta) =", len(infiles)

	for infile in infiles:
		print os.path.basename(infile)
		isoform_dict = read_transcripts(infile)
		print "Number of transcripts = ", len(isoform_dict)
		gene_dict = get_isoform_gene(isoform_dict)
		print "Number of genes = ", len(gene_dict)
		average(gene_dict)


if __name__ == '__main__':
	main()