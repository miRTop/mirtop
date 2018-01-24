import traceback
import os.path as op
import os
import re
import shutil
import pandas as pd
import pysam

from mirtop.libs import do
from mirtop.libs.utils import file_exists
import mirtop.libs.logger as mylog

logger = mylog.getLogger(__name__)


def convert_gff_counts(args):
	""" Reads a GFF file, produces output file containing Expression counts of 
		the format : UID miRNA Variant Sample1 Sample2 ... Sample N 

	"""
	
	print "INFO Reading GFF file ", args.gff
	print "INFO Writing TSV file to directory ", args.out
	
	gff_file = open(args.gff,'r')
	out_file = op.join(args.out,"gff_tsv_expression_counts.tsv")
	#out_file = open(args.out+"gff_tsv_expression_counts.tsv",'w')
	out_file = open(out_file,'w')

	
	for samples_line in gff_file:
		if samples_line.startswith("## COLDATA:"):
			samples = "\t".join(samples_line.strip().split(" ")[2].split(","))
			header = "\t".join(['UID','miRNA','Variant', samples,"\n"])
			out_file.write(header)
			break 

	for mirna_line in gff_file:
		mirna_values = mirna_line.split("\t")[8].split(";")
		UID = mirna_values[1].split(" ")[2]
		mirna = mirna_values[2].split(" ")[2]
		variant = mirna_values[4].split(" ")[2]
		expression = "\t".join(mirna_values[6].split(" ")[2].split(","))
		summary = "\t".join([UID,mirna,variant,expression,"\n"])
		out_file.write(summary)
		
	gff_file.close()
	out_file.close()
	

#convert_gff_counts("/Users/shruthi/Desktop/mirtop/data/examples/gff/2samples.gff","/Users/shruthi/Desktop/outfile.txt")

	
