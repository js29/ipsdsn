#!/usr/bin/env python
import os
import argparse
import fileinput
import subprocess

parser = argparse.ArgumentParser(description = "Merge multiple BAM/CRAM files into one BAM file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--indir", help = "Base directory for the input files")
parser.add_argument("--insuffix", default = ".bam")
parser.add_argument("--outsuffix", default = ".bam")
parser.add_argument("--outdir", help = "Base directory for the output file")
args = parser.parse_args()

# Lines should have two tab-separated fields; the first is the output sample name,
# and the second field is a comma-separated list of sample names to merge.
for line in fileinput.input("-"):
	fields = line.rstrip().split("\t")
	sample_name = fields[0]
	input_names = fields[1].split(",")

	#Construct file names from ids
	outdir = os.path.join(args.outdir, sample_name)
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	output_file = os.path.join(outdir, sample_name + args.outsuffix)
	input_files = [os.path.join(args.indir, in_name, in_name + args.insuffix) for in_name in input_names]
	
	command = " ".join(["samtools merge", output_file] + input_files)
	print(command)
	subprocess.call(['bash','-c',command])

