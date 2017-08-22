#!/usr/bin/env python
import os
import sys
import argparse
import fileinput
import subprocess

parser = argparse.ArgumentParser(description = "Count the number of fragments in BAM file that overlap features in GTF file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--sampleDir", help = "Path to sample directory.")
parser.add_argument("--outputDir", help = "Path to output directory.")
parser.add_argument("--bamSuffix", help = "Full suffix of the bam file.", default = ".Aligned.out.bam")
parser.add_argument("--countsSuffix", help = "Suffix of the counts file.", default = ".counts.txt")
parser.add_argument("--gtf", help = "Path to gene annotations in GTF format")
parser.add_argument("--g", help = "Feature to group by (e.g. exon_id) - default is gene_id")
parser.add_argument("-F", help = "isGTFAnnotationFile - should be 'GTF' or 'SAF'")
parser.add_argument("-f", action = 'store_true', help = "useMetaFeatures")
parser.add_argument("--execute", help = "If True then executes the command, otherwise just prints it out.", default = "True")
parser.add_argument("--strand", help = "0 (unstranded); 1 (stranded); 2(reversely stranded)", default = "0")
parser.add_argument("--multimapping", help = "Count multimapping reads.", default = "False")
parser.add_argument("--unpaired", help = "BAM contains single-end reads.", default = "False")
parser.add_argument("--D", help = "Maximum insert size.", default = "2000")
parser.add_argument("--donotsort", help = "Do not sort the BAM file", default = "False")
parser.add_argument("--O", help = "Assign reads to all overlapping features.", default = "False")
parser.add_argument("--overwrite", help = "Overwrite counts output file if it exists.", default = "False")
args = parser.parse_args()

#Iterate over all ids
for line in fileinput.input("-"):
	line = line.rstrip()
	bam_file = os.path.join(args.sampleDir, line, line + args.bamSuffix)
	if args.outputDir is None:
		args.outputDir = args.sampleDir
	count_file = os.path.join(args.outputDir, line, line + args.countsSuffix)

	if args.overwrite == "False" and os.path.exists(count_file):
		print("".join(["File exists: ", count_file, " and overwrite=False. Exiting."]))
		exit(0)
		
	featureCounts_command = " ".join(["featureCounts -a", args.gtf, "-o", count_file, "-s", args.strand])
	if(args.multimapping == "True"):
		featureCounts_command = featureCounts_command + " -M"
	if(args.donotsort == "True"):
		featureCounts_command = featureCounts_command + " --donotsort"
	if(args.g is not None):
		featureCounts_command = featureCounts_command + " -g " + args.g
	if(args.O == "True"):
		featureCounts_command = featureCounts_command + " -O"
	if(args.F):
		featureCounts_command = featureCounts_command + " -F " + args.F
	if(args.f):
		featureCounts_command = featureCounts_command + " -f"
	if(args.unpaired == "False"):
		featureCounts_command = " ".join([featureCounts_command, "-p -C -D", args.D, "-d 25"])
	command = " ".join([featureCounts_command, bam_file])
	print(command)
	if (args.execute == "True"):
		retval = subprocess.call(['bash','-c',command])
		if (retval != 0):
			exit(retval)
