import sys
import os
import argparse
import fileinput
import subprocess

parser = argparse.ArgumentParser(description = "Add donor id as RG field into the BAM file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--indir", help = "Directory of the input BAM files.")
parser.add_argument("--outdir", help = "Directory of the output BAM files.")
parser.add_argument("--insuffix", help = "Suffix of the input bam file.", default = ".bam")
parser.add_argument("--outsuffix", help = "Suffix of the output bam file.", default = ".RG.bam")
parser.add_argument("--SORT_ORDER", help = "Suffix of the output bam file.", default = "unsorted") # unsorted, queryname, coordinate, duplicate, unknown
parser.add_argument("--execute", help = "Execute the script", default = "False")
args = parser.parse_args()

#Construct file names
for line in fileinput.input("-"):
		line = line.rstrip()
		fields = line.split("\t")
		sample_name = fields[0]
		path_in = os.path.join(args.indir, sample_name, sample_name + args.insuffix)
		path_out = os.path.join(args.outdir, sample_name, sample_name + args.outsuffix)
		picard_path = "/software/java/bin/java -jar -Xmx900m /nfs/users/nfs_k/ka8/software/picard-tools-1.134/picard.jar"
		command = " ".join([picard_path, "AddOrReplaceReadGroups", "I="+path_in, "O="+path_out , "RGSM="+sample_name, "RGLB=1", "RGPU=1", "RGPL=Illumina"])
		if args.SORT_ORDER != "unsorted":
			command = command + " SORT_ORDER=" + args.SORT_ORDER
		print(command)
		if args.execute == "True":
			subprocess.call(['bash','-c',command])