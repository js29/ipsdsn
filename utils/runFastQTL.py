import sys
import os
import argparse
import fileinput
import subprocess

parser = argparse.ArgumentParser(description = "Test for association between genotype and phenotype.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--vcf", help = "Path to the genotype VCF file.")
parser.add_argument("--bed", help = "Path to the phenotype bed file.")
parser.add_argument("--cov", help = "Path to the covariates txt file.")
parser.add_argument("--include-covariates", help = "File listing covariates to include")
parser.add_argument("--window", help = "Size of the cis window.")
parser.add_argument("--out", help = "Prefix of the output files.")
parser.add_argument("--permute", help = "Number of permutations to perform.")
parser.add_argument("--threshold", help = "Filter out gene-snp pairs with p value above this threshold.")
parser.add_argument("--execute", help = "Execute the script", default = "False")
args = parser.parse_args()

#Check for required paramters
if args.out == None:
	sys.exit("--out is a required parameter.")
if args.vcf == None:
	sys.exit("--vcf is a required parameter.")
if args.bed == None:
	sys.exit("--bed is a required parameter.")
if args.window == None:
	sys.exit("--window is a required parameter.")

#Construct file names
for line in fileinput.input("-"):
	chunk = line.rstrip()
	chunks = chunk.split(" ")
	print(chunks)
	outfile = args.out + ".chunk_" + chunks[0] + "_" + chunks[1] + ".txt.gz"
	command = " ".join(["fastQTL --seed 37 --vcf", args.vcf, "--bed", args.bed, "--chunk", chunk, "--out", outfile, "--window", args.window])

	if (args.cov != None):
		command = " ".join([command, "--cov", args.cov])
		if (args.include_covariates != None):
			command = " ".join([command, "--include-covariates", args.include_covariates])
	
	if (args.threshold != None): 
		command  = " ".join([command, "--threshold", args.threshold])

	if (args.permute != None): 
		command  = " ".join([command, "--permute", args.permute])
	print(command)
	if args.execute == "True":
		subprocess.call(['bash','-c',command])