import sys
import os
import argparse
import fileinput
import subprocess
import gzip

parser = argparse.ArgumentParser(description = "Add snp coordinates to fastQTL p-values file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--vcf", help = "Path to the genotype VCF file.")
parser.add_argument("--snpcol", default=2, help = "rsID column in the fastqtl output file.")
parser.add_argument("--fastqtl", help = "Path to the phenotype bed file.")
args = parser.parse_args()

vcf_file = gzip.open(args.vcf)
fastqtl_file = gzip.open(args.fastqtl)
args.snpcol = int(args.snpcol) - 1

variant_pos_dict = dict()
for line in vcf_file:
	if line[0] != "#":
		fields = line.split("\t")
		variant_pos_dict[fields[2]] = fields[0:2]

for line in fastqtl_file:
	line = line.rstrip()
	fields = line.split(" ")
	snp_id = fields[args.snpcol]
	if snp_id in variant_pos_dict:
		coords = variant_pos_dict[snp_id]
	else:
		coords = ["NA", "NA"]
	line = " ".join(fields + [coords[0], coords[1]])
	print(line)

