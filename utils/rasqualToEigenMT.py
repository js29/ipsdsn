import sys
import os
import gzip
import argparse
import fileinput
import subprocess
from scipy import stats

def main():
	parser = argparse.ArgumentParser(description = "Convert RASQUAL output into a format suitable for eigenMT. Extract relevant columns and convert chisq statistic into p-value", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--rasqualOut", type=file, help = "Path to the merged RASQUAL output file")
	args = parser.parse_args()

	rasqual_file = openg(args.rasqualOut)
	print "\t".join(["snps","gene","statistic","pvalue","FDR","beta"])
	for line in rasqual_file:
		line = line.rstrip()
		fields = line.split("\t")
		gene_id = fields[0]
		snp_id = fields[1]
		pi = fields[11]
		#Calculate p-value:
		chi_stat = float(fields[10])
		p_value = stats.chi2.sf(chi_stat, 1)
		snp = "\t".join([snp_id, gene_id, str(chi_stat), str(p_value), str(p_value), pi])
		if (snp_id != "SKIPPED"): #Ignore skipped genes
			print(snp)


def isgzfile(fname):
	namelen = len(fname)
	return fname[namelen-3:namelen] == ".gz"

def openg(f, mode=None):
	if type(f) is file:
		if isgzfile(f.name):
			return gzip.GzipFile(fileobj=f)
		return f
	elif type(f) is str:
		if isgzfile(f):
			if mode is None:
				mode = 'rb'
			return gzip.open(f, mode)
		if mode is None:
			mode = 'r'
		return open(f, mode)
	else:
		die("openg: unrecognized type for parameter 'f': {0}".format(type(f)))
 

main()
