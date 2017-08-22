import sys
import os
import argparse
import fileinput
import subprocess
#This avoids the Broken pipe error when output is piped into head
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 

parser = argparse.ArgumentParser(description = "Runs bam2junc.sh for each sample piped to stdin.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--indir", required=True, help = "Directory of the input files.")
parser.add_argument("--insuffix", required=True, help = "Suffix of the input BAM file.")
parser.add_argument("--juncfilelist", required=True, help = "Text file to add path of new junction file to.")
args = parser.parse_args()

#Load sample names from disk
for line in fileinput.input("-"):
	sample_name = line.rstrip()
	path_in = os.path.join(args.indir, sample_name, sample_name + args.insuffix)
	path_out = path_in + ".junc"
	command = " ".join(["bam2junc.sh", path_in, path_out, "; echo", path_out, ">>", args.juncfilelist])
	print(command)
	subprocess.call(['bash', '-c', command])

