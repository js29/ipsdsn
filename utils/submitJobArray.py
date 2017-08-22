#!/usr/bin/env python
import sys
import os
import argparse
import fileinput

parser = argparse.ArgumentParser(description = "Submit jobs to the farm.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-j", "--jobname", help = "Name of the job")
parser.add_argument("--arraymaxjobs", type=int, help = "Max simultaneous array jobs")
parser.add_argument("--startat", type=int, default=1, help = "Array index to start at")
parser.add_argument("-c", "--command", help = "Excact command to be submitted.")
parser.add_argument("-m", "--MEM", help = "Amount of memory required.")
parser.add_argument("-o", "--farmout", help = "Folder for FARM log files", default = "FarmOut")
parser.add_argument("-n", "--ncores", help = "Number of cores to use", default = "1")
parser.add_argument("-q", "--queue", help = "Queue for submitting the jobs into.", default = "normal")
args = parser.parse_args()

memory_string = ' -R"span[hosts=1] select[mem>' + args.MEM + '] rusage[mem=' + args.MEM + ']" -M ' +args.MEM

arrayinputDir = os.path.join("arrayinput_" + args.jobname)
if not os.path.exists(arrayinputDir):
    os.makedirs(arrayinputDir)

lines = lines = sys.stdin.readlines()

for i,line in enumerate(lines):
	#line = line.rstrip()
	inputFname = os.path.join(arrayinputDir, 'input.' + str(i+args.startat))
	with open(inputFname, 'w') as f:
		f.write(line)

os.chdir(arrayinputDir)

if (args.farmout == "FarmOut"):
	args.farmout = "../FarmOut"
if not os.path.exists(args.farmout):
    os.makedirs(args.farmout)

jobname = args.jobname + "[{0}-{1}]".format(args.startat, args.startat+len(lines)-1)
output_file = os.path.join(args.farmout, jobname + '.%I.txt')

if args.arraymaxjobs:
	jobname = jobname + "%{0}".format(args.arraymaxjobs)

command = "".join(["bsub -G team170 -o ", output_file, " -q " + args.queue, " -n " + args.ncores, memory_string, " -J ", jobname, " -i ", "input.%I ", args.command])
print(command) 
os.system(command)

