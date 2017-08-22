#!/usr/bin/env python
import sys
import os
import argparse
import fileinput
import subprocess

parser = argparse.ArgumentParser(description = "Submit jobs to the farm.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-j", "--jobname", required=True, help = "Name of the job")
parser.add_argument("--addjobid", action='store_true', help = "Adds first field of input to job ID")
parser.add_argument("-c", "--command", help = "Exact command to be submitted.")
parser.add_argument("-m", "--MEM", type=str, default="100", help = "Amount of memory required.")
parser.add_argument("-o", "--farmout", help = "Folder for FARM log files", default = "FarmOut")
parser.add_argument("-n", "--ncores", help = "Number of cores to use", default = "1")
parser.add_argument("-q", "--queue", help = "Queue for submitting the jobs into.", default = "normal")
parser.add_argument("-b", "--blocking", action='store_true', help = "Wait for job to complete before returning.")
args = parser.parse_args()

memory_string = ' -R"span[hosts=1] select[mem>' + args.MEM +'] rusage[mem=' + args.MEM +']" -M '+args.MEM

if not os.path.exists(args.farmout):
    os.makedirs(args.farmout)

bsubcmd = "bsub"
bsub_bg = ""
if args.blocking:
	bsubcmd = "bsub -K"
	bsub_bg = " &"

# If no stdin, just run the command
if sys.stdin.isatty():
	output_file = os.path.join(args.farmout, args.jobname + '.%J.txt')
	command = "".join([bsubcmd, " -G team170 -o ", output_file, " -q " + args.queue, " -n " + args.ncores, memory_string, " -J ", args.jobname, " \"", args.command, "\"", bsub_bg])
	print(command) 
	#os.system(command)
	subprocess.call(command, shell=True, executable="/usr/local/bin/bash")
else:
	# Run one command for each line on stdin
	for line in fileinput.input("-"):
		line = line.rstrip()
	
		jobname = args.jobname
		if args.addjobid:
			fields = line.split("\t")
			jobname = jobname + "." + fields[0]
		
		output_file = os.path.join(args.farmout, jobname + '.%J.txt')
		
		command = "".join([bsubcmd, " -G team170 -o ", output_file, " -q " + args.queue, " -n " + args.ncores, memory_string, " -J ", jobname, " \"echo '", line,"'", " | ", args.command, "\"", bsub_bg])
		print(command) 
		#os.system(command)
		subprocess.call(command, shell=True, executable="/usr/local/bin/bash")

