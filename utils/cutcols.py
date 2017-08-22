#!/usr/bin/env python
import argparse
import sys
import os
import os.path
import subprocess

parser = argparse.ArgumentParser(description="Cuts columns from a file by name, assuming that the first line is a header.")
parser.add_argument("--file", type=file, required=True, metavar='FILE', help="input file name")
parser.add_argument("--colsFile", type=file, metavar='FILE', help="file containing column names (one per line)")
parser.add_argument("--cols", type=str, metavar='colname1,colname2,...', help="comma-separated list of column names")
parser.add_argument("--sep", type=str, metavar='SEP[\\t]', default="\t", help="field separator (default: tab)")
parser.add_argument("--ordered", action='store_true', help="should cols be output in same order as specified? (slower)")
parser.add_argument("--debug", action='store_true')

args = parser.parse_args()

def IsInteger(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False
   
cols = []
if args.cols != None:
	cols.extend(args.cols.split(','))
if args.colsFile != None:
	cols.extend(args.colsFile.read().splitlines())
if len(cols) <= 0:
	sys.stderr.write("No column names provided. You must specify either --cols or --colsFile\n")
	exit(1)

if args.debug:
	sys.stderr.write(str(cols) + "\n")

headerStr = args.file.readline()
headerCols = headerStr.strip().split(args.sep)
headerColDict = {}
for i, colname in enumerate(headerCols):
	headerColDict[colname] = i
	if args.debug:
		sys.stderr.write(":".join([str(i+1), colname]) + "\n")

colIndexes = []
for colname in cols:
	if IsInteger(colname):
		headerColDict[colname] = int(colname)-1
	if (not colname in headerColDict):
		sys.stderr.write("Column '" + colname + "' not found in header\n")
		exit(1)
	colIndexes.append(headerColDict[colname])

if (not args.ordered):
	cutCmd = "cut -f " + ",".join(str(x+1) for x in colIndexes) + " " + args.file.name
else:
	cutCmd = "paste"
	for x in colIndexes:
		cutCmd += " <(cut -f {0} {1})".format(x+1, args.file.name)

if args.debug:
	sys.stderr.write("CMD:" + cutCmd + "\n")
subprocess.call(cutCmd, shell=True, executable='/bin/bash')
