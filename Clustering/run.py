# import os
# import sys
# import inspect
import argparse

# currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
# sys.path.insert(0, parentdir)

import runClustering

parser = argparse.ArgumentParser()
parser.add_argument("--filename", type=str,
                    help="Root file to run clustering on")
parser.add_argument("--filenamePhoton", type=str,
                    help="Root file to run clustering on")
parser.add_argument("--output", type=str, help="outfile name")
parser.add_argument("--min", type=float,
                    help="Emin direct damage, default 5eV")
parser.add_argument("--max", type=float,
                    help="Emax direct damage, default 37.5eV")
parser.add_argument("--indirect", type=float,
                    help="Indirect damage percentage, default 0.405")
parser.add_argument("--sugar", type=str, help="sugar geometry file")
parser.add_argument("--sepR", type=bool, help="save to separate files for R")


args = parser.parse_args()

if args.filename:
    filename = args.filename
else:
    print("Missing filename")

if args.filenamePhoton:
    filenamePhoton = args.filenamePhoton
else:
    filenamePhoton = False

if args.output:
    outputFilename = args.output
else:
    print("Missing output filename")

if args.min:
    fEMinDamage = args.min
else:
    fEMinDamage = 5

if args.max:
    fEMaxDamage = args.max
else:
    fEMaxDamage = 37.5

if args.indirect:
    probIndirect = args.indirect
else:
    probIndirect = 0.405

if args.sugar:
    sugarPosFilename = args.sugar
if args.sepR:
    separate_r = args.sepR
else:
    print("Missing sugar geometry file")


runClustering.runClustering(filename, outputFilename, fEMinDamage, fEMaxDamage,
                            probIndirect, sugarPosFilename, filenamePhoton=filenamePhoton, separate_r=separate_r)
