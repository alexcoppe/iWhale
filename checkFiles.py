#!/usr/bin/env python

import os
from subprocess import call
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description="Check fastq files in directory")
    parser.add_argument('-s', '--sconsdir', action='store', type=str, help="Path to directory containing SConstruct and settings.py files", default=".")
    parser.add_argument('-p', '--pairedend', action='store', type=bool, help="Paired-end or Single-end analysis (paired-end default)", default=True)
    args = parser.parse_args()
    samples_directories = []
    empty_directories = []
    print args.sconsdir
    for f in os.listdir(os.getcwd()):
        if os.path.isdir(f):
            os.chdir(f)
            if args.pairedend == True:
                if not (os.path.isfile("1.fastq.gz") and os.path.isfile("2.fastq.gz")):
                    empty_directories.append(f)
                else:
                    samples_directories.append(f)
                os.chdir("..")
            else:
                print "single"
    if len(empty_directories) > 0:
        nameEmpytDir = ",".join(empty_directories)
        sys.stderr.write("{} directories are empty or don't have paired fastq files\n".format(nameEmpytDir))
        sys.exit(1)
    else:
        for directory in samples_directories:
            os.chdir(directory)
            os.system("scons -f {}/SConstruct".format(args.sconsdir))
            os.chdir("..")
        sys.exit(0)



if __name__ == "__main__":
    main()
