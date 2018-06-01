#!/usr/bin/env python

import os
from subprocess import call
import argparse
import sys

def checkTumorControlMatches(txt):
    if os.path.exists(txt):
        with open(txt,'r') as f:
            for line in f:
                tumorSample = line.strip().split(" ")[0]
                controlSample = line.strip().split(" ")[1]
                if os.path.isdir("{}/{}".format(os.getcwd(),tumorSample)) and os.path.isdir("{}/{}".format(os.getcwd(),controlSample)):
                    pass
                else:
                    return "Error: {} or {} not found".format(tumorSample,controlSample)
    else:
        return "Error: {} file not found".format(txt)
    return True

def main():
    parser = argparse.ArgumentParser(description="Check fastq files in directory")
    parser.add_argument('-s', '--sconsdir', action='store', type=str, help="Path to directory containing SConstruct and settings.py files", default=".")
    parser.add_argument('-p', '--pairedend', action='store', type=bool, help="Paired-end or Single-end analysis (paired-end default)", default=True)
    parser.add_argument('-j', '--processors', action='store', type=int, help="Number of processors to be used", default=1)
    args = parser.parse_args()
    samples_directories = []
    empty_directories = []
    for f in os.listdir(os.getcwd()):
        if os.path.isdir(f):
            os.chdir(f)
            if args.pairedend == True:
                if not (os.path.isfile("1.fastq.gz") and os.path.isfile("2.fastq.gz")):
                    if not f == "Variants":
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
            os.system("scons -j {} -f {}/SConstruct".format(args.processors,args.sconsdir))
            os.chdir("..")
        goodSampleMatches = checkTumorControlMatches("tumor_control_samples.txt")
        if goodSampleMatches == True:
            os.system("mkdir Variants")
            os.system("mkdir Variants/Mutect Variants/Mutect2 Variants/Strelka2 Variants/VarScan")
            pairFile = open("tumor_control_samples.txt")
            for pair in pairFile:
                tumor,normal = pair.split(" ")[0],pair.split(" ")[1]
                pairName = tumor+"_"+normal.strip()
                os.system("mkdir Variants/Mutect/{0} Variants/Mutect2/{0} Variants/Strelka2/{0} Variants/VarScan/{0}".format(pairName))
                os.system("scons -j {} --debug=explain -f {}/Scons_variant_calling tumor={} normal={}".format(args.processors,args.sconsdir,tumor,normal))
        else:
            sys.stderr.write(goodSampleMatches+"\n")
            sys.exit(1)
        sys.exit(0)


if __name__ == "__main__":
    main()
