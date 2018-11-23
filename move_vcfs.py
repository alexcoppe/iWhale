#!/usr/bin/env python

import argparse
import os
import shutil

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Moves VCFs files from Variants folder to VCF folder")
    parser.add_argument("-f","--folder",help = "Enter Variants folder path",required = True)
    parser.add_argument("-o","--output",help = "Enter destination folder path",required = True)
    args = parser.parse_args()

    VariantsDir = args.folder
    OutputDir = args.output

    outputPath = os.path.abspath(OutputDir)
    os.chdir(outputPath)
    destinationFolder = outputPath + "/VCF/"
    if not os.path.isdir(destinationFolder):
        os.mkdir("VCF")

    VariantsPath = os.path.abspath(VariantsDir)
    os.chdir(VariantsDir)
    for directory in os.listdir(VariantsPath):
        if os.path.isdir(directory) and directory in ["Mutect","Mutect2","Strelka2","VarScan"]:
            currentDir = os.getcwd()
            vcDir = currentDir + "/" + directory
            os.chdir(vcDir)
            for sample in os.listdir(vcDir):
                if os.path.isdir(sample):
                    os.chdir(sample)
                    fileList = []
                    for vcf in os.listdir(vcDir + "/" + sample):
                        if os.path.isfile(vcf) and vcf.startswith("0") and vcf.endswith(".vcf"):
                            fileList.append(vcf)
                    fileList.sort()
                    if directory in ["Mutect2","Strelka2","VarScan"]:
                        vcfToMoveList = fileList[-2:]
                        for vcfToMove in vcfToMoveList:
                            shutil.copy(os.getcwd() + "/" + vcfToMove,destinationFolder + vcfToMove)
                            newFileName = vcfToMove.replace(vcfToMove[0:2],sample)
                            os.rename(destinationFolder+vcfToMove,destinationFolder + newFileName)
                    else:
                        vcfToMove = fileList[-1]
                        shutil.copy(os.getcwd() + "/" + vcfToMove,destinationFolder + vcfToMove)
                        newFileName = vcfToMove.replace(vcfToMove[0:2],sample)
                        os.rename(destinationFolder+vcfToMove,destinationFolder + newFileName)
                    os.chdir(vcDir)
            os.chdir(VariantsPath)








