#!/usr/bin/env python

import argparse
import gzip

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Add GT field in VCF file made by Strelka2")
    parser.add_argument("-v","--vcf",help="Enter a vcf file made by Strelka2")
    parser.add_argument("-t","--type",help = "Enter the type of variants contained in vcf file (snps or indels)")
    args = parser.parse_args()

    def parse_info(info,type,ref,alt):
        splittedInfo = info.split(";")
        for field in splittedInfo:
            if field.startswith("NT="):
                normalGenotype = field.split("=")[1]
                normalGT = genotype_converter(normalGenotype)
            if field.startswith("SGT="):
                sgtValue = field.split("=")[1]
                tumorGenotype = sgtValue.split("->")[1]
                tumorGT = genotype_converter_type(type,tumorGenotype,ref,alt)
        return normalGT,tumorGT

    def genotype_converter(genotype):
        if genotype == "ref":
            return "0/0"
        elif genotype == "het":
            return "0/1"
        elif genotype == "hom":
            return "1/1"
        else:
            return "./."

    def genotype_converter_type(type,genotype,ref,alt):
        if type == "indels":
            return genotype_converter(genotype)
        if type == "snvs" or type == "snps":
            if genotype[0]==genotype[1]:
                if genotype[0] == alt:
                    return "1/1"
                if genotype[0] == ref:
                    return "0/0"
            else:
                if genotype[0] == ref and genotype[1] == alt:
                    return "0/1"
                if genotype[0] == alt and genotype[1] == ref:
                    return "1/0"


    vcf = args.vcf
    if vcf.endswith(".gz"):
        inputVCF = gzip.open(args.vcf)
    else:
        inputVCF = open(args.vcf)
    GT_header = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\""
    for line in inputVCF:
        if not line.startswith("#"):
            splittedLine = line.strip().split("\t")
            ref,alt,info,format,normal,tumor = splittedLine[3],splittedLine[4],splittedLine[7],splittedLine[8],splittedLine[9],splittedLine[10]
            typeVariants = args.type
            GTs_tuple = parse_info(info,typeVariants,ref,alt)
            updated_format = "GT:"+format
            updated_normal = GTs_tuple[0]+":"+normal
            updated_tumor = GTs_tuple[1]+":"+tumor
            updated_field_list = [updated_format,updated_normal,updated_tumor]
            new_line_list = splittedLine[0:8] + updated_field_list
            new_line = "\t".join(new_line_list)
            print new_line
        else:
            if args.type == "snvs" or args.type == "snps":
                if line.startswith("##FORMAT=<ID=AU"):
                    print GT_header
                    print line,
                else:
                    print line,
            if args.type == "indels":
                if line.startswith("##FORMAT=<ID=BCN50"):
                    print GT_header
                    print line,
                else:
                    print line,

