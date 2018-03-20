import os

#### Quick documentation ####
##Files needed by pipeline:
## - GRChxx.xx.fa
## - GRChxx.xx.fa.fai
## - GRChxx.xx.dict
## - dbSNP_vxxx_xxxx_noCHR.vcf - use remove_chr.py on vcf downloaded from dbSNP FTP site 
## - CosmicCodingMuts_vxx.vcf
## - Bed file containing exome captured regions - *.bed
## - config.py - python file containing variables to input in the pipeline. All variables must be declared before
#running the pipeline
##Program needed by pipeline:
## - bwa - for alignment
## - picard.jar - for bam processing
## - GenomeAnalysisTK.jar - for indel realignment
## - Bedtools - for coverage statistics computation

#########################
##### SCons Settings ####
#########################

vars = Variables("config.py")
#TOCHANGE
vars.Add('referenceDir', 'The path to the directory containing the reference file',
"/home/andrea/Andrea/Script/Docker/wes_analysis/test/genome") 
vars.Add('processors', 'Number of CPUs to be used', "2")
#TOCHANGE
vars.Add('picard','The path to picard.jar',"/home/andrea/Andrea/Tools/picard_2.12.0/picard.jar")
#TOCHANGE
vars.Add('gatk4','The path to gatk4',"/home/andrea/Andrea/Script/Docker/wes_analysis/test/gatk-4.0.2.1/gatk")


env = Environment(ENV = os.environ, SHELL = '/bin/bash', variables = vars)
env.AppendENVPath('PATH', os.getcwd())
Decider('timestamp-newer')

####################
##### Arguments ####
####################

referenceDir = ARGUMENTS.get("referenceDir", env["referenceDir"])
if not referenceDir.endswith("/"): referenceDir = referenceDir + "/"

processors = ARGUMENTS.get("processors", env["processors"])

picard = ARGUMENTS.get("picard",env["picard"])

gatk4 = ARGUMENTS.get("gatk4",env["gatk4"])

sampleName = os.path.basename(os.getcwd())

##################################################################
##### Alignment of the reads with bwa and sorting with picard ####
##################################################################

bwaCMD = "bwa mem -M -R \"@RG\\tID:{}\\tLB:{}\\tSM:{}\\tPL:ILLUMINA\"  -t {} {}reference.fa".format(sampleName,"exome",sampleName, processors, referenceDir) + " <(zcat 1.fastq.gz) <(zcat 2.fastq.gz) | "

sortSamCMD = "java -jar {} SortSam INPUT=/dev/stdin OUTPUT=$TARGET SORT_ORDER=coordinate CREATE_INDEX=true".format(picard)

buildBamCMD = bwaCMD + sortSamCMD

bam = env.Command(["01_{}.bam".format(sampleName)],[],buildBamCMD)

################################################
##### Removal of PCR duplicates with Picard ####
################################################

pcrRemovalCMD = "java -Xmx4g -jar {} MarkDuplicates I=$SOURCE O=$TARGET M=metrics.txt REMOVE_DUPLICATES=true CREATE_INDEX=true".format(picard)
pcrRemoval = env.Command(["02_mapping-rmdup.bam"], [bam], pcrRemovalCMD)

######################################################################################################
##### Remove low mapping quality reads: bad CIGAR, unmapped reads, not primary aligned, ##############
##### failing vendor quality check, duplicated and mapping quality unavailable          ##############
######################################################################################################

filteringBamCMD = "{} PrintReads -R {}reference.fa -I $SOURCE -O $TARGET -RF MappingQualityNotZeroReadFilter -RF GoodCigarReadFilter -RF MappedReadFilter -RF PrimaryLineReadFilter -RF PassesVendorQualityCheckReadFilter -RF MappingQualityAvailableReadFilter -RF MateOnSameContigOrNoMappedMateReadFilter -RF PairedReadFilter".format(gatk4,referenceDir)
filteringBam = env.Command(["03_mapping-rmdup-cleaned.bam"],[pcrRemoval],filteringBamCMD)


