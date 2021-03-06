import os

#### Quick documentation ####
##Files needed by pipeline:
## - GRChxx.xx.fa
## - GRChxx.xx.fa.fai
## - GRChxx.xx.dict
## - dbSNP_vxxx_xxxx_noCHR.vcf - use remove_chr.py on vcf downloaded from dbSNP FTP site 
## - Bed file containing exome captured regions - *.bed
## - config.py - python file containing variables to input in the pipeline. All variables must be declared before
#running the pipeline
##Program needed by pipeline:
## - bwa - for alignment
## - picard.jar - for bam processing
## - GenomeAnalysisTK.jar - for indel realignment
## - gatk - GATK4 for read filtering and base recalibration
## - Bedtools - for coverage statistics computation

#########################
##### SCons Settings ####
#########################

vars = Variables("/working/configuration.py")

vars.Add('reference', 'The path to the directory containing the reference file',"reference.fa")
vars.Add('processors', 'Number of CPUs to be used', "2")
vars.Add('picard','The path to picard.jar',"/tmp/picard.jar")
vars.Add('gatk4','The path to gatk4',"/tmp/gatk4/gatk")
vars.Add('gatk3','The path to gatk3',"/tmp/gatk3/GenomeAnalysisTK.jar")
vars.Add('dbsnpVCF','The path to dbSNP VCF',"All_20180418.vcf.gz")
vars.Add('exomeRegions',"The path to bed file containing exome regions","exome_regions.bed")

env = Environment(ENV = os.environ, SHELL = '/bin/bash', variables = vars)
env.AppendENVPath('PATH', os.getcwd())
Decider('timestamp-newer')

####################
##### Arguments ####
####################

reference = ARGUMENTS.get("reference", env["reference"])

processors = ARGUMENTS.get("processors", env["processors"])

picard = ARGUMENTS.get("picard",env["picard"])

gatk4 = ARGUMENTS.get("gatk4",env["gatk4"])

gatk3 = ARGUMENTS.get("gatk3",env["gatk3"])

dbsnpVCF = ARGUMENTS.get("dbsnpVCF",env["dbsnpVCF"])

exomeRegions = ARGUMENTS.get("exomeRegions",env["exomeRegions"])

sampleName = os.path.basename(os.getcwd())

##################################################################
##### Alignment of the reads with bwa and sorting with picard ####
##################################################################

bwaCMD = "bwa mem -M -R \"@RG\\tID:{}\\tLB:{}\\tSM:{}\\tPL:ILLUMINA\"  -t {} /annotations/{}".format(sampleName,"exome",sampleName, processors, reference) + " <(zcat ${SOURCES[0]}) <(zcat ${SOURCES[1]}) | "

sortSamCMD = "java -jar {} SortSam INPUT=/dev/stdin OUTPUT=$TARGET SORT_ORDER=coordinate CREATE_INDEX=true".format(picard)

buildBamCMD = bwaCMD + sortSamCMD

bam = env.Command(["01_{}.bam".format(sampleName)],["1.fastq.gz","2.fastq.gz"],buildBamCMD)

################################################
##### Removal of PCR duplicates with Picard ####
################################################

pcrRemovalCMD = "java -Xmx4g -jar {} MarkDuplicates I=$SOURCE O=$TARGET M=metrics.txt REMOVE_DUPLICATES=true CREATE_INDEX=true".format(picard)
pcrRemoval = env.Command(["02_mapping-rmdup.bam"], [bam], pcrRemovalCMD)

######################################################################################################
##### Remove low mapping quality reads: bad CIGAR, unmapped reads, not primary aligned, ##############
##### failing vendor quality check, duplicated and mapping quality unavailable          ##############
######################################################################################################

filteringBamCMD = "{} PrintReads -R /annotations/{} -I $SOURCE -O $TARGET -RF MappingQualityNotZeroReadFilter -RF GoodCigarReadFilter -RF MappedReadFilter -RF PrimaryLineReadFilter -RF PassesVendorQualityCheckReadFilter -RF MappingQualityAvailableReadFilter -RF MateOnSameContigOrNoMappedMateReadFilter -RF PairedReadFilter".format(gatk4,reference)
filteringBam = env.Command(["03_mapping-rmdup-cleaned.bam"],[pcrRemoval],filteringBamCMD)

#########################################
#### Local realignment around indels ####
#########################################

#Table of putative indels

putativeIndelsTableCMD = "java -Xmx4g -jar {} -T RealignerTargetCreator -R /annotations/{} -o $TARGET  -I $SOURCE -nt {}".format(gatk3,reference,processors)
putativeIndelsTable = env.Command(["04_realigning.intervals"], [filteringBam], putativeIndelsTableCMD)

#Local realignment around indels

indelRealignerCmd = "java -Xmx4g -jar {} -R /annotations/{}".format(gatk3, reference) + " -I ${SOURCES[0]} -T IndelRealigner -targetIntervals ${SOURCES[1]} -o $TARGET"
indelRealignment = env.Command(["05_realigned.bam"], [filteringBam, putativeIndelsTable],indelRealignerCmd)

######################################
#### Quality score recalibration #####
######################################

recalibrationTableCMD = "{} BaseRecalibrator -I $SOURCE -R /annotations/{}  --known-sites /annotations/{} -O $TARGET".format(gatk4,reference,dbsnpVCF)
recalibrationTable = env.Command(["06_{}.grp".format(sampleName)], [indelRealignment], recalibrationTableCMD)

renderReadsCMD = "{} ApplyBQSR -R /annotations/{}".format(gatk4,reference) + " -I ${SOURCES[0]} -bqsr ${SOURCES[1]} -O $TARGET"
renderReads = env.Command(["07_recalibrated.bam"], [indelRealignment, recalibrationTable], renderReadsCMD)

######################
##### Statistics #####
######################

coverageHistCMD = "bedtools coverage -hist -abam $SOURCE -b /annotations/{}".format(exomeRegions) + " | grep ^all > $TARGET"
coverageHist = env.Command(["08_{}-coverage-hist.txt".format(sampleName)], [renderReads], coverageHistCMD)
