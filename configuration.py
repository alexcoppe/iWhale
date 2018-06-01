#Computational resources

processors = 3

#reference = "/home/andrea/Andrea/Script/Docker/wes_analysis/test/genome/reference.fa"
reference = "reference.fa"

#gatk4 = "/home/andrea/Andrea/Script/Docker/wes_analysis/test/gatk-4.0.2.1/gatk"

#mutect = "/home/andrea/Andrea/Script/Docker/wes_analysis/test/mutect-1.1.7.jar"

#cosmic = "/home/andrea/Andrea/Script/Docker/wes_analysis/test/annotation/Cosmic/CosmicCodingMuts_v82.vcf"
cosmic = "CosmicCodingMuts_v82.vcf"

#dbsnpVCF = "/home/andrea/Andrea/Script/Docker/wes_analysis/test/annotation/dbSNP_hg19/dbSNP_v150_20170710_noCHR.vcf.gz"
dbsnpVCF = "dbSNP_v150_20170710_noCHR.vcf.gz"

#exomeRegions = "/home/andrea/Andrea/Script/Docker/wes_analysis/test/annotation/bed/S07604514_Padded_noChr.bed"
exomeRegions = "small.bed"

gnomadVCF = "gnomAD_lite_docker.vcf"

contamination = "y"

javaOptions = 4

maxPopulationAlleleFrequency = 0.2

minPopulationAlleleFrequency = 0

#gnomadVCF = "/home/andrea/Andrea/Script/Docker/wes_analysis/test/annotation/gnomAD/gnomAD_lite_docker.vcf"

variantCallers = "mutect,mutect2,strelka2,varscan"

#### MuTect2 parameters ####
mutect2Parameters = "--dont-use-soft-clipped-bases --output-mode EMIT_VARIANTS_ONLY"

############## VarScan parameters ###############
# VarScan --min-coverage
varScanMinCoverage = 10
# VarScan --pvalue
varScanPvalue = 0.05
# Varscan --min-var-freq 0.02 
varScanMinVarFreq = 0.02
