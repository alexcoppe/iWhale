#Computational resources

processors = 3

#reference = "/home/andrea/Andrea/Script/Docker/wes_analysis/test/genome/reference.fa"
reference = "reference.fa"


cosmic = "CosmicCodingMuts_v82.vcf"

dbsnpVCF = "dbSNP_v150_20170710_noCHR.vcf.gz"

# ClinVar vcf file
clinvar = "clinvar.vcf"

exomeRegions = "small.bed"
#exomeRegions = "SureSelect_Human_All_Exon_UTR_V5.bed"

gnomadVCF = "gnomAD_lite_docker.vcf"

contamination = "y"

javaOptions = 4


snpeffDir = "/annotations"
snpeffGenomeVersion = "GRCh37.75"

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

# Custom variant callers parameters

mutectParameters = ""
mutect2Parameters = ""
varScanParameters = ""
