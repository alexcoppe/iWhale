# Variant callers to be used, separated by ,
variantCallers = "mutect,mutect2,strelka2,varscan"

# Computational resources
processors = 3

# Human reference
reference = "reference.fa"

# Cosmic database in VCF format
cosmic = "CosmicCodingMuts_v82.vcf"

# dbSNP database in VCF format
#dbsnpVCF = "dbSNP_v150_20170710_noCHR.vcf.gz"
dbsnpVCF = "All_20180418.vcf.gz"

# ClinVar in VCF compressed format
clinvar = "clinvar_20190311.vcf.gz"

# Exome regions in bed format
exomeRegions = "exome_regions.bed"

# Genome aggregation database gnomAD in VCF format
gnomadVCF = "af-only-gnomad.raw.sites.b37.vcf.gz"

# Computation of predicted contamination in tumor samples by CalculateContamination from GATK4 (y/n)
contamination = "y"
# Predicted contamination in tumor samples maxPopulationAlleleFrequency parameter
maxPopulationAlleleFrequency = 0.2
# Predicted contamination in tumor samples minPopulationAlleleFrequency parameter
minPopulationAlleleFrequency = 0

# Options to run Java
javaOptions = 4

# snpeff directory location
snpeffDir = "/annotations"
# snpeff Genome Version
snpeffGenomeVersion = "GRCh37.75"


#### MuTect2 parameters ##############
mutect2Parameters = "--dont-use-soft-clipped-bases --output-mode EMIT_VARIANTS_ONLY"


### VarScan parameters ###############
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
