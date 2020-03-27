# Variant callers to be used, separated by ,
variantCallers = "mutect2,strelka2,varscan"

# Computational resources
processors = 4

# Human reference
reference = "GRCh38.fa"

# Cosmic database in VCF format
cosmic = "CosmicCodingMuts.vcf.gz"

# dbSNP database in VCF format
#dbsnpVCF = "dbSNP_v150_20170710_noCHR.vcf.gz"
dbsnpVCF = "All_20180418.vcf.gz"

# ClinVar in VCF compressed format
clinvar = "clinvar_20200316.vcf.gz"

# Exome regions in bed format
exomeRegions = "exome_regions.bed"

# Genome aggregation database gnomAD in VCF format
gnomadVCF = "gnomad.exomes.r2.1.sites.grch38.vcf.gz"

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
snpeffGenomeVersion = "GRCh38.86"


#### MuTect2 parameters ##############
mutect2Parameters = "--dont-use-soft-clipped-bases --output-mode EMIT_VARIANTS_ONLY"


### VarScan parameters ###############
# VarScan --min-coverage
varScanMinCoverage = 10
# VarScan --pvalue
varScanPvalue = 0.05
# Varscan --min-var-freq 0.02
varScanMinVarFreq = 0.02


