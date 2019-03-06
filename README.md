# :whale: iWhale
### Dockerized Whole Exome Sequencing (WES) pipeline.

A pipeline for Whole Exome Sequencing (WES) variants identification in mathced tumor-normal samples. It runs into a [Docker](https://www.docker.com) container, ready to be downloaded and used on any operating system 
supported by [Docker](https://www.docker.com/).

All the steps of the pipeline and their dependencies are controlled by [SCons](https://scons.org/) so that in case of any stop, like killing by error or even shutting down the computer, it will automatically resume the analysis from the last run process.

Three variant calling softwares are used by the pipeline: [Mutect2](https://software.broadinstitute.org/gatk/gatk4) , [VarScan2](http://dkoboldt.github.io/varscan/), and [Strelka2](https://github.com/Illumina/strelka) and the user is allowed to choose which to use and change their default settings.


### Software versions currently used

| Program        | Description| Version |
| ------------- |:-------------| :-------------| 
[Burrows-Wheeler Aligner (BWA)](http://bio-bwa.sourceforge.net/) |BWA is a software package for mapping low-divergent sequences against a large reference genome.| 0.7.17 |
[Picard](https://broadinstitute.github.io/picard/) | A set of command line tools for manipulating high-throughput sequencing (HTS) data and formats.| 2.17.11 |
[GATK4](https://software.broadinstitute.org/gatk/gatk4) | GATK4 is the first and only open-source software package that covers all major variant classes for both germline and cancer genome analysis. | 4.0.6.0 |
[GATK3](https://gatkforums.broadinstitute.org/gatk/categories/gatk-guide) | A variety of tools with a primary focus on variant discovery and genotyping. | 3.8-1 |
[Strelka2](https://github.com/Illumina/strelka) | Strelka2 is a fast and accurate small variant caller optimized for analysis of germline variation in small cohorts and somatic variation in tumor/normal sample pairs. | 2.9.2 |
[VarScan](http://dkoboldt.github.io/varscan/)| VarScan is a platform-independent software tool developed at the Genome Institute at Washington University to detect variants in NGS data. | 2.4.2 |
[SnpEff](http://snpeff.sourceforge.net/) | Genomic variant annotations and functional effect prediction toolbox. | 4_3t |
[bedtools](https://bedtools.readthedocs.io/en/latest/)|Collectively, the bedtools utilities are a swiss-army knife of tools for a wide-range of genomics analysis tasks. | 2.17.0 |
