# :whale: iWhale
### Dockerized Whole Exome Sequencing (WES) pipeline

A pipeline for Whole Exome Sequencing (WES) variants identification in mathced tumor-normal samples. It runs into a [Docker](https://www.docker.com) container, ready to be downloaded and used on any operating system 
supported by [Docker](https://www.docker.com/).

All the steps of the pipeline and their dependencies are controlled by [SCons](https://scons.org/) so that in case of any stop, like killing by error or even shutting down the computer, it will automatically resume the analysis from the last run process.

Three variant calling softwares are used by the pipeline: [Mutect2](https://software.broadinstitute.org/gatk/gatk4) , [VarScan2](http://dkoboldt.github.io/varscan/), and [Strelka2](https://github.com/Illumina/strelka) and the user is allowed to choose which to use and change their default settings.

# iWhale usage

The first thing to do is to download iWhale from [Docker Hub](https://hub.docker.com/):

```
docker pull alexcoppe/iWhale
```
The working directory has to contain the following elements:
- two directories for each sample (one for tumor and one for control sample) including the two fastq files
- a text file called **tumor_control_samples.txt**
- a python file called **configuration.py**

The working directory **must not contain** other directories except the ones of the samples indicated above

#### Sample directories structure
Each sample must be in its own directory containing the two paired-end gz-compressed fastq files. The files **must** be called **1.fastq.gz** and **2.fastq.gz** 

#### tumor_control_samples.txt file structure
This is a simple text file organized by two columns separated by tab: in the first column there are tumor directories names and in the second one the matched control directories names 

```
tumor_sample1 control_sample1
tumor_sample2 control_sample2
...
```

#### configuration.py file structure
This file is essential and can be empty. It can be used to set parameters of the tools used by iWhale. All the possible parameters that you can set are gathered and explained in this file: [configuration.py](https://raw.githubusercontent.com/alexcoppe/iWhale/master/configuration.py?token=AV00ooKqBg4p7Us1Lnsk6GLpgw3fL3mUks5cnOYCwA%3D%3D)


# Software versions currently used

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


# Databases currently used

iWhale uses databases and sequences indicated in the table below. Many of these sequence files should be indexed. We provide a bash script that do all the steps, from download to index.

| Sequences or Databases | Description| Version |
| ------------- |:-------------| :-------------| 
| [Human genome database](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/) |Feb. 2009 assembly of the human genome (hg19,GRCh37 Genome Reference Consortium Human Reference 37|  hg19 or GRCh37|
|[dbSNP](https://www.ncbi.nlm.nih.gov/snp)| dbSNP contains human single nucleotide variations, microsatellites, and small-scale insertions and deletions | All_20180418.vcf.gz|
|[gnomAD](https://gnomad.broadinstitute.org/)| The Genome Aggregation Database (gnomAD), developed by an international coalition of investigators, with the goal of aggregating and harmonizing both exome and genome sequencing data from a wide variety of large-scale sequencing projects |2018-05-22|
|[SnpEff GRCh37](http://snpeff.sourceforge.net/index.html)| SnpEff annotation for the human genome reference genome GRCh37)| GRCh37 |
|[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)| ClinVar aggregates information about genomic variation and its relationship to human health |20190311|


# Getting database files and indexing by yourself

All commands are launch from the directory containing the downloaded data. **Many of the command take a LOT OF TIME to conclude**.

### BWA indexing of Human genome.

It produces many files:
- reference.fa.amb
- reference.fa.ann
- reference.fa.bwt
- reference.fa.pac
- reference.fa.sa

This step takes a lot of time:

```
bwa index reference.fa
```

### Index of the FASTA file with human genome data for picard.

The produced files is:
 - reference.dict

```
java -jar ~/local/picard.jar CreateSequenceDictionary R=reference.fa O=reference.dict
```

### Creation of the reference.fa.fai index.

The produced files is:
 - reference.fa.fai

```
samtools faidx reference.fa
```
 
### dbSNP download:

```
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz 
```

### Removing chr from dbSNP downloaded file (from chr1 to 1)

```
gunzip All_20180418.vcf.gz
```


### Indexing dbSNP VCF with tabix.

Need to install tabix in your computer.
Produced file:
 - All_20180418.vcf.gz.tbi
 
 This step takes a lot of time
 
```
tabix -fp vcf  All_20180418.vcf.gz
```

### Download of gnomAD data

Download gnomAD data from [http://bioinfo5pilm46.mit.edu/software/GATK/resources/](http://bioinfo5pilm46.mit.edu/software/GATK/resources/)
```
wget http://bioinfo5pilm46.mit.edu/software/GATK/resources/af-only-gnomad.raw.sites.b37.vcf
wget http://bioinfo5pilm46.mit.edu/software/GATK/resources/af-only-gnomad.raw.sites.b37.vcf.tbi
```

### Download of ClinVar data

Download of ClinVar data from [ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/](ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/).
Obtained files:
 - clinvar_20190311.vcf.gz
 - clinvar_20190311.vcf.gz.tbi 
 
Their date portion could be different, if so, remember to put the right filename in the configuration.py file:

```
clinvar = "clinvar_20190311.vcf.gz"
```

