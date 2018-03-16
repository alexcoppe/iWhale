import os

vars = Variables("config.py")
#TOCHANGE
vars.Add('referenceDir', 'The path to the directory containing the reference file',
"/home/andrea/Andrea/Script/Docker/wes_analysis/test/genome") 
vars.Add('processors', 'Number of CPUs to be used', "2")
#TOCHANGE
vars.Add('picard','The path to picard.jar',"/home/andrea/Andrea/Tools/picard_2.12.0/picard.jar")

env = Environment(ENV = os.environ, SHELL = '/bin/bash', variables = vars)
env.AppendENVPath('PATH', os.getcwd())
Decider('timestamp-newer')

referenceDir = ARGUMENTS.get("referenceDir", env["referenceDir"])
if not referenceDir.endswith("/"): referenceDir = referenceDir + "/"

processors = ARGUMENTS.get("processors", env["processors"])

picard = ARGUMENTS.get("picard",env["picard"])

sampleName = os.path.basename(os.getcwd())

bwaCMD = "bwa mem -M -R \"@RG\\tID:{}\\tLB:{}\\tSM:{}\\tPL:ILLUMINA\"  -t {} {}reference.fa".format(sampleName,"exome",sampleName, processors, referenceDir) + " <(zcat 1.fastq.gz) <(zcat 2.fastq.gz) | "

sortSamCMD = "java -jar {} SortSam INPUT=/dev/stdin OUTPUT=$TARGET SORT_ORDER=coordinate CREATE_INDEX=true".format(picard)

buildBamCMD = bwaCMD + sortSamCMD

bwa = env.Command(["01_{}.bam".format(sampleName)],[],buildBamCMD)


