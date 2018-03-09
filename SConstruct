import os

vars = Variables("config.py")
vars.Add('referenceDir', 'The path to the directory containing the reference file', "~/annotation")
vars.Add('processors', 'Number of CPUs to be used', "2")

env = Environment(ENV = os.environ, SHELL = '/bin/bash', variables = vars)
env.AppendENVPath('PATH', os.getcwd())
Decider('timestamp-newer')

referenceDir = ARGUMENTS.get("referenceDir", env["referenceDir"])
if not referenceDir.endswith("/"): referenceDir = referenceDir + "/"

processors = ARGUMENTS.get("processors", env["processors"])

sampleName = os.path.basename(os.getcwd())

bwaCMD = "bwa mem -M -R \"@RG\\tID:{}\\tLB:{}\\tSM:{}\\tPL:ILLUMINA\"  -t {} {}reference.fa ".format(sampleName,
"exome", sampleName, processors, referenceDir) + "<(zcat 1.fastq.gz) <(zcat 2.fastq.gz)"
bwa = env.Command(["01_{}.bam".format(sampleName)],[],bwaCMD)
