######################################
# wrapper for rule: alignment_DNA
######################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell


f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: alignment_DNA \n##\n")
f.close()

TOOL = "bwa"
SAMTOOLS = "samtools"
bwa_ref_prefix = re.sub(".bwt$","",snakemake.input.ref)

shell.executable("/bin/bash")

version = str(subprocess.Popen(TOOL+" 2>&1 | grep \"[Vv]ersion\" | cut -f 2 -d \" \"", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## VERSION: BWA "+version+"\n")
f.close()

version = str(subprocess.Popen(SAMTOOLS+" --version 2>&1 | grep \"samtools\" ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

if hasattr(snakemake.input, 'fastq2'):
    fastq2_file = snakemake.input.fastq2 + " "
else:
    fastq2_file = ""

command = TOOL+" mem -t "+str(snakemake.threads)+" -R '@RG\\tID:"+snakemake.params.lib_name+"_"+snakemake.wildcards.sample+"\\tSM:"+snakemake.wildcards.sample+"\\tPL:illumina' \
    -v 1 "+bwa_ref_prefix+" " + " ".join(snakemake.input.fastqs) + " 2>> "+snakemake.log.run+" | " \
    +SAMTOOLS+" sort -@ "+str(snakemake.threads)+" -o "+snakemake.output.bam+" /dev/stdin 2>> "+snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = SAMTOOLS+" index -b "+snakemake.output.bam+" 2>> "+snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
