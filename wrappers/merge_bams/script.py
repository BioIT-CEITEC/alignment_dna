######################################
# wrapper for rule: merge_bams
######################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

shell.executable("/bin/bash")

print(snakemake.input)

if isinstance(snakemake.input,list) and len(snakemake.input) > 1:

    f = open(snakemake.log.run, 'a+')
    f.write("\n##\n## RULE: merge_bams \n##\n")
    f.close()

    command = "sambamba merge -t "+ str(snakemake.threads) +\
                   " " + snakemake.output.bam + ".temp.bam" +\
                   " " + " ".join(snakemake.input) +\
                   " >> " + snakemake.log.run + " 2>&1 "

    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

    command = " samtools view -H " + snakemake.output.bam + ".temp.bam" +\
              " | sed 's/SM:[^\\]*/SM:" + snakemake.params.sample_name + "/' | samtools reheader /dev/stdin " +\
              snakemake.output.bam + ".temp.bam >> " + snakemake.output.bam

    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

    command = "samtools index " + snakemake.output.bam

    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

    shell("rm " + snakemake.output.bam + ".temp.bam")
    shell("rm " + snakemake.output.bam + ".temp.bam.bai")

else:

    shell("mv -T " + snakemake.input + " " + snakemake.output.bam)
    shell("mv -T " + snakemake.input + ".bai " + snakemake.output.bam + ".bai")
