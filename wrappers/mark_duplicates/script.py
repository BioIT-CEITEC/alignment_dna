######################################
# wrapper for rule: mark_duplicates
######################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell


f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: mark_duplicates \n##\n")
f.close()

shell.executable("/bin/bash")

if snakemake.params.primer_based != "yes":
    os.makedirs(os.path.dirname(snakemake.params.mtx),exist_ok=True)
    if snakemake.params.umi == "no":

        f = open(snakemake.log.run, 'at')
        version = str(subprocess.Popen("picard MarkDuplicates --version 2>&1",shell=True,stdout=subprocess.PIPE).communicate()[0], 'utf-8')
        f.write("## VERSION: Picard "+version+"\n")
        f.close()

        command = "picard MarkDuplicates INPUT="+snakemake.input.bam+" OUTPUT="+snakemake.output.bam+" METRICS_FILE="+snakemake.params.mtx+" REMOVE_DUPLICATES="+snakemake.params.rmDup+" \
            ASSUME_SORTED=true PROGRAM_RECORD_ID=null VALIDATION_STRINGENCY=LENIENT -Xmx"+str(snakemake.resources.mem)+"g 2>> "+snakemake.log.run+" "
        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)

        command = "samtools index "+snakemake.output.bam
        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)

        command = "rm "+snakemake.input.bam
        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)

        command = "rm "+snakemake.input.bam + ".bai"
        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)

    else:


        shell.executable("/bin/bash")

        version = str(subprocess.Popen("je markdupes --version 2>&1 | tail -n 1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
        f = open(snakemake.log.run, 'at')
        f.write("## VERSION: je-suite "+version+"\n")
        f.close()
        java_opts = "export _JAVA_OPTIONS='-Xmx"+str(snakemake.resources.mem)+"g'"

        command = java_opts + "&& je markdupes " + \
                  "INPUT=" + snakemake.input.bam +\
                  " OUTPUT="+snakemake.output.bam+\
                  " METRICS_FILE="+snakemake.params.mtx+\
                  " REMOVE_DUPLICATES="+snakemake.params.rmDup+ \
                  " ASSUME_SORTED=true" + \
                  " PROGRAM_RECORD_ID=null"+ \
                  " VALIDATION_STRINGENCY=LENIENT"+ \
                  " SPLIT=_" + \
                  " MM=1" + \
                  " 2>> "+snakemake.log.run+" "

        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)

        command = "samtools index -@" + str(snakemake.threads) + " " + snakemake.output.bam
        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)

else:

    shell("mv -T " + snakemake.input.bam + " " + snakemake.output.bam)
    shell("mv -T " + snakemake.input.bai + " " + snakemake.output.bai)
