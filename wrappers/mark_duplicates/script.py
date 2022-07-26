######################################
# wrapper for rule: mark_duplicates
######################################
import os
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: mark_duplicates \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

if snakemake.params.mark_duplicates == True:
    os.makedirs(os.path.dirname(snakemake.output.mtx),exist_ok=True)

    if snakemake.params.UMI == "no_umi" or snakemake.params.umi_usage == "no":

        command = "picard MarkDuplicates INPUT="+snakemake.input.bam+" OUTPUT="+snakemake.output.bam+" METRICS_FILE="+snakemake.output.mtx+" REMOVE_DUPLICATES="+str(snakemake.params.rmDup)+ " \
            ASSUME_SORTED=true PROGRAM_RECORD_ID=null VALIDATION_STRINGENCY=LENIENT -Xmx" +str(snakemake.resources.mem)+"g 2>> "+log_filename+" "
        f = open(log_filename, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)

        # command = "samtools index "+snakemake.output.bam
        # f = open(log_filename, 'at')
        # f.write("## COMMAND: "+command+"\n")
        # f.close()
        # shell(command)

    else:

        java_opts = "export _JAVA_OPTIONS='-Xmx"+str(snakemake.resources.mem)+"g'"

        command = java_opts + "&& je markdupes " + \
                  "INPUT=" + snakemake.input.bam +\
                  " OUTPUT="+snakemake.output.bam+\
                  " METRICS_FILE="+snakemake.output.mtx+\
                  " REMOVE_DUPLICATES="+str(snakemake.params.rmDup)+ \
                  " ASSUME_SORTED=true" + \
                  " PROGRAM_RECORD_ID=null"+ \
                  " VALIDATION_STRINGENCY=LENIENT"+ \
                  " SPLIT=_" + \
                  " MM=1" + \
                  " 2>> "+log_filename+" "

        f = open(log_filename, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)

        # command = "samtools index -@" + str(snakemake.threads) + " " + snakemake.output.bam
        # f = open(log_filename, 'at')
        # f.write("## COMMAND: "+command+"\n")
        # f.close()
        # shell(command)



    # if snakemake.params.keep_not_markDups_bam == False:
    #
    #     command = "rm " + snakemake.input.bam
    #     f = open(log_filename, 'at')
    #     f.write("## COMMAND: " + command + "\n")
    #     f.close()
    #     shell(command)


        # command = "rm " + snakemake.input.bam + ".bai"
        # f = open(log_filename, 'at')
        # f.write("## COMMAND: " + command + "\n")
        # f.close()
        # shell(command)

else:

    shell("mv -T " + snakemake.input.bam + " " + snakemake.output.bam)
    #shell("mv -T " + snakemake.input.bai + " " + snakemake.output.bai)
