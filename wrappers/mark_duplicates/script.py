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

shell.executable("/bin/bash")

if snakemake.params.mark_duplicates == True:
    os.makedirs(os.path.dirname(snakemake.params.mtx),exist_ok=True)

    if snakemake.params.UMI == "no_umi" or snakemake.params.umi_usage == "no":

        f = open(log_filename, 'at')
        version = str(subprocess.Popen("picard MarkDuplicates --version 2>&1",shell=True,stdout=subprocess.PIPE).communicate()[0], 'utf-8')
        f.write("## VERSION: Picard "+version+"\n")
        f.close()

        command = "picard MarkDuplicates INPUT="+snakemake.input.bam[0]+" OUTPUT="+snakemake.output.bam[0]+" METRICS_FILE="+snakemake.params.mtx+" REMOVE_DUPLICATES="+str(snakemake.params.rmDup)+ " \
            ASSUME_SORTED=true PROGRAM_RECORD_ID=null VALIDATION_STRINGENCY=LENIENT -Xmx" +str(snakemake.resources.mem)+"g 2>> "+log_filename+" "
        f = open(log_filename, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)

        command = "samtools index "+snakemake.output.bam[0]
        f = open(log_filename, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)

    else:

        if snakemake.params.umi_usage == "mark_duplicates":
            version = str(subprocess.Popen("je markdupes --version 2>&1 | tail -n 1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
            f = open(log_filename, 'at')
            f.write("## VERSION: je-suite "+version+"\n")
            f.close()
            java_opts = "export _JAVA_OPTIONS='-Xmx"+str(snakemake.resources.mem)+"g'"

            command = java_opts + "&& je markdupes " + \
                      "INPUT=" + snakemake.input.bam[0] +\
                      " OUTPUT="+snakemake.output.bam[0]+\
                      " METRICS_FILE="+snakemake.params.mtx+\
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

            command = "samtools index -@" + str(snakemake.threads) + " " + snakemake.output.bam[0]
            f = open(log_filename, 'at')
            f.write("## COMMAND: "+command+"\n")
            f.close()
            shell(command)

        else:
            version = str(subprocess.Popen("gencore --version 2>&1", shell=True,stdout=subprocess.PIPE).communicate()[0], 'utf-8')
            f = open(log_filename, 'at')
            f.write("## VERSION: gencore "+version+"\n")
            f.close()

            command ="gencore -i "+snakemake.input.bam[0]+" -o "+snakemake.output.bam[0]+" --ref "+str(snakemake.input.fa)+" --bed "+str(snakemake.input.lib_ROI)+" --supporting_reads "+str(snakemake.params.umi_consensus_min_support)+\
                     " --html "+snakemake.params.report_path+"umi_concensus.html"+\
                     " --json "+snakemake.params.report_path+"umi_concensus.json"+\
                     " 2>> "+log_filename+" "

            f = open(log_filename, 'at')
            f.write("## COMMAND: "+command+"\n")
            f.close()
            shell(command)


            command = "samtools index "+snakemake.output.bam[0]
            f = open(log_filename, 'at')
            f.write("## COMMAND: "+command+"\n")
            f.close()
            shell(command)

    if snakemake.params.keep_not_markDups_bam == False:
        command = "rm " + snakemake.input.bam[0]
        f = open(log_filename, 'at')
        f.write("## COMMAND: " + command + "\n")
        f.close()
        shell(command)


        command = "rm " + snakemake.input.bam[0] + ".bai"
        f = open(log_filename, 'at')
        f.write("## COMMAND: " + command + "\n")
        f.close()
        shell(command)

else:

    shell("mv -T " + snakemake.input.bam[0] + " " + snakemake.output.bam[0])
    shell("mv -T " + snakemake.input.bai[0] + " " + snakemake.output.bai[0])
