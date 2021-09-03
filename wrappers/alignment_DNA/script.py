######################################
# wrapper for rule: alignment_DNA
######################################
import re
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)
print(log_filename)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: alignment_DNA \n##\n")
f.close()

version = str(subprocess.Popen("bwa 2>&1 | grep \"[Vv]ersion\" | cut -f 2 -d \" \"", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## VERSION: BWA "+version+"\n")
f.close()

version = str(subprocess.Popen("samtools --version 2>&1 | grep \"samtools\" ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## VERSION: samtools "+version+"\n")
f.close()

bwa_ref_prefix = snakemake.input.ref[0].split('.')[0]
#bwa_ref_prefix = re.sub(".bwt$","",snakemake.input.ref[0])

command = "bwa mem -t " + str(snakemake.threads) + " -R '@RG\\tID:" + snakemake.params.lib_name + "_" + snakemake.wildcards.sample + "\\tSM:" + snakemake.wildcards.sample + "\\tPL:illumina' -v 1 "  + bwa_ref_prefix + " " + snakemake.input.r1[0] + " " + snakemake.input.r2[0] +  " 2>> " + log_filename + " | " + "samtools sort -@ " + str(snakemake.threads) + " -o " + snakemake.output.bam[0] + " /dev/stdin 2>> " + log_filename

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "samtools index -b "+snakemake.output.bam[0]+" 2>> "+log_filename
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
