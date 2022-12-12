######################################
# wrapper for rule: umi_concensus
######################################
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: umi_concensus \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()


command = "export TMPDIR=" + snakemake.params.tmpd + " TMP=" + snakemake.params.tmpd + " && gencore" + \
          " -i " + snakemake.input.bam + \
          " -o " + snakemake.output.bam + \
          " --ref " + str(snakemake.input.fa) + \
          " --bed " + str(snakemake.input.lib_ROI) + \
          " --supporting_reads " + str(snakemake.params.umi_consensus_min_support) + \
          " --html " + snakemake.output.html + \
          " --json " + snakemake.output.json + \
          " 2>> " + log_filename + " "


f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

if not snakemake.params.keep_not_markDups_bam:
    command = "rm " + snakemake.input.bam
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    command = "rm " + snakemake.input.bai
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)
