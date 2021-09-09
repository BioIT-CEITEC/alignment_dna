######################################
# wrapper for rule: trim_adapters
######################################
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: trim_adapters \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

if snakemake.params.paired == "PE":
        paired_flag = " --paired"
else:
        paired_flag = ""

command = "trim_galore --fastqc " + paired_flag + " --basename "+ snakemake.input + " " + snakemake.input +" -o "+snakemake.params.outdir+ " 2>> "+log_filename
with open(log_filename, 'at') as f:
        f.write("## COMMAND: " + command + "\n")
shell(command)
