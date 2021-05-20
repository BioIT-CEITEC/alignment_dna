######################################
# wrapper for rule: umi_concensus
######################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell


f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: umi_concensus \n##\n")
f.close()

FGBIO = "fgbio"
SAMTOOLS = "samtools"
PICARD = "picard"
BWA = "bwa"

shell.executable("/bin/bash")
bwa_ref_prefix = re.sub(".bwt$","",snakemake.input.bwa_ref)

f = open(snakemake.log.run, 'at')
umi_version = str(subprocess.Popen(FGBIO+" --version 2>&1 | grep \"[Vv]ersion\" | cut -f 2 -d \" \"",shell=True,stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f.write("## VERSION: Fgbio "+umi_version+"\n")

version = str(subprocess.Popen(BWA+" 2>&1 | grep \"[Vv]ersion\" | cut -f 2 -d \" \"", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f.write("## VERSION: BWA "+version+"\n")

version = str(subprocess.Popen(SAMTOOLS+" --version 2>&1 | grep \"samtools\" ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f.write("## VERSION: samtools "+version+"\n")

f.close()

pair_reads = int(re.sub("[^0-9]","",str(subprocess.Popen(SAMTOOLS+" view -c -f 1 "+snakemake.input.bam, shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')))
print(pair_reads)

os.makedirs(os.path.dirname(snakemake.params.umi_histogram),exist_ok=True)

if pair_reads > 1000:
    bwa_opt = " -p "
    bam_to_fastq_opt = " INTERLEAVE=true "
else:
    bwa_opt = ""
    bam_to_fastq_opt = ""

command = FGBIO+" -Xmx10g SortBam -i "+snakemake.input.bam+" -o /dev/stdout -s Queryname 2>> "+snakemake.log.run+" | " +\
          FGBIO+" -Xmx10g SetMateInformation -i /dev/stdin -o /dev/stdout 2>> "+snakemake.log.run+" | " +\
          FGBIO+" -Xmx10g AnnotateBamWithUmis -i /dev/stdin -f "+snakemake.params.umi_fastq+" -o /dev/stdout -t RX 2>> "+snakemake.log.run +" | "+ \
          SAMTOOLS+" sort -@ "+str(snakemake.threads)+" -o "+snakemake.params.umi_temp_file+" /dev/stdin 2>> "+snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = SAMTOOLS+" index -b "+snakemake.params.umi_temp_file+" 2>> "+snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mkdir -p " + os.path.dirname(snakemake.params.umi_histogram)
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = FGBIO+" -Xmx10g GroupReadsByUmi --input "+snakemake.params.umi_temp_file \
                  +" --output /dev/stdout" \
                  +" --family-size-histogram "+snakemake.params.umi_histogram \
                  +" --strategy=adjacency 2>> "+snakemake.log.run +" | "+ \
          FGBIO+" -Xmx10g CallMolecularConsensusReads --input /dev/stdin" \
                  +" --output "+snakemake.params.umi_temp_file2 \
                  +" --min-reads="+str(snakemake.params.min_umi_size) \
                  +" --min-input-base-quality="+str(snakemake.params.min_input_base_quality) \
                  +" --error-rate-pre-umi="+str(snakemake.params.error_rate_pre_umi) \
                  +" --error-rate-post-umi="+str(snakemake.params.error_rate_post_umi) \
                  +" >> "+snakemake.log.run + " 2>&1 "
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command =  PICARD +" -Xmx10g SamToFastq INPUT="+snakemake.params.umi_temp_file2 \
                  +" F=/dev/stdout QUIET=true " + bam_to_fastq_opt \
                  +" 2>> "+snakemake.log.run + " | " + \
           BWA    +" mem -t "+str(snakemake.threads)\
                  + bwa_opt\
                  +" " + bwa_ref_prefix+" /dev/stdin 2>> "+snakemake.log.run + " | " +\
           PICARD+" -Xmx10g MergeBamAlignment VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true QUIET=true" \
                  +" UNMAPPED=" + snakemake.params.umi_temp_file2 \
                  +" ALIGNED=/dev/stdin O=" + snakemake.output.bam \
                  +" R="+snakemake.input.fa_ref \
                  +" ATTRIBUTES_TO_RETAIN=X0 ATTRIBUTES_TO_RETAIN=ZS ATTRIBUTES_TO_RETAIN=ZI ATTRIBUTES_TO_RETAIN=ZM ATTRIBUTES_TO_RETAIN=ZC ATTRIBUTES_TO_RETAIN=ZN \
                     RV=cd RV=ce ORIENTATIONS=FR MAX_GAPS=-1 SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=false"\
                  +" >> " + snakemake.log.run + " 2>&1 "
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = SAMTOOLS+" index -b "+snakemake.output.bam+" 2>> "+snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)


command = "rm "+re.sub("\.bam$","",snakemake.output.bam) + ".bai"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

#
#
# shell.executable("/bin/bash")
# bwa_ref_prefix = re.sub(".bwt$","",snakemake.input.bwa_ref[0])
#
# f = open(snakemake.log.run, 'at')
# umi_version = str(subprocess.Popen("fgbio --version 2>&1 | grep \"[Vv]ersion\" | cut -f 2 -d \" \"",shell=True,stdout=subprocess.PIPE).communicate()[0], 'utf-8')
# f.write("## VERSION: Fgbio "+umi_version+"\n")
#
# version = str(subprocess.Popen("umitools 2>&1 | grep \"[Vv]ersion\" | cut -f 2 -d \" \"", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
# f.write("## VERSION: umitools "+version+"\n")
#
# version = str(subprocess.Popen("samtools --version 2>&1 | grep \"samtools\" ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
# f.write("## VERSION: samtools "+version+"\n")
#
# f.close()
#
# # pair_reads = int(re.sub("[^0-9]","",str(subprocess.Popen(SAMTOOLS+" view -c -f 1 "+snakemake.input.bam, shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')))
# # print(pair_reads)
# #
# # if pair_reads > 1000:
# #     bwa_opt = " -p "
# #     bam_to_fastq_opt = " INTERLEAVE=true "
# # else:
# #     bwa_opt = ""
# #     bam_to_fastq_opt = ""
#
# # command = FGBIO+" -Xmx10g SortBam -i "+snakemake.input.bam+" -o /dev/stdout -s Queryname 2>> "+snakemake.log.run+" | " +\
# #           FGBIO+" -Xmx10g SetMateInformation -i /dev/stdin -o /dev/stdout 2>> "+snakemake.log.run+" | " +\
# #           FGBIO+" -Xmx10g AnnotateBamWithUmis -i /dev/stdin -f "+snakemake.params.umi_fastq+" -o /dev/stdout -t RX 2>> "+snakemake.log.run +" | "+ \
# #           SAMTOOLS+" sort -@ "+str(snakemake.threads)+" -o "+snakemake.params.umi_temp_file+" /dev/stdin 2>> "+snakemake.log.run
# # f = open(snakemake.log.run, 'at')
# # f.write("## COMMAND: "+command+"\n")
# # f.close()
# # shell(command)
# #
# # command = SAMTOOLS+" index -b "+snakemake.params.umi_temp_file+" 2>> "+snakemake.log.run
# # f = open(snakemake.log.run, 'at')
# # f.write("## COMMAND: "+command+"\n")
# # f.close()
# # shell(command)
# #
# # command = "mkdir -p " + os.path.dirname(snakemake.params.umi_histogram)
# # f = open(snakemake.log.run, 'at')
# # f.write("## COMMAND: "+command+"\n")
# # f.close()
# # shell(command)
#
#
#
# # command = "umi_tools group -I "+snakemake.input.bam \
# #                   +" --output-bam --stdout " + snakemake.params.umi_temp_file \
# #                   +" --paired" \
# #                   +" --log="+snakemake.log.run +" | " \
# #                   + " fgbio -Xmx10g --tmp-dir=/mnt/ssd/ssd_1/tmp SortBam --input /dev/stdin" \
# #                   + " --sort-order=Coordinate --output /dev/stdout" \
# #                   + " fgbio -Xmx10g CallMolecularConsensusReads --input /dev/stdin" \
# #                   +" --output /dev/stdout" \
# #                   +" --min-reads="+str(snakemake.params.min_umi_size) \
# #                   +" --min-input-base-quality="+str(snakemake.params.min_input_base_quality) \
# #                   +" --error-rate-pre-umi="+str(snakemake.params.error_rate_pre_umi) \
# #                   +" --error-rate-post-umi="+str(snakemake.params.error_rate_post_umi) \
# #                   +" 2>> "+snakemake.log.run + " |"\
# #                   +" samtools sort -@ "+ str(snakemake.threads) +" -o "+snakemake.output.bam+" /dev/stdin 2>> "+snakemake.log.run
#
# os.makedirs(os.path.dirname(snakemake.params.umi_histogram),exist_ok=True)
#
# command = "fgbio -Xmx10g GroupReadsByUmi --input "+snakemake.params.umi_temp_file2 \
#                   +" --output " + snakemake.params.umi_temp_file2 + ".test.bam" \
#                   +" --family-size-histogram "+snakemake.params.umi_histogram \
#                   +" -t BX --strategy=adjacency"
#
# # command = " fgbio -Xmx50g CallMolecularConsensusReads --input " + snakemake.params.umi_temp_file2 \
# #                   +" --output " + snakemake.params.umi_temp_file2 + ".test.bam" \
# #                   +" --min-reads="+str(snakemake.params.min_umi_size) \
# #                   +" --min-input-base-quality="+str(snakemake.params.min_input_base_quality) \
# #                   +" --error-rate-pre-umi="+str(snakemake.params.error_rate_pre_umi) \
# #                   +" --error-rate-post-umi="+str(snakemake.params.error_rate_post_umi)
#
#
# # command = "fgbio -Xmx10g --tmp-dir=/mnt/ssd/ssd_1/tmp SortBam --input " + snakemake.params.umi_temp_file \
# #                   + " -s TemplateCoordinate --output " + snakemake.params.umi_temp_file2
#
#
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)
#
# command = SAMTOOLS+" index -b "+snakemake.output.bam+" 2>> "+snakemake.log.run
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)
#
# # command =  PICARD +" -Xmx10g SamToFastq INPUT="+snakemake.params.umi_temp_file2 \
# #                   +" F=/dev/stdout QUIET=true " + bam_to_fastq_opt \
# #                   +" 2>> "+snakemake.log.run + " | " + \
# #            BWA    +" mem -t "+str(snakemake.threads)\
# #                   + bwa_opt\
# #                   +" " + bwa_ref_prefix+" /dev/stdin 2>> "+snakemake.log.run + " | " +\
# #            PICARD+" -Xmx10g MergeBamAlignment VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true QUIET=true" \
# #                   +" UNMAPPED=" + snakemake.params.umi_temp_file2 \
# #                   +" ALIGNED=/dev/stdin O=" + snakemake.output.bam \
# #                   +" R="+snakemake.input.fa_ref[0] \
# #                   +" ATTRIBUTES_TO_RETAIN=X0 ATTRIBUTES_TO_RETAIN=ZS ATTRIBUTES_TO_RETAIN=ZI ATTRIBUTES_TO_RETAIN=ZM ATTRIBUTES_TO_RETAIN=ZC ATTRIBUTES_TO_RETAIN=ZN \
# #                      RV=cd RV=ce ORIENTATIONS=FR MAX_GAPS=-1 SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=false"\
# #                   +" >> " + snakemake.log.run + " 2>&1 "
# # f = open(snakemake.log.run, 'at')
# # f.write("## COMMAND: "+command+"\n")
# # f.close()
# # shell(command)
#
#
#
#
# # command = "rm "+re.sub("\.bam$","",snakemake.output.bam) + ".bai"
# # f = open(snakemake.log.run, 'at')
# # f.write("## COMMAND: "+command+"\n")
# # f.close()
# # shell(c
