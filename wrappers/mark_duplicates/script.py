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

if snakemake.params.mark_duplicates == True:
    os.makedirs(os.path.dirname(snakemake.params.mtx),exist_ok=True)

    if snakemake.params.UMI == "no_umi" or snakemake.params.umi_usage == "no":

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

        # command = "samtools index "+snakemake.output.bam
        # f = open(snakemake.log.run, 'at')
        # f.write("## COMMAND: "+command+"\n")
        # f.close()
        # shell(command)



    else:

        if snakemake.params.umi_usage == "mark_duplicates":
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

            # command = "samtools index -@" + str(snakemake.threads) + " " + snakemake.output.bam
            # f = open(snakemake.log.run, 'at')
            # f.write("## COMMAND: "+command+"\n")
            # f.close()
            # shell(command)



        else:

            f = open(snakemake.log.run, 'a+')
            f.write("\n##\n## RULE: umi_concensus \n##\n")
            f.close()

            FGBIO = "fgbio"
            SAMTOOLS = "samtools"
            PICARD = "picard"
            BWA = "bwa"

            bwa_ref_prefix = re.sub(".bwt$", "", snakemake.input.bwa_ref)

            f = open(snakemake.log.run, 'at')
            umi_version = str(
                subprocess.Popen(FGBIO + " --version 2>&1 | grep \"[Vv]ersion\" | cut -f 2 -d \" \"", shell=True,
                                 stdout=subprocess.PIPE).communicate()[0], 'utf-8')
            f.write("## VERSION: Fgbio " + umi_version + "\n")

            version = str(subprocess.Popen(BWA + " 2>&1 | grep \"[Vv]ersion\" | cut -f 2 -d \" \"", shell=True,
                                           stdout=subprocess.PIPE).communicate()[0], 'utf-8')
            f.write("## VERSION: BWA " + version + "\n")

            version = str(subprocess.Popen(SAMTOOLS + " --version 2>&1 | grep \"samtools\" ", shell=True,
                                           stdout=subprocess.PIPE).communicate()[0], 'utf-8')
            f.write("## VERSION: samtools " + version + "\n")

            f.close()

            pair_reads = int(re.sub("[^0-9]", "", str(
                subprocess.Popen(SAMTOOLS + " view -c -f 1 " + snakemake.input.bam, shell=True,
                                 stdout=subprocess.PIPE).communicate()[0], 'utf-8')))
            print(pair_reads)

            os.makedirs(os.path.dirname(snakemake.params.umi_histogram), exist_ok=True)

            if pair_reads > 1000:
                bwa_opt = " -p "
                bam_to_fastq_opt = " INTERLEAVE=true "
            else:
                bwa_opt = ""
                bam_to_fastq_opt = ""

            command = FGBIO + " -Xmx10g SortBam -i " + snakemake.input.bam + " -o /dev/stdout -s Queryname 2>> " + snakemake.log.run + " | " + \
                      FGBIO + " -Xmx10g SetMateInformation -i /dev/stdin -o /dev/stdout 2>> " + snakemake.log.run + " | " + \
                      FGBIO + " -Xmx10g AnnotateBamWithUmis -i /dev/stdin -f " + snakemake.params.umi_fastq + " -o /dev/stdout -t RX 2>> " + snakemake.log.run + " | " + \
                      SAMTOOLS + " sort -@ " + str(
                snakemake.threads) + " -o " + snakemake.params.umi_temp_file + " /dev/stdin 2>> " + snakemake.log.run
            f = open(snakemake.log.run, 'at')
            f.write("## COMMAND: " + command + "\n")
            f.close()
            shell(command)

            command = SAMTOOLS + " index -b " + snakemake.params.umi_temp_file + " 2>> " + snakemake.log.run
            f = open(snakemake.log.run, 'at')
            f.write("## COMMAND: " + command + "\n")
            f.close()
            shell(command)

            command = "mkdir -p " + os.path.dirname(snakemake.params.umi_histogram)
            f = open(snakemake.log.run, 'at')
            f.write("## COMMAND: " + command + "\n")
            f.close()
            shell(command)

            command = FGBIO + " -Xmx10g GroupReadsByUmi --input " + snakemake.params.umi_temp_file \
                      + " --output /dev/stdout" \
                      + " --family-size-histogram " + snakemake.params.umi_histogram \
                      + " --strategy=adjacency 2>> " + snakemake.log.run + " | " + \
                      FGBIO + " -Xmx10g CallMolecularConsensusReads --input /dev/stdin" \
                      + " --output " + snakemake.params.umi_temp_file2 \
                      + " --min-reads=" + str(snakemake.params.min_umi_size) \
                      + " --min-input-base-quality=" + str(snakemake.params.min_input_base_quality) \
                      + " --error-rate-pre-umi=" + str(snakemake.params.error_rate_pre_umi) \
                      + " --error-rate-post-umi=" + str(snakemake.params.error_rate_post_umi) \
                      + " >> " + snakemake.log.run + " 2>&1 "
            f = open(snakemake.log.run, 'at')
            f.write("## COMMAND: " + command + "\n")
            f.close()
            shell(command)

            command = PICARD + " -Xmx10g SamToFastq INPUT=" + snakemake.params.umi_temp_file2 \
                      + " F=/dev/stdout QUIET=true " + bam_to_fastq_opt \
                      + " 2>> " + snakemake.log.run + " | " + \
                      BWA + " mem -t " + str(snakemake.threads) \
                      + bwa_opt \
                      + " " + bwa_ref_prefix + " /dev/stdin 2>> " + snakemake.log.run + " | " + \
                      PICARD + " -Xmx10g MergeBamAlignment VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true QUIET=true" \
                      + " UNMAPPED=" + snakemake.params.umi_temp_file2 \
                      + " ALIGNED=/dev/stdin O=" + snakemake.output.bam \
                      + " R=" + snakemake.input.fa_ref \
                      + " ATTRIBUTES_TO_RETAIN=X0 ATTRIBUTES_TO_RETAIN=ZS ATTRIBUTES_TO_RETAIN=ZI ATTRIBUTES_TO_RETAIN=ZM ATTRIBUTES_TO_RETAIN=ZC ATTRIBUTES_TO_RETAIN=ZN \
                                 RV=cd RV=ce ORIENTATIONS=FR MAX_GAPS=-1 SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=false" \
                      + " >> " + snakemake.log.run + " 2>&1 "
            f = open(snakemake.log.run, 'at')
            f.write("## COMMAND: " + command + "\n")
            f.close()
            shell(command)

            command = SAMTOOLS + " index -b " + snakemake.output.bam + " 2>> " + snakemake.log.run
            f = open(snakemake.log.run, 'at')
            f.write("## COMMAND: " + command + "\n")
            f.close()
            shell(command)

            command = "rm " + re.sub("\.bam$", "", snakemake.output.bam) + ".bai"
            f = open(snakemake.log.run, 'at')
            f.write("## COMMAND: " + command + "\n")
            f.close()
            shell(command)


    if keep_not_markDups_bam == False:
        command = "rm " + snakemake.input.bam
        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: " + command + "\n")
        f.close()
        shell(command)

        command = "rm " + snakemake.input.bam + ".bai"
        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: " + command + "\n")
        f.close()
        shell(command)

else:

    shell("mv -T " + snakemake.input.bam + " " + snakemake.output.bam)
    shell("mv -T " + snakemake.input.bai + " " + snakemake.output.bai)
