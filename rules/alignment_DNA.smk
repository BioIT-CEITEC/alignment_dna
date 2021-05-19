# ALIGNMENT RULES
#

"""

rule mark_duplicates:
    input:
        bam = "mapped/{sample}.not_markDups.bam",
        bai = "mapped/{sample}.not_markDups.bam.bai"
    output:
        bam = "mapped/marked_duplicates/{sample}.markDups.bam",
        bai = "mapped/marked_duplicates/{sample}.markDups.bam.bai"
    log: run = "sample_logs/{sample}/mark_duplicates.log"
    threads: 8
    resources: mem = 10
    params:
        mtx = "map_qc/{sample}/MarkDuplicates/{sample}.markDups_metrics.txt",
        primer_based = config["primer_based"],
        rmDup = "false", # allow possibility for rm duplicates true
        umi = config["UMI"],
        umi_fastq = "raw_fastq/{sample}.UMI.fastq",
        umi_temp_file = "mapped/{sample}.not_markDups.UMIannot.bam"
    conda:  "../wrappers/mark_duplicates/env.yaml"
    script: "../wrappers/mark_duplicates/script.py"




rule umi_concensus:
    input:
        bam = "mapped/{sample}.not_markDups.bam",
        bai = "mapped/{sample}.not_markDups.bam.bai",
        bwa_ref = expand("{ref_dir}/index/BWA/{ref}.bwt",ref_dir=reference_directory,ref = config["reference"])[0],
        fa_ref = expand("{ref_dir}/seq/{ref}.fa",ref_dir=reference_directory,ref = config["reference"])[0],
    output:
        bam = "mapped/{sample}.umi_concensus.bam",
        bai = "mapped/{sample}.umi_concensus.bam.bai"
    log: run = "sample_logs/{sample}/umi_concensus.log"
    threads: 20
    resources: mem = 10
    params: umi_fastq = "raw_fastq/{sample}.UMI.fastq",
            umi_temp_file = "mapped/{sample}.not_markDups.UMIannot.bam",
            umi_histogram = "map_qc/{sample}.umi_concensus.histogram",
            umi_temp_file2 = "mapped/{sample}.not_markDups.UMIcons.bam",
            min_umi_size = 1,
            min_input_base_quality = 20,
            error_rate_pre_umi = 50,
            error_rate_post_umi = 30
    conda:  "../wrappers/umi_concensus/env.yaml"
    script: "../wrappers/umi_concensus/script.py"




"""

"""
"""
def merge_bams_input(wildcards):
    if config["replicates"] == True:
        return expand("mapped/reseq/{sample}.{reseq}.bam",sample = sample_tab.sample_name,reseq = config["replicates"])
    else:
        return expand("mapped/{sample}.not_markDups.bam",sample = sample_tab.sample_name)


rule merge_bams:
    input: merge_bams_input
    output:
        bam = "mapped/merged/{sample}.not_markDups.bam",
        bai = "mapped/merged/{sample}.not_markDups.bam.bai",
    log: run = "sample_logs/{sample}/merge_bams.log"
    threads: 10
    params: sample_name =  sample_tab.sample_name
    conda:  "../wrappers/merge_bams/env.yaml"
    script: "../wrappers/merge_bams/script.py"


def alignment_DNA_input(wildcards):
    if config["preprocess"] == True:
        preprocessed = "cleaned_fastq"
    else:
        preprocessed = "raw_fastq"

    #if config["lib_reverse_read_length"] == 0:
    #    return expand(os.path.join(preprocessed,"{sample}_SE.fastq.gz"),sample=sample_tab.sample_name)
    #else:
    #    return [os.path.join(preprocessed,"{sample}_R1.fastq.gz"),os.path.join(preprocessed,"{sample}_R2.fastq.gz")]

    if read_pair_tags == "":
        return expand(os.path.join(preprocessed,"{sample}.fastq.gz"),sample=sample_tab.sample_name)
    else:
        return [os.path.join(preprocessed,"{sample}_R1.fastq.gz"),os.path.join(preprocessed,"{sample}_R2.fastq.gz")]

rule alignment_DNA:
    input:
        fastqs = alignment_DNA_input,
        ref = expand("{ref_dir}/index/BWA/{ref}.bwt",ref_dir=reference_directory,ref = config["reference"])[0],
    output:
        bam = "mapped/{sample}.not_markDups.bam",
        bai = "mapped/{sample}.not_markDups.bam.bai"
    log: run = "sample_logs/{sample}/alignment_DNA.log"
    params: lib_name = config["library_name"]
    threads: 40
    conda:  "../wrappers/alignment_DNA/env.yaml"
    script: "../wrappers/alignment_DNA/script.py"

