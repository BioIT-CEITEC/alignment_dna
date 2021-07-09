# ####################################
# # REFERENCE INFO PREPARATION
# #
#
# rule ref_info_copy:
#     input:  expand("{ref_dir}/info.txt", zip , ref_dir=reference_directory)[0],
#     output: "genomic_reference_info.txt",
#     run:
#         shell(" cp {input} {output} ")



def alignment_DNA_input(wildcards):
    if config["trim_adapters"] == True or config["trim_quality"] == True:
        preprocessed = "cleaned_fastq"
    else:
        preprocessed = "raw_fastq"
    if read_pair_tags == [""]:
        return os.path.join(preprocessed,"{sample}.fastq.gz")
    else:
        return [os.path.join(preprocessed,"{sample}_R1.fastq.gz"),os.path.join(preprocessed,"{sample}_R2.fastq.gz")]

rule alignment_DNA:
    input:
        fastqs = alignment_DNA_input,
        ref = expand("{ref_dir}/index/BWA/{ref}.bwt",ref_dir=reference_directory,ref = config["reference"])[0],
    output:
        bam = "mapped/{sample}.not_markDups.bam",
        bai = "mapped/{sample}.not_markDups.bam.bai"
    log: "logs/{sample}/alignment_DNA.log"
    params: lib_name=config["library_name"]
    threads: 40
    conda: "../wrappers/alignment_DNA/env.yaml"
    script: "../wrappers/alignment_DNA/script.py"


def mark_duplicates_ref(wildcards):
    input = {
        'bam': "mapped/{sample}.not_markDups.bam",
        'bai': "mapped/{sample}.not_markDups.bam.bai",
    }
    if config["umi_usage"] == "umi_concensus":
        input['ref'] = expand("{ref_dir}/index/BWA/{ref}.bwt",ref_dir=reference_directory,ref=config["reference"])[0],
        input['lib_ROI'] = expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
        input['fa'] = expand("{ref_dir}/seq/{ref}.fa",ref_dir=reference_directory, ref = config["reference"])[0],
    return input

rule mark_duplicates:
    input: unpack(mark_duplicates_ref)
    output:
        bam="mapped/{sample}.bam",
        bai="mapped/{sample}.bam.bai"
    log: "logs/{sample}/mark_duplicates.log"
    threads: 8
    resources: mem = 10
    params:
        mtx = "map_qc/{sample}/MarkDuplicates/{sample}.markDups_metrics.txt",
        mark_duplicates = config["mark_duplicates"],
        rmDup = config["remove_duplicates"], # allow possibility for rm duplicates true
        UMI = config["UMI"],
        umi_usage = config["umi_usage"],
        keep_not_markDups_bam= config["keep_not_markDups_bam"],
        umi_consensus_min_support=config["umi_consensus_min_support"],
        report_path = "map_qc/{sample}/MarkDuplicates/"
    conda:  "../wrappers/mark_duplicates/env.yaml"
    script: "../wrappers/mark_duplicates/script.py"