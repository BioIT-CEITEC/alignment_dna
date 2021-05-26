# ALIGNMENT RULES
#

#
# def pull_merge_after_alignment_input(wildcards):
#     if config["replicates"] != False:
#         if expand(os.path.isfile("mapped/{sample}{rep}.bam"),rep = config["replicates"]):
#             return "mapped/{sample}{rep}.bam"
#         else:
#             return expand("mapped/{sample}{rep}.{umi_usage}.bam",umi_usage = sample_tab.umi_usage)
#     else:
#         return "mapped/{sample}.{umi_usage}.bam"
#
#
# rule pull_merge_after_alignment:
#     input:  pull_merge_after_alignment_input
#     output: bam = "mapped/{sample}.bam",
#             bai = "mapped/{sample}.bam.bai",
#     log:    run = "sample_logs/{sample}/pull_merge_after_alignment.log"
#     shell: "touch {output}"

def mark_duplicates_ref(wildcards):
    input = {
        'bam': expand("mapped/{sample}.not_markDups.bam",sample=sample_tab.sample_name)[0],
        'bai': expand("mapped/{sample}.not_markDups.bam.bai",sample=sample_tab.sample_name)[0],
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
    log: run = "sample_logs/{sample}/mark_duplicates.log"
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


def alignment_DNA_input(wildcards):
    if config["preprocess"] == True:
        preprocessed = "cleaned_fastq"
    else:
        preprocessed = "raw_fastq"
    if read_pair_tags == "":
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
    log: "sample_logs/{sample}/alignment_DNA.log"
    params: lib_name=config["library_name"]
    threads: 40
    conda: "../wrappers/alignment_DNA/env.yaml"
    script: "../wrappers/alignment_DNA/script.py"
