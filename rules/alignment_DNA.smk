
# def trim_adapters_input(wildcards):
#     if config["trim_adapters"]:
#         if read_pair_tags == [""]:
#             return "raw_fastq/{sample}.fastq.gz"
#         else:
#             return ["raw_fastq/{sample}_R1.fastq.gz","raw_fastq/{sample}_R2.fastq.gz"]
#
# rule trim_adapters:
#     input:  trim_adapters_input,
#     output: fastq = expand("cleaned_fastq/{{sample}}{read_pair_tag}.fastq.gz",read_pair_tag = read_pair_tags),
#             trim_stats = expand("qc_reports/{{sample}}/trim_galore/trim_stats{read_pair_tag}.log",read_pair_tag=read_pair_tags)
#     log:    "logs/{sample}/trim_adapters.log"
#     params: paired = config["is_paired"],
#             outdir = "cleaned_fastq",
#     threads: 8
#     conda: "../wrappers/trim_adapters/env.yaml"
#     script: "../wrappers/trim_adapters/script.py"


def alignment_DNA_input(wildcards):
    if config["trim_adapters"]:
        preprocessed = "cleaned_fastq"
    else:
        preprocessed = "raw_fastq"
    if read_pair_tags == [""]:
        return [os.path.join(preprocessed,"{sample}.fastq.gz")]
    else:
        return [os.path.join(preprocessed,"{sample}_R1.fastq.gz"),os.path.join(preprocessed,"{sample}_R2.fastq.gz")]


rule alignment_DNA:
    input:  fastqs = alignment_DNA_input,
            ref = expand("{ref_dir}/index/BWA/{ref}.bwt",ref_dir=reference_directory,ref = config["reference"])[0],
    output: bam = "mapped/{sample}.not_markDups.bam",
            bai = "mapped/{sample}.not_markDups.bam.bai"
    log:    "logs/{sample}/alignment_DNA.log"
    params: entity_name=config["entity_name"]
    threads: 40
    conda: "../wrappers/alignment_DNA/env.yaml"
    script: "../wrappers/alignment_DNA/script.py"


rule mark_duplicates:
    input:  bam = "mapped/{sample}.not_markDups.bam",
            bai = "mapped/{sample}.not_markDups.bam.bai",
    output: bam = "mapped/{sample}.markDups.bam",
    log:    "logs/{sample}/mark_duplicates.log"
    resources: mem=10
    params: mark_duplicates=config["mark_duplicates"],
            rmDup=config["remove_duplicates"],
            UMI=config["UMI"],
            umi_usage=config["umi_usage"],
            keep_not_markDups_bam=config["keep_not_markDups_bam"],
            tmpd = GLOBAL_TMPD_PATH,
            mtx = "qc_reports/{sample}/MarkDuplicates/{sample}.markDups_metrics.txt"
    conda: "../wrappers/mark_duplicates/env.yaml"
    script: "../wrappers/mark_duplicates/script.py"


rule umi_concensus:
    input:  bam = "mapped/{sample}.not_markDups.bam",
            bai = "mapped/{sample}.not_markDups.bam.bai",
            ref = expand("{ref_dir}/index/BWA/{ref}.bwt",ref_dir=reference_directory,ref=config["reference"])[0],
            lib_ROI = expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
            fa = expand("{ref_dir}/seq/{ref}.fa",ref_dir=reference_directory,ref=config["reference"])[0],
    output: bam = "mapped/{sample}.concensus.bam",
            html = "qc_reports/{sample}/umi_concensus/umi_concensus.html",
            json = "qc_reports/{sample}/umi_concensus/umi_concensus.json",
    log: "logs/{sample}/umi_concensus.log"
    params: umi_consensus_min_support = config["umi_consensus_min_support"],
            keep_not_markDups_bam=config["keep_no_umi_consensus_bam"],
            tmpd = GLOBAL_TMPD_PATH
    conda: "../wrappers/umi_concensus/env.yaml"
    script: "../wrappers/umi_concensus/script.py"


def index_and_stats_input(wildcards):
    if not config["umi_usage"] == "umi_concensus":
        return "mapped/{sample}.markDups.bam"
    else:
        return "mapped/{sample}.concensus.bam"

rule index_and_stats:
    input:  index_and_stats_input,
    output: bam = "mapped/{sample}.bam",
            bai = "mapped/{sample}.bam.bai",
            idxstats = "qc_reports/{sample}/index_and_stats/{sample}.idxstats.tsv",
            flagstats = "qc_reports/{sample}/index_and_stats/{sample}.flagstat.tsv",
    log:    "logs/{sample}/index_and_stats.log"
    threads:    8
    conda: "../wrappers/index_and_stats/env.yaml"
    script: "../wrappers/index_and_stats/script.py"


rule alignment_DNA_multiqc:
    input:  bam = expand("mapped/{sample}.bam",sample = sample_tab.sample_name),
            idxstats = expand("qc_reports/{sample}/index_and_stats/{sample}.idxstats.tsv",sample = sample_tab.sample_name),
    output: html= "qc_reports/all_samples/alignment_DNA_multiqc/multiqc.html"
    log:    "logs/all_samples/alignment_DNA_multiqc.log"
    params: mark_duplicates=config["mark_duplicates"],
            umi_usage = config["umi_usage"]
    conda: "../wrappers/alignment_DNA_multiqc/env.yaml"
    script: "../wrappers/alignment_DNA_multiqc/script.py"

