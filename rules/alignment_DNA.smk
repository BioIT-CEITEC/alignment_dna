

def alignment_DNA_multiqc_input(wildcards):
    input = {
        "bam": BR.remote(expand("mapped/{sample}.bam",sample=sample_tab.sample_name)),
        "samtools": BR.remote(expand("qc_reports/{sample}/index_and_stats/{sample}.{stat}.tsv",sample=sample_tab.sample_name,stat=["flagstat", "idxstats"])),
    }
    if config["trim_adapters"]:
        input["trim_galore"] = BR.remote(expand("qc_reports/{sample}/trim_galore/trim_stats{read_pair_tag}.log",sample=sample_tab.sample_name,read_pair_tag=read_pair_tags))
    if config["mark_duplicates"]:
        input["mark_duplicates"] = BR.remote(expand("qc_reports/{sample}/MarkDuplicates/{sample}.markDups_metrics.txt",sample=sample_tab.sample_name))
    return input

rule alignment_DNA_multiqc:
    input: unpack(alignment_DNA_multiqc_input)
    output: html=BR.remote("qc_reports/all_samples/alignment_DNA_multiqc/multiqc.html")
    log: BR.remote("logs/all_samples/alignment_DNA_multiqc.log")
    params: trim_adapters=config["trim_adapters"],
            mark_duplicates=config["mark_duplicates"]
    conda: "../wrappers/alignment_DNA_multiqc/env.yaml"
    script: "../wrappers/alignment_DNA_multiqc/script.py"


def index_and_stats_input(wildcards):
    if not config["umi_usage"] == "umi_concensus":
        return BR.remote("mapped/{sample}.markDups.bam")
    else:
        return BR.remote("mapped/{sample}.concensus.bam")

rule index_and_stats:
    input:  index_and_stats_input
    output: bam = BR.remote("mapped/{sample}.bam"),
            bai = BR.remote("mapped/{sample}.bam.bai"),
            idxstats = BR.remote("qc_reports/{sample}/index_and_stats/{sample}.idxstats.tsv"),
            flagstats = BR.remote("qc_reports/{sample}/index_and_stats/{sample}.flagstat.tsv"),
    log:    BR.remote("logs/{sample}/index_and_stats.log")
    threads:    8
    conda: "../wrappers/index_and_stats/env.yaml"
    script: "../wrappers/index_and_stats/script.py"


rule umi_concensus:
    input:  bam = BR.remote("mapped/{sample}.not_markDups.bam"),
            ref = BR.remote(expand("{ref_dir}/tool_data/BWA/{ref}.{end}",ref_dir=reference_directory,ref=config["reference"],end=["amb", "ann", "bwt", "pac", "sa"])),
            lib_ROI = BR.remote(expand("{ref_dir}/DNA_panel/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])),
            fa = BR.remote(expand("{ref_dir}/seq/{ref}.fa",ref_dir=reference_directory,ref=config["reference"])),
    output: bam = BR.remote("mapped/{sample}.concensus.bam"),
            html = BR.remote("qc_reports/{sample}/umi_concensus/umi_concensus.html"),
            json = BR.remote("qc_reports/{sample}/umi_concensus/umi_concensus.json"),
    log: BR.remote("logs/{sample}/umi_concensus.log")
    params: umi_consensus_min_support = config["umi_consensus_min_support"],
    conda: "../wrappers/umi_concensus/env.yaml"
    script: "../wrappers/umi_concensus/script.py"


rule mark_duplicates:
    input: bam = BR.remote("mapped/{sample}.not_markDups.bam"),
    output: bam = BR.remote("mapped/{sample}.markDups.bam"),
            mtx= BR.remote("qc_reports/{sample}/MarkDuplicates/{sample}.markDups_metrics.txt")
    log:    BR.remote("logs/{sample}/mark_duplicates.log")
    threads: 8
    resources: mem=10
    params: mark_duplicates=config["mark_duplicates"],
            rmDup=config["remove_duplicates"],
            UMI=config["UMI"],
            umi_usage=config["umi_usage"],
            keep_not_markDups_bam=config["keep_not_markDups_bam"],
    conda: "../wrappers/mark_duplicates/env.yaml"
    script: "../wrappers/mark_duplicates/script.py"


def alignment_DNA_input(wildcards):
    if config["trim_adapters"]:
        preprocessed = "cleaned_fastq"
    else:
        preprocessed = "raw_fastq"
    if read_pair_tags == [""]:
        return [BR.remote(os.path.join(preprocessed,"{sample}.fastq.gz"))]
    else:
        return [BR.remote(os.path.join(preprocessed,"{sample}_R1.fastq.gz")),
                BR.remote(os.path.join(preprocessed,"{sample}_R2.fastq.gz"))]


rule alignment_DNA:
    input:  fastqs = alignment_DNA_input,
            ref= BR.remote(expand("{ref_dir}/tool_data/BWA/{ref}.{end}",ref_dir=reference_directory,ref=config["reference"],org=config["organism"],end=["amb", "ann", "bwt", "pac", "sa"])),
    output: bam = BR.remote("mapped/{sample}.not_markDups.bam"),
            bai = BR.remote("mapped/{sample}.not_markDups.bam.bai")
    log:    BR.remote("logs/{sample}/alignment_DNA.log")
    params: entity_name=config["entity_name"]
    threads: 40
    conda: "../wrappers/alignment_DNA/env.yaml"
    script: "../wrappers/alignment_DNA/script.py"


def trim_adapters_input(wildcards):
    if read_pair_tags == [""]:
        return [BR.remote("raw_fastq/{sample}.fastq.gz")]
    else:
        return [BR.remote("raw_fastq/{sample}_R1.fastq.gz"),
                BR.remote("raw_fastq/{sample}_R2.fastq.gz")]

rule trim_adapters:
    input:  trim_adapters_input,
    output: fastq = BR.remote(expand("cleaned_fastq/{{sample}}{read_pair_tag}.fastq.gz",read_pair_tag=read_pair_tags)),
            trim_stats = BR.remote(expand("qc_reports/{{sample}}/trim_galore/trim_stats{read_pair_tag}.log",read_pair_tag=read_pair_tags)),
    log:    BR.remote("logs/{sample}/trim_adapters.log")
    params: paired = config["is_paired"],
    threads: 8
    conda: "../wrappers/trim_adapters/env.yaml"
    script: "../wrappers/trim_adapters/script.py"


