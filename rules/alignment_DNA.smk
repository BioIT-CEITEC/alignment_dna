

rule alignment_DNA_multiqc:
    input:  bam = expand("mapped/{sample}.bam",sample = sample_tab.sample_name),
            idxstats = expand("qc_reports/{sample}/qc_samtools/{sample}.idxstats.tsv",sample = sample_tab.sample_name),
    output: html= "qc_reports/all_samples/alignment_DNA_multiqc/multiqc.html"
    log:    "logs/all_samples/alignment_DNA_multiqc.log"
    params: trim_adapters=config["trim_adapters"],
            mark_duplicates=config["mark_duplicates"]
    conda: "../wrappers/alignment_DNA_multiqc/env.yaml"
    script: "../wrappers/alignment_DNA_multiqc/script.py"


rule qc_samtools:
    input:  bam = "mapped/{sample}.bam"
    output: idxstats = "qc_reports/{sample}/qc_samtools/{sample}.idxstats.tsv",
            flagstats = "qc_reports/{sample}/qc_samtools/{sample}.flagstat.tsv",
    log:    "logs/{sample}/qc_samtools.log"
    threads:    1
    conda: "../wrappers/qc_samtools/env.yaml"
    script: "../wrappers/qc_samtools/script.py"


def mark_duplicates_ref(wildcards):
    input = {
        'bam': "mapped/{sample}.not_markDups.bam",
        'bai': "mapped/{sample}.not_markDups.bam.bai",
    }
    if config["umi_usage"] == "umi_concensus":
        input['ref'] = expand("{ref_dir}/index/BWA/{ref}.bwt",ref_dir=reference_directory,ref=config["reference"])[0],
        input['lib_ROI'] = expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
        input['fa'] = expand("{ref_dir}/seq/{ref}.fa",ref_dir=reference_directory,ref=config["reference"])[0],
    return input


rule mark_duplicates:
    input:  unpack(mark_duplicates_ref)
    output: bam = "mapped/{sample}.bam",
            bai = "mapped/{sample}.bam.bai"
    log:    "logs/{sample}/mark_duplicates.log"
    threads: 8
    resources: mem=10
    params: mtx= "qc_reports/{sample}/MarkDuplicates/{sample}.markDups_metrics.txt",
            mark_duplicates=config["mark_duplicates"],
            rmDup=config["remove_duplicates"],# allow possibility for rm duplicates true
            UMI=config["UMI"],
            umi_usage=config["umi_usage"],
            keep_not_markDups_bam=config["keep_not_markDups_bam"],
            umi_consensus_min_support=config["umi_consensus_min_support"],
            report_path="qc_reports/{sample}/MarkDuplicates/"
    conda: "../wrappers/mark_duplicates/env.yaml"
    script: "../wrappers/mark_duplicates/script.py"


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


def trim_adapters_input(wildcards):
    if config["trim_adapters"]:
        if read_pair_tags == [""]:
            return "raw_fastq/{sample}.fastq.gz"
        else:
            return ["raw_fastq/{sample}_R1.fastq.gz","raw_fastq/{sample}_R2.fastq.gz"]


rule trim_adapters:
    input:  trim_adapters_input,
    output: fastq = expand("cleaned_fastq/{{sample}}{read_pair_tag}.fastq.gz",read_pair_tag = read_pair_tags),
            trim_stats = expand("qc_reports/{{sample}}/trim_galore/trim_stats{read_pair_tag}.log",read_pair_tag=read_pair_tags)
    log:    "logs/{sample}/trim_adapters.log"
    params: paired = paired,
            outdir = "cleaned_fastq",
    threads: 8
    conda: "../wrappers/trim_adapters/env.yaml"
    script: "../wrappers/trim_adapters/script.py"
