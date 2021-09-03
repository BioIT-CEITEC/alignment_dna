# ####################################
# # REFERENCE INFO PREPARATION
# #
#
# rule ref_info_copy:
#     input:  expand("{ref_dir}/info.txt", zip , ref_dir=reference_directory)[0],
#     output: "genomic_reference_info.txt",
#     run:
#         shell(" cp {input} {output} ")


# multiqc_search_paths = "./qc_reports/*/MarkDuplicates/*" + " ./mapped/*" + " ./qc_reports/*/qc_samtools/*"

rule alignment_DNA_multiqc:
    input:  bam = S3.remote(expand("{bucket}/mapped/{sample}.bam",sample = sample_tab.sample_name,bucket=S3_BUCKET)),
            idxstats = S3.remote(expand("{bucket}/qc_reports/{sample}/qc_samtools/{sample}.idxstats.tsv",sample = sample_tab.sample_name,bucket=S3_BUCKET)),
            flagstats = S3.remote(expand("{bucket}/qc_reports/{sample}/qc_samtools/{sample}.flagstat.tsv",sample = sample_tab.sample_name,bucket=S3_BUCKET)),
    output: html=S3.remote(expand("{bucket}/qc_reports/all_samples/alignment_DNA_multiqc/multiqc.html",bucket=S3_BUCKET))
    log: S3.remote(expand("{bucket}/logs/all_samples/alignment_DNA_multiqc.log",bucket=S3_BUCKET))
    conda: "../wrappers/alignment_DNA_multiqc/env.yaml"
    script: "../wrappers/alignment_DNA_multiqc/script.py"


rule qc_samtools:
    input:  bam = S3.remote(expand("{bucket}/mapped/{{sample}}.bam",bucket=S3_BUCKET))
    output: idxstats = S3.remote(expand("{bucket}/qc_reports/{{sample}}/qc_samtools/{{sample}}.idxstats.tsv",bucket=S3_BUCKET)),
            flagstats = S3.remote(expand("{bucket}/qc_reports/{{sample}}/qc_samtools/{{sample}}.flagstat.tsv",bucket=S3_BUCKET))
    log:    S3.remote(expand("{bucket}/logs/{{sample}}/qc_samtools.log",bucket=S3_BUCKET))
    threads: 1
    conda: "../wrappers/qc_samtools/env.yaml"
    script: "../wrappers/qc_samtools/script.py"


def mark_duplicates_ref(wildcards):
    input = {
        'bam': S3.remote(expand("{bucket}/mapped/{{sample}}.not_markDups.bam",bucket=S3_BUCKET)),
        'bai': S3.remote(expand("{bucket}/mapped/{{sample}}.not_markDups.bam.bai",bucket=S3_BUCKET)),
    }
    if config["umi_usage"] == "umi_concensus":
        input['ref'] = S3.remote(expand("cerit_sc_test/references/homsap/{ref}/index/BWA/{ref}.bwt",ref=config["reference"])[0]),
        input['lib_ROI'] = \
        S3.remote(expand("cerit_sc_test/references/homsap/GRCh37-p13/intervals/{lib_ROI}/{lib_ROI}.bed",ref=config["reference"],lib_ROI=config["lib_ROI"])[0]),
        input['fa'] = S3.remote(expand("cerit_sc_test/references/homsap/{ref}/seq/{ref}.fa",ref=config["reference"])[0]),
    return input


rule mark_duplicates:
    input: unpack(mark_duplicates_ref)
    output:
        bam=S3.remote(expand("{bucket}/mapped/{{sample}}.bam",bucket=S3_BUCKET)),
        bai=S3.remote(expand("{bucket}/mapped/{{sample}}.bam.bai",bucket=S3_BUCKET))
    log: S3.remote(expand("{bucket}/logs/{{sample}}/mark_duplicates.log",bucket=S3_BUCKET))
    threads: 8
    resources: mem=10
    params:
        mtx="qc_reports/{sample}/MarkDuplicates/{{sample}}.markDups_metrics.txt",
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
    if config["trim_adapters"] == True or config["trim_quality"] == True:
        preprocessed = "cleaned_fastq"
    else:
        preprocessed = "raw_fastq"
    input = { 'ref' : S3.remote(expand("cerit_sc_test/references/homsap/GRCh37-p13/index/BWA/GRCh37-p13.{ref}",ref=["amb","ann","bwt","pac","sa"]))}
    if read_pair_tags == [""]:
        input['sam'] = S3.remote(expand("cerit_sc_test/alignment_dna/{preprocessed}/{{sample}}.fastq.gz",preprocessed=preprocessed))
    else:
        input['r1'] = S3.remote(expand("cerit_sc_test/alignment_dna/{preprocessed}/{{sample}}_R1.fastq.gz",preprocessed=preprocessed))
        input['r2'] = S3.remote(expand("cerit_sc_test/alignment_dna/{preprocessed}/{{sample}}_R2.fastq.gz",preprocessed=preprocessed))
    return input

rule alignment_DNA:
    input: unpack(alignment_DNA_input),
    output:
        bam = S3.remote(expand("{bucket}/mapped/{{sample}}.not_markDups.bam",bucket=S3_BUCKET)),
        bai = S3.remote(expand("{bucket}/mapped/{{sample}}.not_markDups.bam.bai",bucket=S3_BUCKET))
    log: S3.remote(expand("{bucket}/logs/{{sample}}/alignment_DNA.log",bucket=S3_BUCKET))
    params: lib_name=config["library_name"]
    threads: 40
    conda: "../wrappers/alignment_DNA/env.yaml"
    script: "../wrappers/alignment_DNA/script.py"
