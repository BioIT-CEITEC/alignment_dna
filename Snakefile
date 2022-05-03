from snakemake.utils import min_version
import os

min_version("7.2.1")

config["computing_type"] = "kubernetes"

configfile: "config.json"

module BR:
    snakefile: gitlab("bioroots/bioroots_utilities", path="bioroots_utilities.smk",branch="main")
    config: config

use rule * from BR as other_*

##### Config processing #####
#
GLOBAL_REF_PATH = config["globalResources"]
sample_tab = BR.load_sample()

read_pair_tags = BR.set_read_pair_tags()[0]
paired = BR.set_read_pair_tags()[1] # nahradit if not "is_paired in config:

# ##### Reference processing #####
# #

BR.load_ref()
BR.load_organism()
reference_directory = BR.reference_directory()



wildcard_constraints:
    sample = "|".join(sample_tab.sample_name),
    read_pair_tag = "(_R.)?"

#### Target rules #####
# rule all:
#     input:  BR.remote(expand("mapped/{sample}.bam",sample = sample_tab.sample_name)),
#             BR.remote(expand("mapped/{sample}.bam.bai", sample = sample_tab.sample_name)),
#             BR.remote("qc_reports/all_samples/alignment_DNA_multiqc/multiqc.html")
rule all:
    input:  bam = BR.remote(expand("mapped/{sample}.not_markDups.bam",sample = sample_tab.sample_name)),
            bai = BR.remote(expand("mapped/{sample}.not_markDups.bam.bai",sample = sample_tab.sample_name)),



##### Modules #####

include: "rules/alignment_DNA.smk"


