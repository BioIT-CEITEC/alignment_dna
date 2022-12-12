from snakemake.utils import min_version

min_version("7.2.1")

configfile: "config.json"

module BR:
    snakefile: gitlab("bioroots/bioroots_utilities", path="bioroots_utilities.smk",branch="main")
    config: config

use rule * from BR as other_*

##### Config processing #####
#
GLOBAL_TMPD_PATH = config["globalTmpdPath"]
os.makedirs(GLOBAL_TMPD_PATH, exist_ok=True)

sample_tab = BR.load_sample()
read_pair_tags = BR.set_read_pair_tags()

# ##### Reference processing #####
# #

reference = BR.load_ref()
organism = BR.load_organism()
reference_directory = BR.reference_directory()


wildcard_constraints:
    sample = "|".join(sample_tab.sample_name),
    read_pair_tag = "(_R.)?"

#### Target rules #####

rule all:
    input:  BR.remote(expand("mapped/{sample}.bam",sample = sample_tab.sample_name)),
            BR.remote(expand("mapped/{sample}.bam.bai", sample = sample_tab.sample_name)),
            BR.remote("qc_reports/all_samples/alignment_DNA_multiqc/multiqc.html")
    output: BR.remote("completed.txt")
    shell: "touch {output}"



### Modules #####

include: "rules/alignment_DNA.smk"
