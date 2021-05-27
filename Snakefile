import os
import pandas as pd
import json
from snakemake.utils import min_version

min_version("5.18.0")

GLOBAL_REF_PATH = "/mnt/references/"

##### Config processing #####
# Folders
#
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

# Reference processing
#
if config["lib_ROI"] != "wgs":
    # setting reference from lib_ROI
    f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","lib_ROI"),)
    lib_ROI_dict = json.load(f)
    f.close()
    config["reference"] = [ref_name for ref_name in lib_ROI_dict.keys() if type(lib_ROI_dict[ref_name]) == "dictionary" and config["lib_ROI"] in lib_ROI_dict[ref_name].values()][0]


# setting organism from reference
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference"),)
reference_dict = json.load(f)
f.close()
config["organism"] = [organism_name for organism_name in reference_dict.keys() if type(reference_dict[organism_name]) == "dictionary" and config["reference"] in reference_dict[organism_name].values()][0]


# Samples
#
sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")

if config["lib_reverse_read_length"] == 0:
    read_pair_tags = [""]
else:
    read_pair_tags = ["_R1","_R2"]

wildcard_constraints:
    sample = "|".join(sample_tab.sample_name),
    read_pair_tag = "(_R.)?"


##### Target rules #####
rule all:
    input:  expand("mapped/{sample}.bam",sample = sample_tab.sample_name),
            expand("mapped/{sample}.bam.bai", sample = sample_tab.sample_name),

##### Modules #####

include: "rules/alignment_DNA.smk"