from snakemake.utils import min_version
import json
import os

min_version("7.2.1")

configfile: "config.json"
#
# module BR:
#     snakefile: gitlab("bioroots/bioroots_utilities", path="bioroots_utilities.smk",branch="main")
#     config: config
#
# use rule * from BR as other_*
#
# ##### Config processing #####
# #
# GLOBAL_REF_PATH = config["globalResources"]
# sample_tab = BR.load_sample()
# read_pair_tags = BR.set_read_pair_tags()[0]
# paired = BR.set_read_pair_tags()[1] # nahradit if not "is_paired in config:
#
# # ##### Reference processing #####
# # #
# # from boto3 import client
# # AWS_ID="acgt"
# # AWS_KEY="P84RsiL5TmHu0Ijd"
# # S3_BUCKET = "acgt"
# # ###############################x
#
# # FILE_TO_READ = "/resources/resources_info/lib_ROI.json"
# # client = client('s3',aws_access_key_id="acgt",aws_secret_access_key="P84RsiL5TmHu0Ijd",region_name="", endpoint_url="https://storage-elixir1.cerit-sc.cz")
# #
# # f = client.get_object(Bucket=S3_BUCKET, Key=FILE_TO_READ)
# # lib_ROI_dict = json.loads(f["Body"].read().decode('utf-8'))
# # config["reference"] = [ref_name for ref_name in lib_ROI_dict.keys() if isinstance(lib_ROI_dict[ref_name],dict) and config["lib_ROI"] in lib_ROI_dict[ref_name].keys()][0]
# #
# #
# #
# # FILE_TO_READ2 = "/resources/resources_info/reference.json"
# # ff = client.get_object(Bucket=S3_BUCKET, Key=FILE_TO_READ2)
# # reference_dict = json.loads(ff["Body"].read().decode('utf-8'))
# # config["organism"] = [organism_name.lower().replace(" ","_") for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]
# #
# #
# # reference_directory = os.path.join("/resources","organisms",config["organism"],config["reference"])
#
#
#
# BR.load_ref()
# BR.load_organism()
# reference_directory = BR.reference_directory()
#
#
# wildcard_constraints:
#     sample = "|".join(sample_tab.sample_name),
#     read_pair_tag = "(_R.)?"

##### Target rules #####
# rule all:
#     input:  BR.remote(expand("mapped/{sample}.bam",sample = sample_tab.sample_name)),
#             BR.remote(expand("mapped/{sample}.bam.bai", sample = sample_tab.sample_name)),
            #BR.remote("qc_reports/all_samples/alignment_DNA_multiqc/multiqc.html")

rule all:
    input: "text.txt",
           "text2.txt"
    # output: BR.remote("text.txt")
    # shell:"touch {output}"

##### Modules #####

include: "rules/alignment_DNA.smk"
