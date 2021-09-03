import os
import pandas as pd
import json
import boto3
from snakemake.utils import min_version
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider

min_version("5.18.0")
configfile: "config.json"

AWS_ID="DP2K6W7WY0WJGJCD7OGU"
AWS_KEY=[REDIGED]
S3_BUCKET = "cerit_sc_test"

S3 = S3RemoteProvider(host="https://s3.cl2.du.cesnet.cz",access_key_id=AWS_ID, secret_access_key=AWS_KEY)
client = boto3.client('s3', aws_access_key_id=AWS_ID, aws_secret_access_key=AWS_KEY, region_name="", endpoint_url="https://s3.cl2.du.cesnet.cz") 

GLOB_REF_PATH= "references/"
REF_INFO_PATH = GLOB_REF_PATH + "reference_info/"


# Reference processing
#
if config["lib_ROI"] != "wgs":
    path_to_read = REF_INFO_PATH + "lib_ROI.json"
    obj = client.get_object(Bucket=S3_BUCKET,Key=path_to_read)
    obj_data = obj["Body"].read()
    # setting reference from lib_ROI
    lib_ROI_dict = json.loads(obj_data)
    config["reference"] = [ref_name for ref_name in lib_ROI_dict.keys() if isinstance(lib_ROI_dict[ref_name],dict) and config["lib_ROI"] in lib_ROI_dict[ref_name].keys()][0]


# setting organism from reference
#
path_to_read = REF_INFO_PATH + "reference_fix.json"
obj = client.get_object(Bucket=S3_BUCKET,Key=path_to_read)
obj_data = obj["Body"].read()
reference_dict = json.loads(obj_data)
config["organism"] = [organism_name.lower().replace(" ","_") for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].values()][0]


##### Config processing #####
# Folders
#
reference_directory = os.path.join("cerit_sc_test/",GLOB_REF_PATH,config["organism"],config["reference"])

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
    input:  S3.remote(expand("{bucket}/mapped/{sample}.bam",sample = sample_tab.sample_name,bucket=S3_BUCKET)),
            S3.remote(expand("{bucket}/mapped/{sample}.bam.bai", sample = sample_tab.sample_name,bucket=S3_BUCKET)),
            S3.remote(expand("{bucket}/qc_reports/all_samples/alignment_DNA_multiqc/multiqc.html",bucket=S3_BUCKET))

##### Modules #####

include: "rules/alignment_DNA.smk"
