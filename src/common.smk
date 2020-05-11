configfile: "config.yml"
from Bio.Seq import reverse_complement
from os.path import join as opj
from src.parse_samples import get_samples

wildcard_constraints:
    run = "\d+"

if config["ITS"]["FWD"]:
    config["ITS"]["FWD_RC"] = reverse_complement(config["ITS"]["FWD"])
if config["ITS"]["REV"]:
    config["ITS"]["REV_RC"] = reverse_complement(config["ITS"]["REV"])

# Use only a small set of reads if testing the workflow
if config["test"]:
    config["maxSpotId"] = "-X 1000"
else:
    config["maxSpotId"] = ""
samples, pools, config = get_samples(config["sample_file"], config)

if config["source"] == "sra":
    ruleorder: fastq_dump > download_sample
elif config["source"] == "url":
    ruleorder: download_sample > fastq_dump

def get_files(wildcards):
    suffix = "_{}.fastq.gz".format(wildcards.R)
    files = []
    input_dir = config["data_dir"]
    for runID in samples[wildcards.sample].keys():
        gene = samples[wildcards.sample][runID]['gene']
        files.append(opj(input_dir,"{}{}".format(samples[wildcards.sample][runID][gene], suffix)))
    return sorted(files)