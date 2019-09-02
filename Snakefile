configfile: "config.yml"

from Bio.Seq import reverse_complement
from os.path import join as opj
from src.parse_samples import get_samples

if config["ITS"]["FWD"]:
    config["ITS"]["FWD_RC"] = reverse_complement(config["ITS"]["FWD"])
if config["ITS"]["REV"]:
    config["ITS"]["REV_RC"] = reverse_complement(config["ITS"]["REV"])

samples = get_samples(config["sample_file"], config["pool"])

def workflow_results(samples):
    return []
#    for sampleName in samples.keys():
#        for

rule all:
    """
    Collect the main outputs of the workflow
    """
    input:
        workflow_results(samples)


rule download_samples:
    output:
        R1=opj(config["data_dir"],"{sample}_R1.fastq.gz"),
        R2=opj(config["data_dir"],"{sample}_R2.fastq.gz")
    params:
        data_dir=opj(config["data_dir"])
    shell:
        """
        fastq-dump --split-3 --gzip -O {params.data_dir} {wildcards.sample}
        """

rule demultiplex:
    input:
        R1=opj(config["data_dir"],"{sample}_R1.fastq.gz"),
        R2=opj(config["data_dir"],"{sample}_R2.fastq.gz")
    output:
        R1=opj(config["results_dir"],"preprocess","deML","{sample}_R1.fastq.gz"),
        R2=opj(config["results_dir"],"preprocess","deML","{sample}_R2.fastq.gz")

def get_files(wildcards):
    suffix = "_{}.fastq.gz".format(wildcards.R)
    files = []
    if config["demultiplex"]:
        input_dir = opj(config["results_dir"],"preprocess","deML")
    else:
        input_dir = config["data_dir"]
    for runID in samples[wildcards.sample].keys():
        gene = list(samples[wildcards.sample][runID].keys())[0]
        files.append(opj(input_dir,"{}{}".format(samples[wildcards.sample][runID][gene], suffix)))
    return sorted(files)

rule handle_replicates:
    """Collect all samples with same sampleName but different runID"""
    input:
        get_files
    output:
        opj(config["results_dir"],"preprocess","staged","{sample}_{R}.fastq.gz")
    run:
        if len(input) == 1:
            shell("ln -s {input[0]} {output[0]}")
        else:
            shell("gunzip -c {input} | gzip -c > {output}")


rule prefilter:
    """Runs the ITS prefiltering step for ambiguous bases"""
    input:
        R1=opj(config["results_dir"],"preprocess","staged","{sample}_R1.fastq.gz"),
        R2=opj(config["results_dir"],"preprocess","staged","{sample}_R2.fastq.gz")
    output:
        R1=opj(config["results_dir"],"preprocess","ITS","filtN","{sample}_R1.fastq.gz"),
        R2=opj(config["results_dir"],"preprocess","ITS","filtN","{sample}_R2.fastq.gz")
    shell:
        """
        Rscript src/R/prefilter.R {input.R1} {input.R2} {output.R1} {output.R2}
        """

rule cut_ITS_primers:
    """Removes primers from ITS data using cutadapt"""
    input:
        R1=opj(config["results_dir"],"preprocess","ITS","filtN","{sample}_R1.fastq.gz"),
        R2=opj(config["results_dir"],"preprocess","ITS","filtN","{sample}_R2.fastq.gz")
    output:
        R1=opj(config["results_dir"],"preprocess","ITS","R1","{sample}_R1.trimmed.fastq.gz"),
        R2=opj(config["results_dir"],"preprocess","ITS","R2","{sample}_R2.trimmed.fastq.gz")
    log:
        "logs/cutadapt/{sample}.cutadapt.log"
    params:
        FWD=config["ITS"]["FWD"],
        REV=config["ITS"]["REV"],
        FWD_RC=config["ITS"]["FWD_RC"],
        REV_RC=config["ITS"]["REV_RC"]
    shell:
        """
        cutadapt -g {params.FWD} -a {params.REV_RC} -G {params.REV} -A {params.FWD_RC} \
         -n 2 -o {output.R1} -p {output.R2} {input.R1} {input.R2} > {log} 2>&1
        """

# rule filter_trim:
#     input:
#         config["input_16S"] if config["16S"] else "OUTPUT from cutting"
#     output:
#         ""
#     log:
#         "logs/"
#     script:
#         """
#         filtertrim.R
#         """

rule dada2:
    """Calls DADA2 Rscript with trimmed input"""
    input:
        R1=expand("results/preprocess/{{gene}}/R1/{sample}_R1.trimmed.fastq.gz", sample=samples.keys()),
        R2=expand("results/preprocess/{{gene}}/R2/{sample}_R2.trimmed.fastq.gz", sample=samples.keys())
    output:
        expand(opj(config["results_dir"],"dada2","{{gene}}","dada_{x}.rds"), x=["f,r"]),
        expand(opj(config["results_dir"],"dada2","{{gene}}","{y}.rds"), y=["mergers","seqtab"])
    log:
        "logs/dada2/{gene}.dada2.log"
    params:
        fw_dir=opj(config["results_dir"],"preprocess","{gene}","R1"),
        rv_dir=opj(config["results_dir"],"preprocess","{gene}","R2"),
        out_dir=opj(config["results_dir"],"dada2","{gene}")
    resources:
        runtime = lambda wildcards, attempt: attempt**2*60*4
    threads: config["threads"]
    shell:
        """
        rundada2.R {input.fw_dir} {input.rv_dir} {params.out_dir}} {threads} > {log} 2>&1
        """
