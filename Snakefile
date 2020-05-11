include: "src/common.smk"

rule all:
    """
    Collect the main outputs of the workflow
    """
    input:
        expand(opj(config["results_dir"],"dada2","ITS","dada_{x}.rds"), x=["f","r"]),
        expand(opj(config["results_dir"],"dada2","ITS","{y}.rds"), y=["mergers","seqtab"])


rule fastq_dump:
    output:
        R1 = opj(config["data_dir"], "{sample_id}_1.fastq.gz"),
        R2 = opj(config["data_dir"], "{sample_id}_2.fastq.gz")
    log:
        opj("logs","{sample_id}.fastq_dump.log")
    params:
        data_dir = lambda wildcards, output: os.path.dirname(output.R1),
        acc = lambda wildcards: samples[wildcards.sample_id],
        spots = config["maxSpotId"]
    shell:
        """
        fastq-dump {params.spots} --split-3 --gzip -O {params.data_dir} \
            {params.acc} > {log} 2>&1
        mv {params.data_dir}/{params.acc}_1.fastq.gz {output.R1}
        mv {params.data_dir}/{params.acc}_2.fastq.gz {output.R2}
        """

rule download_sample:
    output:
        R1 = opj(config["data_dir"], "{sample_id}_1.fastq.gz"),
        R2 = opj(config["data_dir"], "{sample_id}_2.fastq.gz")
    log:
        opj("logs","{sample_id}.curl.log")
    params:
        data_dir = lambda wildcards, output: os.path.dirname(output.R1),
        url = lambda wildcards: samples[wildcards.sample_id]
    shell:
        """
        curl -L -o {params.data_dir}/{wildcards.sample_id}.tar.gz \
            {params.url} > {log} 2>&1
        tar -C {params.data_dir} xzf {params.data_dir}/{wildcards.sample_id}.tar.gz
        rm {params.data_dir}/{wildcards.sample_id}.tar.gz
        """

rule pool:
    """Collect all samples with same sampleName but different runID"""
    input:
        R1 = lambda wildcards: pools[wildcards.pool_id]["R1"],
        R2 = lambda wildcards: pools[wildcards.pool_id]["R2"]
    output:
        R1 = opj(config["results_dir"], "intermediate", "pooled", "{pool_id}_R1.fastq.gz"),
        R2 = opj(config["results_dir"], "intermediate", "pooled", "{pool_id}_R2.fastq.gz")
    message:
        "Pooling samples for {wildcards.pool_id}"
    shell:
        """
        gunzip -c {input.R1} | gzip -c > {output.R1}
        gunzip -c {input.R2} | gzip -c > {output.R2}
        """


rule prefilter:
    """Runs the ITS prefiltering step for ambiguous bases"""
    input:
        R1 = opj(config["results_dir"], "intermediate", "pooled",
                 "{pool_id}_R1.fastq.gz"),
        R2 = opj(config["results_dir"], "intermediate", "pooled",
                 "{pool_id}_R2.fastq.gz")
    output:
        R1 = opj(config["results_dir"], "intermediate", "preprocess", "ITS",
                 "filtN", "{pool_id}_R1.fastq.gz"),
        R2 = opj(config["results_dir"], "intermediate", "preprocess", "ITS",
                 "filtN", "{pool_id}_R2.fastq.gz"),
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
        R1=opj(config["results_dir"],"preprocess","ITS","R1","{sample}_R1.fastq.gz"),
        R2=opj(config["results_dir"],"preprocess","ITS","R2","{sample}_R2.fastq.gz")
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

rule stage_16S:
    input:
        opj(config["results_dir"],"preprocess","staged","{sample}_{R}.fastq.gz")
    output:
        opj(config["results_dir"],"preprocess","16S","{R}","{sample}_{R}.fastq.gz")
    shell:
        """
        ln -s {input} {output}
        """

def get_dada_input(wildcards):
    files = []
    for sample in samples.keys():
        for runID in samples[sample].keys():
            if list(samples[sample][runID].keys())[0] == wildcards.gene:
                files.append(opj(config["results_dir"],"preprocess",wildcards.gene,"R1",sample+"_R1.fastq.gz"))
                files.append(opj(config["results_dir"],"preprocess",wildcards.gene,"R2",sample+"_R2.fastq.gz"))
    return files

rule dada2:
    """Calls DADA2 Rscript with trimmed input"""
    input:
         get_dada_input
    output:
        expand(opj(config["results_dir"],"dada2","{{gene}}","dada_{x}.rds"), x=["f","r"]),
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
        rundada2.R {params.fw_dir} {params.rv_dir} {params.out_dir} {threads} > {log} 2>&1
        """
