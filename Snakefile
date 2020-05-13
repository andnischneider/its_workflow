include: "src/common.smk"

rule all:
    """
    Collect the main outputs of the workflow
    """
    input:
        #expand(opj(config["results_dir"],"dada2","ITS","dada_{x}.rds"), x=["f","r"]),
        #expand(opj(config["results_dir"],"dada2","ITS","{y}.rds"), y=["mergers","seqtab"])
        expand(opj(config["results_dir"], "intermediate", "preprocess",
                 "prefilter", "{pool_id}_R{i}.fastq.gz"),
               pool_id = pools.keys(), i = [1,2])

rule fastq_dump:
    output:
        R1 = opj(config["data_dir"], "{sample_id}_R1.fastq.gz"),
        R2 = opj(config["data_dir"], "{sample_id}_R2.fastq.gz")
    log:
        opj("logs","{sample_id}.fastq_dump.log")
    params:
        data_dir = lambda wildcards, output: os.path.dirname(output.R1),
        acc = lambda wildcards: samples[wildcards.sample_id],
        spots = config["maxSpotId"]
    conda:
        "envs/sratools.yml"
    shell:
        """
        fastq-dump {params.spots} --split-3 --gzip -O {params.data_dir} \
            {params.acc} > {log} 2>&1
        mv {params.data_dir}/{params.acc}_1.fastq.gz {output.R1}
        mv {params.data_dir}/{params.acc}_2.fastq.gz {output.R2}
        """

rule pool:
    """Collect all samples with same sampleName but different runID"""
    input:
        lambda wildcards: pools[wildcards.pool_id][wildcards.R]
    output:
        opj(config["results_dir"], "intermediate", "pooled", "{pool_id}_{R}.fastq.gz"),
    message:
        "Pooling {wildcards.pool_id}_R{wildcards.R}.fastq.gz"
    run:
        if len(input) == 1:
            src = os.path.abspath(input[0])
            dst = os.path.abspath(output[0])
            shell("ln -s {src} {dst}")
        else:
            shell("gunzip -c {input} | gzip -c > {output}")

rule prefilter:
    """Runs the ITS prefiltering step for ambiguous bases"""
    input:
        R1 = opj(config["results_dir"], "intermediate", "pooled",
                 "{pool_id}_R1.fastq.gz"),
        R2 = opj(config["results_dir"], "intermediate", "pooled",
                 "{pool_id}_R2.fastq.gz")
    output:
        R1 = opj(config["results_dir"], "intermediate", "preprocess",
                 "prefilter", "{pool_id}_R1.fastq.gz"),
        R2 = opj(config["results_dir"], "intermediate", "preprocess",
                 "prefilter", "{pool_id}_R2.fastq.gz")
    log:
        opj("logs", "{pool_id}", "prefilter.log")
    threads: 4
    params:
        truncLen = config["dada2"]["truncLen"],
        maxN = config["dada2"]["maxN"],
        maxEE = config["dada2"]["maxEE"],
        rm_phix = config["dada2"]["rm_phix"],
        minLen = config["dada2"]["minLen"],
        truncQ = config["dada2"]["truncQ"]
    conda:
        "envs/dada2.yml"
    script:
        "src/R/prefilter.R"

rule generate_revcomp:
    output:
        fwd_rc = temp(opj(config["results_dir"], "preprocess", "fwd_rc")),
        rev_rc = temp(opj(config["results_dir"], "preprocess", "rev_rc"))
    params:
        fwd = config["cutadapt"]["FWD"],
        rev = config["cutadapt"]["REV"]
    conda:
        "envs/biopython.yml"
    shell:
        """
        python src/revcomp.py {params.fwd} > {output.fwd_rc}
        python src/revcomp.py {params.rev} > {output.rev_rc}
        """

rule cut_ITS_primers:
    """Removes primers from ITS data using cutadapt"""
    input:
        R1 = opj(config["results_dir"], "preprocess", "prefilter", "{pool_id}_R1.fastq.gz"),
        R2 = opj(config["results_dir"], "preprocess", "prefilter", "{pool_id}_R2.fastq.gz"),
        fwd_rc = opj(config["results_dir"], "preprocess", "fwd_rc"),
        rev_rc = opj(config["results_dir"], "preprocess", "rev_rc")
    output:
        R1 = opj(config["results_dir"], "preprocess", "cutadapt", "R1","{sample}_R1.fastq.gz"),
        R2 = opj(config["results_dir"], "preprocess", "cutadapt", "R2","{sample}_R2.fastq.gz")
    log:
        "logs/cutadapt/{sample}.cutadapt.log"
    params:
        FWD=config["cutadapt"]["FWD"],
        REV=config["cutadapt"]["REV"]
    shell:
        """
        A=$(cat {input.fwd_rc})
        a=$(cat {input.rev_rc})
        
        cutadapt -g {params.FWD} -a $a -G {params.REV} -A $A \
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
