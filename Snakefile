rule all:
    """
    Collect the main outputs of the workflow
    """
    input:
        ""

rule cut_ITS_primers:
    input:
        config["input_ITS"]
    output:
        "data/ITS/cutadapt"
    log:
        "logs/..."
    script:
        """
        cutITS.R config["Forward"] config["Reverse"]
        """

rule filter_trim:
    input:
        config["input_16S"] if config["16S"] else "OUTPUT from cutting"
    output:
        ""
    log:
        "logs/"
    script:
        """
        filtertrim.R
        """

rule dada2:
    """Calls DADA2 Rscript with trimmed input"""
    input:
        R1=expand("results/preprocess/{{gene}}/R1/{sample}_R1.trimmed.fastq.gz", sample=samples.keys()),
        R2=expand("results/preprocess/{{gene}}/R2/{sample}_R1.trimmed.fastq.gz", sample=samples.keys())
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
    threads: 8
    shell:
        """
        rundada2.R {input.fw_dir} {input.rv_dir} {params.out_dir}} {threads} > {log} 2>&1
        """
