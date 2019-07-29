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
    input:
        ""
    output:
        ""
    log:
        "logs/dada2/"
    script:
        """
        rundada2.R
        """
