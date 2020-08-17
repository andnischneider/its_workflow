# Amplicon data preprocessing and analysis workflow for Schneider et al.

This repository contains code to analyse the ITS amplicon sequencing data with DADA2 for further comparisons with RNA-seq data from the same samples. It uses data collected for Haas et al (2018), but using dada2 and Swarm clustering instead of OTU clustering.

## Repository contents

The repository consists of two units that are run separately:

### 1 demultiplex_wf -> download and demultiplex raw data

The folder contains instructions on how to run the workflow as a docker container. It will download the raw data from the ENA and demultiplex the sequences into files per sample, as well as concatenate technical replicates.

### 2 workflow -> ITS amplicon sequencing data preprocessing workflow

The folder contains a snakemake workflow that will reproduce the preprocessing of the demultiplexed ITS amplicon sequencing data as used in the study.

#### How to use it

The workflow is run through snakemake, from the root folder of the repository (where this readme sits). To continue 
with the demultiplexed data, we need to move it from the _demultiplex_wf_ subfolder.

```bash
mkdir $(pwd)/data
mv $(pwd)/demultiplex_wf/data/ $(pwd)/data
```
Next we will create a conda environment ("its_wf") needed to execute the snakemake workflow,
and then activate it:

```bash
conda env create -n its_wf -f $(pwd)/environment.yml
conda activate its_wf
```

Once the conda environment has been created successfully, we can execute the workflow with the 
following command:

```bash
snakemake -s $(pwd)/workflow/Snakefile -pr -j 4 --use-conda
```

This will output the final count matrix and other results (such as sequences for every Swarm OTU and taxonomic assignments)
into the results/ folder.

### 3 analysis -> R scripts to analyse ITS and RNA data

The folder contains scripts to reproduce the figures in the publication.
