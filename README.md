# Amplicon data preprocessing and analysis workflow for Schneider et al.

This repository contains code to analyse the ITS amplicon sequencing data with DADA2 for further comparisons with RNA-seq data from the same samples. It uses data collected for Haas et al (2018), but using dada2 and Swarm clustering instead of OTU clustering.

## Repository contents

The repository consists of two units:

### 1 demultiplex_wf -> download and demultiplex raw data

The folder contains instructions on how to run the workflow as a docker container. It will download the raw data from the ENA and demultiplex the sequences into files per sample.

### 2 workflow -> ITS amplicon sequencing data preprocessing workflow

The folder contains a snakemake workflow that will reproduce the preprocessing of the demultiplexed ITS amplicon sequencing data as used in the study.
