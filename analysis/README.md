# Analysis scripts

## Contents

This folder contains R scripts to reproduce figures in the publication "Comparative fungal community analyses using metatranscriptomics and ITS-amplicon sequencing from Norway spruce" by Schneider et al.

## Instructions

Each of the scripts can be run independently, once **Prepare_ALL.R** has been fully executed, but the starting point for every *Figure_X.R* script is a clean R environment.

The name of the folders and scripts indicate which figure will be reproduced with the contained code. Before any of the figures can be plotted, the data has to be processed using the **Prepare_ALL.R** script

Additionally, some of the output files from the RNA-Seq metatranscriptomic workflow have to be copied into *data/RNA*. The following files are needed to reproduce all figures:

Assembly_2012.raw.tsv

gene_taxonomy.tsv

Assembly_2012.tpm.tsv

kingdom.Fungi.kos.raw.tsv

annotation_results.emapper.annotations
