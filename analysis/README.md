# Analysis scripts

## Contents

This folder contains R scripts to reproduce figures in the publication "Comparative fungal community analyses using metatranscriptomics and ITS-amplicon sequencing from Norway spruce" by Schneider et al.

## Instructions

Before the R scripts can be run, these three repositories have to be pulled (this assumes your starting point is in the root of the repository).

```bash
mkdir src

cd src

git clone https://github.com/loalon/Rtoolbox.git

git clone https://github.com/UPSCb/UPSCb-common.git
```


The name of the folders and scripts indicate which figure will be reproduced with the contained code. Before any of the figures can be plotted, the data has to be processed using the **Prepare_ALL.R** script

Additionally, some of the output files from the RNA-Seq metatranscriptomic workflow have to be copied into *data/RNA*. The following files are needed to reproduce all figures:

Assembly_2012.raw.tsv

gene_taxonomy.tsv

Assembly_2012.tpm.tsv

kingdom.Fungi.kos.raw.tsv
