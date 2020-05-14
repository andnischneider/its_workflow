# Demultiplex Workflow
This is a small snakemake workflow used to pull ITS sequences from SRA and 
demultiplex them using [deML](https://github.com/grenaud/deML). Its purpose is
to reproduce the input data for the ITS amplicon workflow used for the paper
by [Schneider et al 2020]().

The raw data is deposited at SRA under experiment accessions:

| accession  | Sample name | 
| ---------- | ----------- |
| ERX2087084 | ITS_roots   |
| ERX2087083 | ITS_needles |

The workflow downloads the sequence data and corresponding index files, then 
runs a specific version of deML to demultiplex the files. Finally, technical 
replicates are pooled.

## How to use it

The full workflow is run as a Docker container. In order to generate the 
demultiplexed and pooled files in folders `data/ITS_root` and `data/ITS_needle`
on your computer you can run:

```bash
docker run -v $(pwd)/data:/analysis/results/DeML_pooled natstreetlab/amplicon_wf:deml
```

If you also want to generate the QC reports of the demultiplexed sequences you 
can run:

```bash
docker run -v $(pwd)/qc:/analysis/results/report -v $(pwd)/data:/analysis/results/DeML_pooled natstreetlab/amplicon_wf:deml
```

This produces two folders `qc/DeML` and `qc/DeML_pooled` with QC reports for
the demultiplexed sequences before and after pooling, respectively.

The full workflow may take several hours depending on your internet connection
and machine hardware. If you want to test the workflow you can have it download 
a subset of sequences from SRA by passing `--config test=True` to 
the `docker run` command, _e.g._:

```bash
docker run -v $(pwd)/qc:/analysis/results/report -v $(pwd)/data:/analysis/results/DeML_pooled natstreetlab/amplicon_wf:deml --config test=True
```