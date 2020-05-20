#!/bin/bash

set -ex

in=input_dir
out=output_dir

module load bioinfo-tools ITSx\1.1.2

ITSx -i $in/seqs.fasta -o $out/prefix -- cpu 4 --multi_thread T --preserve T --partial --minlen 50
