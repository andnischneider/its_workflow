#!/bin/bash

set -ex


#Usage: $0 threads $input_dir $output_dir

threads=$1
in=$2
out=$3


module load bioinfo-tools Swarm\2.2.2

swarm -t $threads -d 3 -z --output-file $out/results.txt --seeds $out/seeds.fasta $in/seqs_sum.fasta
