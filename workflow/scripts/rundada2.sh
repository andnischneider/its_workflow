#!/bin/bash -l

#SBATCH -A u2015005
#SBATCH -p core
#SBATCH -n 8
#SBATCH --mem=64000
#SBATCH -t 2-00:00:00
#SBATCH --mail-type=ALL

set -ex

module load bioinfo-tools R

fwd=$1
rv=$2
out=$3


echo "Running dada2"

Rscript /mnt/picea/projects/metaseq/nstreet/TP_v2/runDada2.R $fwd $rv $out
