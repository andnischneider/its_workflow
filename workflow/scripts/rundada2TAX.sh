#!/bin/bash -l

#SBATCH -A u2015005
#SBATCH -p core
#SBATCH -n 8
#SBATCH --mem=64000
#SBATCH -t 2-00:00:00
#SBATCH --mail-type=ALL

set -ex

module load bioinfo-tools R

echo "assigning taxonomy"

Rscript rundada2TAX.R $1 $2 $3
