#!/bin/bash

set -ex

fw=/mnt/picea/projects/metaseq/nstreet/TP_v2/data/16S/Needles/DeML_pooled/filtered_F
rv=/mnt/picea/projects/metaseq/nstreet/TP_v2/data/16S/Needles/DeML_pooled/filtered_R

mkdir -p /mnt/picea/projects/metaseq/nstreet/TP_v2/logs/dada2/

sbatch -A u2015005 -e /mnt/picea/projects/metaseq/nstreet/TP_v2/logs/dada2/16s_needles.err -o /mnt/picea/projects/metaseq/nstreet/TP_v2/logs/dada2/16s_needles.out /mnt/picea/projects/metaseq/nstreet/TP_v2/submitdada2.sh $fw $rv
