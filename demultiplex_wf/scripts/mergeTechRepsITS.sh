#!/bin/bash

set -ex
​
for f in J_ITS_N J_ITS_R; do
    indir=/mnt/picea/projects/metaseq/nstreet/ITS/DeML_TP/$f
​
    mkdir -p $indir/DeML_pooled
​
    for ff in $(find $indir -name "*_r?.fq.gz"); do
	bn=$(basename $ff)
	id=${bn%_*}
	id2=${id##*_}
	id3=${id2%.*}
	echo $id3; done | sort | uniq | while read line;
    do
	cat $indir/demultiplex_${line}.1_r1.fq.gz $indir/demultiplex_${line}.2_r1.fq.gz $indir/demultiplex_${line}.3_r1.fq.gz > $indir/DeML_pooled/${line}_R1.fastq.gz
	cat $indir/demultiplex_${line}.1_r2.fq.gz $indir/demultiplex_${line}.2_r2.fq.gz $indir/demultiplex_${line}.3_r2.fq.gz > $indir/DeML_pooled/${line}_R2.fastq.gz
    done
done
