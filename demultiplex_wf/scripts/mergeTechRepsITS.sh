#!/bin/bash

indir=$1
outdir=$2

find $indir -name "*_r1.fq.gz" | sed 's/.[1-3]_r1.fq.gz//g' | \
    sort | uniq -c | sed 's/^ \+//g' | awk -F ' ' '{if ($1==3) print $2}' | \
    while read id;
    do
        bn=$(basename $id)
        gunzip -c $id.1_r1.fq.gz $id.2_r1.fq.gz $id.3_r1.fq.gz | gzip -c > $outdir/${bn}_R1.fastq.gz
        gunzip -c $id.1_r2.fq.gz $id.2_r2.fq.gz $id.3_r2.fq.gz | gzip -c > $outdir/${bn}_R2.fastq.gz
    done
