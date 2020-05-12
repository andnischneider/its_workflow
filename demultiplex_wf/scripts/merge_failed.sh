#!/bin/bash

indir=$1
outdir=$2
strip=$3

for f in $(find $indir -name "*_r[1,2].fq.gz"); do
    n=$(basename $f | sed "s/$strip//g")
    mv $f ${f//.fq/.pass.fq}
    cat ${f//.fq/.pass.fq} ${f//.fq/.fail.fq} > $outdir/$n
done