#!/bin/bash


FASTQC=/path/to/software/FastQC/
INPATH=/path/to/01.Fastq/
OUTPATH=/path/to/02.FastQC/

cd /path/workdir/
mkdir $OUTPATH

for j in $INPATH/*.fq.gz
do
    i=$(basename "$j")
    # out=`echo $i| cut -d'.' -f2`
    echo "$i"
    $FASTQC/fastqc $j  -o $OUTPATH/ &
done


