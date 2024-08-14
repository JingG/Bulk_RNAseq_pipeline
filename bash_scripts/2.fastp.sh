#!/bin/bash

## This step is optional, check on the fastQC results first. 
## Run this step if the adapters/polyA/polyG are siginficant.

source activate /path/to/miniconda3/envs/rnaseq

INPATH=/path/to/01.Fastq/
OUTPATH=/path/to/01-2.Trimmed_Fastq/

cd /path/workdir/
mkdir $OUTPATH

for j in $INPATH/*
do
    i=$(basename "$j")
    # out=`echo $i| cut -d'.' -f2`
    echo "$i"
    fastp -w 1 -i $j/*_1.fq.gz -I $j/*_2.fq.gz -o $OUTPATH/$i.R1.fq.gz -O $OUTPATH/$i.R2.fq.gz &
done


