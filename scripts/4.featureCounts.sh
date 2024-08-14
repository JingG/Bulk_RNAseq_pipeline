#!/bin/bash

source activate /path/to/miniconda3/envs/rnaseq


REFPATH=/path/to/reference_genome/
INPATH=/path/to/03.bam/ 
OUTPATH=/path/to/04.fc/

cd /path/workdir/
mkdir $OUTPATH

for j in $INPATH/*
do 
	i=$(basename "$j")
    	# out=`echo $i| cut -d'.' -f2`
    	echo "$i"

	featureCounts -t exon -p -g gene_id --fraction -M \
		-a $REFPATH/Homo_sapiens.GRCh38.110.chr.gtf \
		-o $OUTPATH/$i.count $j/$i.Aligned.sortedByCoord.out.bam &
done


