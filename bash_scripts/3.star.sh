#!/bin/bash

source activate /path/to/miniconda3/envs/rnaseq

REFPATH=/path/to/reference_genome/
INPATH=/path/to/01.Fastq/ #trim or orginal
OUTPATH=/path/to/03.bam/


cd /path/workdir/

# ls $INPATH/ > samplelist.txt

cat samplelist.txt | while read i
do 
    	echo "$i"
	
	STAR --runMode alignReads \
		--genomeDir $REFPATH/star_GRCh38/ \
	       	--readFilesIn $INPATH/$i.R1.fq.gz $INPATH/$i.R2.fq.gz \
       		--sjdbGTFfile $REFPATH/Homo_sapiens.GRCh38.110.chr.gtf  \
		--runThreadN 64 --readFilesCommand zcat \
		--outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix $OUTPATH/$i/$i. 

done


# note: carefully check that the star build on the same version of gtf as used in sjdbGTFfile