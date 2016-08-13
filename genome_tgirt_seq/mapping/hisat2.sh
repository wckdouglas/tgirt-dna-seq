#!/bin/bash

PROJECT_PATH=/scratch/02727/cdw2854/jurkatCells
TRIM_PATH=$PROJECT_PATH/
MAP_PATH=$PROJECT_PATH/bamFiles
REF_PATH=/corral-repl/utexas/2013lambowitz/Ref/hg19/hisat2/grch37_snp
REFERENCE=$REF_PATH/genome_snp
SPLIT_SITES=$REF_PATH/splicesites.tsv
SUFFIX=1P.fq.gz
THREADS=12
mkdir -p $MAP_PATH

for FQ1 in $TRIM_PATH/*$SUFFIX
do
	SAMPLENAME=$(basename ${FQ1/_$SUFFIX})
	FQ2=${FQ1/1P/2P}
#	echo flash --min-overlap=10 \
#			--to-stdout --threads=$THREADS \
#			$FQ1 $FQ2 \
	echo hisat2 --no-mixed --no-discordant \
			--threads $THREADS --known-splicesite-infile $SPLIT_SITES \
			-x $REFERENCE -1 $FQ1 -2 $FQ2 \
		\| samtools view -b -@ $THREADS \
		\| samtools sort -O bam -@ $THREADS  \
		\> $MAP_PATH/${SAMPLENAME}.bam

done
