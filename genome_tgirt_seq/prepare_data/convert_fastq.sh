#!/bin/bash

FQ_PATH=/stor/work/Lambowitz/cdw2854/ecoli_genome/rawData/k12
FASTQ64=$FQ_PATH/raw_sra

for FQ in $FASTQ64/*.fastq.gz
do
	SAMPLENAME=$(basename ${FQ%.fastq.gz})
	END=$(echo $SAMPLENAME | cut -d'_' -f2)
	SAMPLENAME=$(echo $SAMPLENAME | cut -d'_' -f1)
	echo zcat $FQ \
		\| fastq_quality_converter \
	    \| gzip \
		\> $FQ_PATH/${SAMPLENAME}_R${END}_001.fastq.gz	
done
