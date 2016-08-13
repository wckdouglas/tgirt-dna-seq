#!/bin/bash

DATAPATH=/scratch/cdw2854/plasmaDNA/data

for FASTQ in $DATAPATH/SRR*gz
do 
	SAMPLENAME=${FASTQ%.fastq.gz}
	echo zcat $FASTQ \
		\| fastq_quality_converter \
		\| gzip \
		\> $SAMPLENAME'P.fastq.gz'
done
