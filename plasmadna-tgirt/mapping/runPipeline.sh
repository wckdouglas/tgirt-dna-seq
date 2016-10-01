#!/bin/bash

<<<<<<< HEAD
PROJECTPATH=$SCRATCH/plasmaDNA
DATAPATH=$PROJECTPATH/rawData/split
INDEXPATH=$REF/GRCh38/hg38_rDNA
INDEX=$INDEXPATH/genome_rDNA
PYTHON=$(which python)
CORES=8

for FQ in `ls $DATAPATH/*R1*.fastq.gz $DATAPATH/*R1_001.fq.gz | grep 'PD\|SR\|NA'|grep -v 'DB\|NT\|RN' `
do
	SAMPLENAME=$(basename $FQ)
	if [[ $SAMPLENAME == SRR* ]]
	then
		ADAPTORS=TruSeq2-PE.fa
	else
		ADAPTORS=adaptors.fa
	fi
	echo $PYTHON DNAmapping.py \
		--fq1=$FQ \
		--outdir=$PROJECTPATH \
		--index=$INDEX \
		--threads=$CORES \
		--adaptor=$ADAPTORS 
done 
