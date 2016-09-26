#!/bin/bash

PROJECTPATH=$WORK/cdw2854/plasmaDNA
DATAPATH=$PROJECTPATH/rawData
DATAPATH=${DATA}/SA16172/dna_plasma_genome/clipped
INDEXPATH=$REF/GRCh38/hg38_rDNA
INDEX=$INDEXPATH/genome_rDNA
PYTHON=$(which python)
CORES=12

for FQ in `ls $DATAPATH/*R1_001.fastq.gz `
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
