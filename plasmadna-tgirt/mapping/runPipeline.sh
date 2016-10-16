#!/bin/bash

PROJECTPATH=$WORK/cdw2854/plasmaDNA
DATAPATH=$PROJECTPATH/rawData/split
DATAPATH=$DATA/SA16172/dna_plasma_genome/clipped_id
INDEXPATH=$REF/GRCh38/hg38_rDNA
INDEX=$INDEXPATH/genome_rDNA
PYTHON=$(which python)
CORES=12

for FQ in $DATAPATH/*R1*.fq.gz 
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
