#!/bin/bash

PROJECTPATH=$SCRATCH/plasmaDNA
DATAPATH=$PROJECTPATH/rawData
INDEXPATH=$REF/GRCh38
INDEX=$INDEXPATH/genome
PYTHON=$(which python)
CORES=24

for FQ in `ls $DATAPATH/*R1_001.fastq.gz $DATAPATH/*R1_001.fq.gz | grep 'GJ\|RNA\|ref\|error' -v`
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
