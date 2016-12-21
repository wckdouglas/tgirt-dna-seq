#!/bin/bash

PROJECTPATH=$SCRATCH/plasmaDNA
DATAPATH2=$PROJECTPATH/rawData/split
DATAPATH1=$PROJECTPATH/splittedFastq
#DATAPATH=$DATA/SA16172/dna_plasma_genome/clipped_id
INDEXPATH=$SCRATCH/GRCh38/hg38_rDNA
INDEX=$INDEXPATH/genome_rDNA
PYTHON=$(which python)
CORES=3

for FQ in $DATAPATH1/*R1_001.fastq.gz  $DATAPATH2/*R1*.fastq.gz
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
