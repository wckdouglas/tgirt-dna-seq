#!/bin/bash

PROJECTPATH=$SCRATCH/plasmaDNA
PROJECTPATH=$WORK/cdw2854/plasmaDNA
DATAPATH=$PROJECTPATH/splittedFastq
DATAPATH=/stor/work/Lambowitz/Data/NGS/JA17030/umi2id
#DATAPATH=$DATA/SA16172/dna_plasma_genome/clipped_id
INDEXPATH=$REF/GRCh38/hg38_rDNA
INDEX=$INDEXPATH/genome_rDNA
PYTHON=$(which python)
CORES=14

for FQ in $DATAPATH/*R1_001.fastq.gz
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
