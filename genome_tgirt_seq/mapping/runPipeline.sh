#!/bin/bash

PROJECTPATH=$WORK/cdw2854/ecoli_genome
#PROJECTPATH=$WORK/cdw2854/genomeDNA
DATAPATH=$PROJECTPATH/rawData/k12
INDEXPATH=$REF/Ecoli
INDEX=$INDEXPATH/k12_mg1655.fa
#INDEX=$REF//hg19/Sequence/WholeGenomeFasta/genome.fa
PYTHON=$(which python)
ADAPTORS=adaptors.fa
CORES=12

for FQ in $DATAPATH/*R1_001.fastq.gz
do
	SAMPLENAME=$(basename $FQ)
	echo $PYTHON DNAmapping.py --fq1=$FQ  \
		--outdir=$PROJECTPATH --index=$INDEX \
		--threads=$CORES --adaptor=$ADAPTORS 
done 
