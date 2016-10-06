#!/bin/bash

PROJECTPATH=$WORK/cdw2854/ecoli_genome
DATAPATH=$PROJECTPATH/rawData
INDEXPATH=$REF/Ecoli
INDEX=$INDEXPATH/b_strain.fa
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
