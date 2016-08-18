#!/bin/bash

PROJECTPATH=$WORK/$USER/genomeDNA
DATAPATH=$PROJECTPATH/rawData/clipped
INDEXPATH=$REF/hg19/Sequence/WholeGenomeFasta
INDEX=$INDEXPATH/genome.fa
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
