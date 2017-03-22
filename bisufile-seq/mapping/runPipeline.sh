#!/bin/bash

PROJECTPATH=/stor/work/Lambowitz/cdw2854/bisufite_seq
DATAPATH=$PROJECTPATH/raw_Data
INDEXPATH=$REF/hg19/Sequence/WholeGenomeFasta
INDEX=$INDEXPATH/genome.fa
ADAPTORS=adaptors.fa
CORES=24

for FQ in $DATAPATH/P*B*R1_001.fastq.gz
do
	SAMPLENAME=$(basename $FQ)
	echo python DNAmapping.py --fq1=$FQ \
		--outdir=$PROJECTPATH --index=$INDEX \
		--threads=$CORES --adaptor=$ADAPTORS 
done 
