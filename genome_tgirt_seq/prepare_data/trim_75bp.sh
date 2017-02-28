#!/bin/bash

DATAPATH=/stor/work/Lambowitz/cdw2854/ecoli_genome/rawData/k12
OUT_PATH=$DATAPATH/CUT_75BP
mkdir -p $OUT_PATH

for FQ in $DATAPATH/*fastq.gz
do
	SAMPLENAME=$(basename $FQ)
	echo seqtk trimfq -L 75 $FQ \
		\| gzip \
		\> $OUT_PATH/$SAMPLENAME
done
