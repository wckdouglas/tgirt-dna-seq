#!/bin/bash

DATAPATH=${SCRATCH}/plasmaDNA/rawData
SPLIT_PATH=${DATAPATH}/split
PROGRAM=/work/02727/cdw2854/lonestar/src/fastq-tools/bin/splitFastq 

for FQ in $DATAPATH/P1* $DATAPATH/SRR*52_R[12]_001.fq.gz
do
	SAMPLENAME=$(basename ${FQ%.fastq.gz})
	echo ${PROGRAM} -i ${FQ} -o ${SPLIT_PATH}/${SAMPLENAME} -z -n 5000000
done
