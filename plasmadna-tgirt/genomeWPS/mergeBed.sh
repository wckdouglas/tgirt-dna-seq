#!/bin/bash

PROJECT_PATH=/scratch/02727/cdw2854/plasmaDNA
BED_PATH=$PROJECT_PATH/genomeWPS
MERGED_PATH=$BED_PATH/mergedStrandedBED
mkdir -p $MERGED_PATH

SAMPLES=$(ls $BED_PATH/*[de].bed | awk -F'/' {'print $NF'} | cut -d'.' -f1 | sort | uniq)
for SAMPLENAME in $SAMPLES
do
	for LENGTH in Long Short
	do
			cat $BED_PATH/$SAMPLENAME*.$LENGTH.*.bed > $MERGED_PATH/$SAMPLENAME.$LENGTH.bed
	done
done
