#!/bin/bash

PROJECT_PATH=/scratch/02727/cdw2854/plasmaDNA
BED_PATH=$PROJECT_PATH/genomeWPS
MERGED_PATH=$BED_PATH/mergedWPS
mkdir -p $MERGED_PATH

SAMPLES=$(ls $BED_PATH/*.bed | awk -F'/' {'print $NF'} | cut -d'.' -f1 | sort | uniq)
for SAMPLENAME in $SAMPLES
do
	for LENGTH in Long Short
	do
		WGS_FILES=$(ls ${BED_PATH}/${SAMPLENAME}.*.${LENGTH}.bed)
		echo cat ${WGS_FILES} \
			\> ${MERGED_PATH}/${SAMPLENAME}.${LENGTH}.bed
	done
done
