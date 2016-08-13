#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/plasmaDNA
BAM_PATH=$PROJECT_PATH/bamFiles
MERGED_BAM=$BAM_PATH/mergedBam
BEDPATH=$PROJECT_PATH/bedFiles
PYTHON=$(which python)
BAMTOOLS=$(which bamtools)
SAMTOOLS=$(which samtools)
PICARD=$(which picard)
SAMBAMBA=$(which sambamba)
THREADS=12
mkdir -p ${BEDPATH} ${MERGED_BAM}

for SAMPLENAME in `ls $BAM_PATH/[NPS]*bam | rev | cut -d'/' -f1 | rev | cut -d'.' -f1 | cut -d'_' -f1 | sort | uniq`
do
	TEMP_DIR=$BEDPATH/${SAMPLENAME}_tmp
	mkdir -p $TEMP_DIR
	INBAMS=$BAM_PATH/${SAMPLENAME}_*bam
	echo ${SAMTOOLS} cat ${INBAMS} \
		\| ${SAMTOOLS} fixmate -O BAM -c -r  - - \
		\| tee ${MERGED_BAM}/${SAMPLENAME}.nameSorted.bam \
		\| ${PYTHON} bamtobed.py /dev/stdin \
		\| sort -k1,1 -k2,2n -k3,3n -k6,6  --temporary-directory=${TEMP_DIR} -u \
		\> ${BEDPATH}/${SAMPLENAME}.bed  
done
