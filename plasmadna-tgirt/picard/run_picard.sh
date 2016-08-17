#!/bin/bash

PROJECT_PATH=${WORK}/cdw2854/plasmaDNA
BAM_PATH=${PROJECT_PATH}/bamFiles
REF_PATH=${REF}/GRCh38/hg38_rDNA
REF=${REF_PATH}/genome_rDNA.fa


for BAM in ${BAM_PATH}/*bam
do
    SAMPLENAME=$(basename ${BAM%.bam})
	if [[ $SAMPLENAME != G* ]]
	then
		echo bash picard_gc.sh ${BAM} ${REF}
	fi
done
