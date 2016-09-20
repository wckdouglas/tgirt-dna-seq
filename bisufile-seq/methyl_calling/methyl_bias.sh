#!/bin/bash

PROJECT_PATH=${WORK}/cdw2854/bisufite_seq
BAM_PATH=${PROJECT_PATH}/rmdupBAM
REF_PATH=${REF}/GRCh38/hg38_rDNA
REF=${REF_PATH}/genome_rDNA.fa
RESULT_PATH=${PROJECT_PATH}/methyl_bias
mkdir -p ${RESULT_PATH}/logs


for BAM in ${BAM_PATH}/*bam
do
	SAMPLENAME=$(basename ${BAM%.bam})
	echo PileOMeth mbias ${REF} ${BAM} ${RESULT_PATH}/${SAMPLENAME} 
done
