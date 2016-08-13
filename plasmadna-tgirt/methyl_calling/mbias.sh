#!/bin/bash

PROJECT_PATH=/scratch/02727/cdw2854/plasma_project
BAM_PATH=${PROJECT_PATH}/bisulfite
REF_PATH=/corral-repl/utexas/2013lambowitz/Ref/GRCh38/hg38_rDNA
REF=${REF_PATH}/genome_rDNA.fa
RESULT_PATH=${PROJECT_PATH}/methyl_calling
mkdir -p ${RESULT_PATH}/logs


for BAM in ${BAM_PATH}/*bam
do
	SAMPLENAME=$(basename ${BAM%.bam})
	echo PileOMeth mbias ${REF} ${BAM} ${RESULT_PATH}/${SAMPLENAME} 
done
