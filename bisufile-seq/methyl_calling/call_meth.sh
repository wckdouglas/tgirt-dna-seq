#!/bin/bash

PROJECT_PATH=/scratch/02727/cdw2854/plasmaDNA
BAM_PATH=${PROJECT_PATH}/bamFiles
REF_PATH=/corral-repl/utexas/2013lambowitz/Ref/GRCh38/hg38_rDNA
REF=${REF_PATH}/genome_rDNA.fa
RESULT_PATH=${PROJECT_PATH}/methyl_calling
mkdir -p ${RESULT_PATH}/logs

for BAM in ${BAM_PATH}/PDB*bam
do
	SAMPLENAME=$(basename ${BAM%.bam})
	echo PileOMeth extract -q 1 -p 25 -D 10 --fraction ${REF} ${BAM} \
		-o ${RESULT_PATH}/${SAMPLENAME} \&\> ${RESULT_PATH}/logs/${SAMPLENAME}.log
done
