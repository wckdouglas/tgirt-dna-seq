#!/bin/bash

PROJECT_PATH=${WORK}/cdw2854/bisufite_seq
BAM_PATH=${PROJECT_PATH}/rmdupBAM
REF_PATH=${REF}/GRCh38/hg38_rDNA
REF=${REF_PATH}/genome_rDNA.fa
RESULT_PATH=${PROJECT_PATH}/methyl_calling
mkdir -p ${RESULT_PATH}/logs

for BAM in ${BAM_PATH}/PDB*bam
do
	SAMPLENAME=$(basename ${BAM%.bam})
	echo PileOMeth extract -q 1 -p 25 -D 10 --fraction ${REF} ${BAM} \
		-o ${RESULT_PATH}/${SAMPLENAME} \&\> ${RESULT_PATH}/logs/${SAMPLENAME}.log
done
