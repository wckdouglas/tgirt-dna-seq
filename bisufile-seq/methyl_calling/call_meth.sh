#!/bin/bash

PROJECT_PATH=${WORK}/cdw2854/bisufite_seq
BAM_PATH=${PROJECT_PATH}/rmdupBAM
BAM_PATH=${PROJECT_PATH}/bamFiles
REF_PATH=${REF}/hg19/Sequence/WholeGenomeFasta
REF=${REF_PATH}/genome.fa
RESULT_PATH=${PROJECT_PATH}/methyl_calling
mkdir -p ${RESULT_PATH}/logs

for BAM in ${BAM_PATH}/P13*bam
do
	SAMPLENAME=$(basename ${BAM%.bam})
	echo MethylDackel extract \
		-q 1 -p 25 \
		-d 2 -D 2000 \
		--fraction ${REF} ${BAM} \
		-o ${RESULT_PATH}/${SAMPLENAME} \&\> ${RESULT_PATH}/logs/${SAMPLENAME}.log
done
