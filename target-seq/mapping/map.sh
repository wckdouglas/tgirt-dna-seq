#!/bin/bash

PROJECT_PATH=${WORK}/cdw2854/target-seq
SPLIT_DATA=${PROJECT_PATH}/splitted
BAM_PATH=${PROJECT_PATH}/bamFiles
TRIM_PATH=${PROJECT_PATH}/trimmed
INDEX_PATH=${REF}/targeted_gene
ADAPTOR=adaptors.fa
SUFFIX=_R1_001.fastq.gz
CORES=12
mkdir  ${BAM_PATH} ${NAME_SORT_BAM_PATH} ${TRIM_PATH}

for READ1 in ${SPLIT_DATA}/*${SUFFIX}
do
	SAMPLE_NAME=$(basename ${READ1%$SUFFIX})
	READ2=${READ1/R1/R2}
	TRIM_R1=${TRIM_PATH}/${SAMPLE_NAME}_1P.fastq.gz
	TRIM_R2=${TRIM_PATH}/${SAMPLE_NAME}_2P.fastq.gz

	trimmomatic PE \
		-threads ${CORES} \
		-basein ${READ1} \
		-baseout ${TRIM_PATH}/${SAMPLE_NAME}.fastq.gz \
		ILLUMINACLIP:adaptors.fa:2:10:10:2:true MINLEN:20 

	#MAPPING
	bwa mem -t ${CORES} \
		${INDEX_PATH}/hist1h3b.fa \
		${TRIM_R1} \
		${TRIM_R2} \
	| samtools sort \
		-@ ${CORES} \
		-O bam \
		-T ${NAME_SORT_BAM_PATH}/${SAMPLE_NAME} \
	> ${BAM_PATH}/${SAMPLE_NAME}.bam 
	
	#SORT BAM
	samtools index ${BAM_PATH}/${SAMPLE_NAME}.bam
done

python filtering.py
