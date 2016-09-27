#!/bin/bash

PROJECT_PATH=${WORK}/cdw2854/target-seq
SPLIT_DATA=${PROJECT_PATH}/splitted
BAM_PATH=${PROJECT_PATH}/bamFiles
MERGE_PATH=${PROJECT_PATH}/merged_fastq
INDEX_PATH=${REF}/targeted_gene
ADAPTOR=adaptors.fa
SUFFIX=_R1_001.fastq.gz
CORES=12
mkdir -p ${BAM_PATH} ${NAME_SORT_BAM_PATH} ${TRIM_PATH} ${MERGE_PATH}

for READ1 in ${SPLIT_DATA}/*${SUFFIX}
do
	SAMPLE_NAME=$(basename ${READ1%$SUFFIX})
	READ2=${READ1/R1/R2}

	flash --to-stdout ${READ1} ${READ2} \
		--threads ${CORES} \
	| tee ${MERGE_PATH}/${SAMPLE_NAME}.fq \
	| bwa mem -t ${CORES} \
		${INDEX_PATH}/histone_protein.fa \
		- \
	| samtools sort \
		-@ ${CORES} \
		-O bam \
		-T ${NAME_SORT_BAM_PATH}/${SAMPLE_NAME} \
	> ${BAM_PATH}/${SAMPLE_NAME}.bam 
	
	#SORT BAM
	samtools index ${BAM_PATH}/${SAMPLE_NAME}.bam
done
python filtering.py
#		${INDEX_PATH}/hist1h3b.fa \
#		${INDEX_PATH}/hist1h3.fa \
bash ../mismatch/extract_base.sh
