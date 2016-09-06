#!/bin/bash

PROJECT_PATH=${WORK}/cdw2854/target-seq
SPLIT_DATA=${PROJECT_PATH}/splitted
SPLIT_DATA=/stor/work/Lambowitz/Data/NGS/JA16493/splitted
BAM_PATH=${PROJECT_PATH}/bamFiles
SUFFIX=_R1_001.fastq.gz
INDEX_PATH=${REF}/targeted_gene
INDEX_PATH=${REF}/human_transcriptome
mkdir  ${BAM_PATH}

for READ1 in ${SPLIT_DATA}/*${SUFFIX}
do
	READ2=${READ1/R1/R2}
	SAMPLE_NAME=$(basename ${READ1%$SUFFIX})
	echo bowtie2 \
		-p 12 \
		--local \
		-x ${INDEX_PATH}/transcriptome.fa \
		-1 ${READ1} \
		-2 ${READ2} \
		\|samtools sort -@ 12 -O bam -T ${BAM_PATH}/${SAMPLE_NAME}\
		\> ${BAM_PATH}/${SAMPLE_NAME}.bam
done
