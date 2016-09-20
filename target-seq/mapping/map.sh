#!/bin/bash

PROJECT_PATH=${WORK}/cdw2854/target-seq
SPLIT_DATA=${PROJECT_PATH}/splitted
SPLIT_DATA=/stor/work/Lambowitz/Data/NGS/JA16493/splitted
BAM_PATH=${PROJECT_PATH}/bamFiles
TRIM_PATH=${PROJECT_PATH}/trimmed
INDEX_PATH=${REF}/targeted_gene
INDEX_PATH=${REF}/human_transcriptome
NAME_SORT_BAM_PATH=${PROJECT_PATH}/name_sorted
ADAPTOR=adaptors.fa
SUFFIX=_R1_001.fastq.gz
CORES=12
mkdir  ${BAM_PATH} ${NAME_SORT_BAM_PATH} ${TRIM_PATH}

for READ1 in ${SPLIT_DATA}/*${SUFFIX}
do
	SAMPLE_NAME=$(basename ${READ1%$SUFFIX})
	TRIM_FILE=${TRIM_PATH}/${SAMPLE_NAME}
	#TRIMMING
	trimmomatic PE -phred64 -threads ${CORES} \
		-basein ${READ1} -baseout ${TRIM_FILE}.fq.gz \
		ILLUMINACLIP:${ADAPTOR}:2:10:10:1:true MINLEN:15

	#MAPPING
	bowtie2 -p ${CORES} --local \
		-x ${INDEX_PATH}/transcriptome.fa \
		-1 ${TRIM_FILE}_1P.fq.gz \
		-2 ${TRIM_FILE}_2P.fq.gz \
		> ${NAME_SORT_BAM_PATH}/${SAMPLE_NAME}.bam 
	
	#SORT BAM
	samtools sort -@ ${CORES} -O bam -T ${NAME_SORT_BAM_PATH}/${SAMPLE_NAME} ${NAME_SORT_BAM_PATH}/${SAMPLE_NAME}.bam \
		> ${BAM_PATH}/${SAMPLE_NAME}.bam 
	samtools index ${BAM_PATH}/${SAMPLE_NAME}.bam
done
