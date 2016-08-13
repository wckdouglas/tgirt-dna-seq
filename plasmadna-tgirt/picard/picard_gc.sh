#!/bin/bash

PROJECT_PATH=${WORK}/cdw2854/plasmaDNA
BAM_PATH=${PROJECT_PATH}/bamFiles
SORTED_BAM_PATH=${BAM_PATH}/sorted
RESULT_PATH=${PROJECT_PATH}/picard_results
FIGURES_PATH=${PROJECT_PATH}/figures
REF_PATH=${REF}/GRCh38/hg38_rDNA
REF=${REF_PATH}/genome_rDNA.fa
mkdir -p ${RESULT_PATH} ${FIGURES_PATH} ${SORTED_BAM_PATH}

for BAM in ${BAM_PATH}/*bam
do
    SAMPLENAME=$(basename ${BAM%.bam})
	if [[ $SAMPLENAME != G* ]]
	then
		echo Running $SAMPLENAME
		picard SortSam \
			INPUT=${BAM} \
			OUTPUT=/dev/stdout 
			SORT_ORDER=coordinate \
			> ${SORTED_BAM_PATH}/${SAMPLENAME}.sorted.bam 
		echo Sorted $SAMPLENAME
		picard MarkDuplicates \
			INPUT=/dev/stdin \
			OUTPUT=/dev/stdout \
			METRICS_FILE=${RESULT_PATH}/${SAMPLENAME}.duplicate.txt \
			REMOVE_DUPLICATES=True \
			ASSUME_SORTED=true \
			> ${SORTED_BAM_PATH}/${SAMPLENAME}.sorted.rmdup.bam 
		echo Removed duplicates $SAMPLENAME
		echo Collecting Bias
		picard CollectGcBiasMetrics \
			VALIDATION_STRINGENCY=LENIENT \
			SCAN_WINDOW_SIZE=100 \
			INPUT=/dev/stdin \
			CHART_OUTPUT=${FIGURES_PATH}/${SAMPLENAME}.pdf  \
			OUTPUT=${RESULT_PATH}/${SAMPLENAME}.txt \
			SUMMARY_OUTPUT=${RESULT_PATH}/${SAMPLENAME}.GC.summary  \
			REFERENCE_SEQUENCE=$REF \
			ASSUME_SORTED=true 
	fi
done
