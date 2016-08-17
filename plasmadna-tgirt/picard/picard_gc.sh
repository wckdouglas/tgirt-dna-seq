#!/bin/bash

# get variable
BAM=$1
PROJECT_PATH=$(dirname $BAM)
REF=$2

#make directory
SAMPLENAME=$(basename ${BAM%.bam})
BAM_PATH=${PROJECT_PATH}/bamFiles
SORTED_BAM_PATH=${BAM_PATH}/sorted
RESULT_PATH=${PROJECT_PATH}/picard_results
FIGURES_PATH=${PROJECT_PATH}/figures
mkdir -p ${RESULT_PATH} ${FIGURES_PATH} ${SORTED_BAM_PATH}


echo Using project directory: $PROJECT_PATH
echo Running $SAMPLENAME
#
SORTED_BAM=${SORTED_BAM_PATH}/${SAMPLENAME}.sorted.bam 
picard SortSam \
	INPUT=${BAM} \
	OUTPUT=$SORTED_BAM  \
	SORT_ORDER=coordinate 
echo Sorted $SAMPLENAME

RMDUP_BAM=${SORTED_BAM_PATH}/${SAMPLENAME}.sorted.rmdup.bam
picard MarkDuplicates \
	INPUT=$SORTED_BAM \
	OUTPUT=$RMDUP_BAM \
	METRICS_FILE=${RESULT_PATH}/${SAMPLENAME}.duplicate.txt \
	REMOVE_DUPLICATES=False\
	ASSUME_SORTED=true 
echo Marked duplicates $SAMPLENAME

echo Collecting Bias
picard CollectGcBiasMetrics \
	VALIDATION_STRINGENCY=LENIENT \
	SCAN_WINDOW_SIZE=100 \
	INPUT=$BAM \
	CHART_OUTPUT=${FIGURES_PATH}/${SAMPLENAME}.pdf  \
	OUTPUT=${RESULT_PATH}/${SAMPLENAME}.txt \
	SUMMARY_OUTPUT=${RESULT_PATH}/${SAMPLENAME}.GC.summary  \
	REFERENCE_SEQUENCE=$REF \
	ASSUME_SORTED=true 
	INPUT=$RMDUP_BAM \
echo Finished $SAMPLENAME
