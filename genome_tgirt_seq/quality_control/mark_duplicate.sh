#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/ecoli_genome
BAM_PATH=$PROJECT_PATH/bamFiles
MARK_DUP_PATH=$PROJECT_PATH/bamFiles/mark_duplicate
mkdir -p $MARK_DUP_PATH

for BAM in $BAM_PATH/*bam
do
	BAM_NAME=$(basename $BAM)
	echo picard MarkDuplicates \
		INPUT=$BAM \
		OUTPUT=${MARK_DUP_PATH}/${BAM_NAME} \
		METRICS_FILE=${MARK_DUP_PATH}/${BAM_NAME%.bam}.marked_dup_metrics.txt \
		ASSUME_SORT_ORDER=coordinate 
done

