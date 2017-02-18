#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/ecoli_genome
BAM_PATH=$PROJECT_PATH/picard_results
OUT_PATH=$BAM_PATH/no_hiINDEL
REF_PATH=$REF/Ecoli
mkdir -p $OUT_PATH

for BAM in $BAM_PATH/*.MarkDuplicate.[0-9].subsampled.bam
do
	SAMPLENAME=$(basename ${BAM%.bam})
	echo bedtools intersect -v \
			-abam $BAM \
			-b $REF_PATH/k12_mg1655_high_indel.bed \
		\| picard CollectAlignmentSummaryMetrics \
			INPUT=/dev/stdin \
			OUTPUT= $OUT_PATH/${SAMPLENAME}.no_hiINDEL.alignment.metrics \
			IS_BISULFITE_SEQUENCED=false \
			REFERENCE_SEQUENCE=$REF_PATH/k12_mg1655.fa 
done
