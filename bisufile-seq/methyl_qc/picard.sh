#!/bin/bash

PROJECT_PATH=${WORK}/cdw2854/bisufite_seq
BAMPATH=${PROJECT_PATH}/rmdupBAM
METRIC_PATH=${PROJECT_PATH}/methyl_metric
REF_PATH=${REF}/GRCh38/hg38_rDNA
REF_FASTA=${REF_PATH}/genome_rDNA.fa
mkdir -p ${METRIC_PATH}

for BAM in ${BAMPATH}/*bam
do
	SAMPLENAME=$(basename ${BAM%.bam})
	echo picard CollectAlignmentSummaryMetrics \
		INPUT=${BAM} \
		REFERENCE_SEQUENCE=${REF_FASTA} \
		IS_BISULFITE_SEQUENCED=true \
		ASSUME_SORTED=true \
		OUTPUT=${METRIC_PATH}/${SAMPLENAME}.alignment.metrc
done
