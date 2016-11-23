#/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/ecoli_genome
BAM_PATH=$PROJECT_PATH/bamFiles/subsampled_bam
OUT_PATH=$PROJECT_PATH/picard_results
REF_FASTA=$REF/Ecoli/k12_mg1655.fa

for BAM in $BAM_PATH/*bam
do
	SAMPLENAME=$(basename ${BAM%.bam})
	echo picard CollectAlignmentSummaryMetrics \
		REFERENCE_SEQUENCE=$REF_FASTA \
		INPUT=$BAM \
		ASSUME_SORTED=true \
		OUTPUT=$OUT_PATH/${SAMPLENAME}.alignment_summary	
done
