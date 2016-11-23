#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/target-seq
FILTERED_BAM_PATH=$PROJECT_PATH/filtered_bams
BAM_PATH=$PROJECT_PATH/bamFiles
BASE_PATH=$PROJECT_PATH/base_tables
SCRIPT_PATH=${HOME}/tgirt-dna-seq/target-seq/mismatch
INDEX=$REF/targeted_gene/hist1h3b.fa
INDEX=$REF/targeted_gene/histone_protein.fa
mkdir -p $BASE_PATH

for BAM in $FILTERED_BAM_PATH/*bam
do
	SAMPLENAME=$(basename ${BAM%.bam})
	OUT_BAM=$BASE_PATH/${SAMPLENAME}.bam
	echo python ${SCRIPT_PATH}/filter_bam.py ${BAM} ${OUT_BAM}';' \
	samtools mpileup \
		--max-depth 10000000 \
		-f $INDEX \
		${OUT_BAM} \
	\| pileup_to_bed.py \
		--input=- \
		--qual=40 \
		--threads 4 \
	\> $BASE_PATH/${SAMPLENAME}.tsv
done > command.sh
echo Running: 
cat command.sh
parallel :::: command.sh
#	\| bamtools filter -tag '"NM:<4"' -region ENST00000621411  \
