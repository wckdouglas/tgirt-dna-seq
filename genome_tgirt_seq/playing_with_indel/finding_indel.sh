#!/bin/bash

PROJECT_PATH=$WORK/cdw2854/ecoli_genome
BAM_PATH=$PROJECT_PATH/picard_results
RESULT_PATH=$PROJECT_PATH/indel_table
REF_FASTA=$REF/Ecoli/k12_mg1655.fa
SUFFIX=.MarkDuplicate.*.subsampled.bam
mkdir -p $RESULT_PATH

for BAM in $( ls ${BAM_PATH}/*${SUFFIX} | grep -v 'sim')
do
	SAMPLENAME=$(basename ${BAM%.subsampled.bam})
	echo samtools view \
		-b \
		-F 1024 \
		-F 2048 \
		-F 256  $BAM \
	\| $HOME/src/piledriver/bin/bamtools piledriver \
		-fasta $REF_FASTA \
	\| cut -f1-13 \
	\| awk \''$1 == "chrom" || $(NF-1)>0 || $NF >0 '\' \
	\> $RESULT_PATH/${SAMPLENAME}.tsv
done
