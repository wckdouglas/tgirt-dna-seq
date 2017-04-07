#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/ecoli_genome/rawData/k12/umi2id_30
BAM_PATH=$PROJECT_PATH/bamFiles
RESULTPATH=$PROJECT_PATH/bedFiles
mkdir -p $RESULTPATH

for BAM in $BAM_PATH/*bam
do
	SAMPLENAME=$(basename ${BAM%.bam})
	echo samtools sort -n $BAM \
		\| bedtools bamtobed -bedpe -mate1 \
		\| awk \''$1!="."'\' \
		\| bedpe_to_bed.py \
		\| tr \'_\' \'\\t\' \
		\| awk \''{printf "%s_%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,$5,$7,$8,$9 }'\' OFS=\'\\t\' \
		\| sort -k1,1 -k2,2n -k3,3n \
		\| datamash -g 1,2,3,6 collapse 4 \
		\| tee $RESULTPATH/${SAMPLENAME}.collapse.bed \
		\| python ~/TGIRT_UMI/dedup_bed/unique_bed.py -t 1 -i - \
		\> $RESULTPATH/${SAMPLENAME}.unique.bed
done

