#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/ecoli_genome/rawData/k12/umi2id_30
BED_PATH=$PROJECT_PATH/bedFiles
COUNT_PATH=$PROJECT_PATH/frag_counts
mkdir -p $COUNT_PATH

for BED in $BED_PATH/*unique.bed
do
	SAMPLENAME=$(basename ${BED%.bed})
	echo cat $BED \
		\| cut -f1,2,3,6 \
		\| sort -k1,1 -k2,2n -k3,3n -k4,4 \
		\| uniq -c \
		\| awk \''{print $2,$3,$4,$5,$1}'\' OFS=\''\t'\' \
		\| cut -f5 \
		\| sort  \
		\| uniq -c \
		\| awk \''{print $2, $1}'\' \
		\> $COUNT_PATH/${SAMPLENAME}.tsv
done
