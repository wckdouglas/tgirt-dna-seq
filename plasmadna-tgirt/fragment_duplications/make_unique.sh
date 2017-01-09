#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/plasmaDNA
BED_PATH=$PROJECT_PATH/bedFiles
OUT_PATH=$BED_PATH/unique_bed_cluster
mkdir -p $OUT_PATH

for BED in $BED_PATH/*clustered.bed
do
	SAMLENAME=$(basename ${BED%.bed})
	echo cat $BED \
		\| tr \'_\' \''\t'\' \
		\| awk \''{print $1,$2,$3,$4,$7,$8}'\' OFS=\''\t'\' \
		\| sort -k1,1 -k2,2n -k3,3n \
		\| bedtools merge -i - -s -c 4 -o collapse -delim \",\" \
		\| tee $OUT_PATH/${SAMLENAME}.collapse.bed \
		\| python unique_bed.py \
		\> $OUT_PATH/${SAMLENAME}_unique.bed
done


