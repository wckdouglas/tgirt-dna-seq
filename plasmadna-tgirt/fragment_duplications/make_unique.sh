#!/bin/bash

PROJECT_PATH=$WORK/cdw2854/plasmaDNA
#PROJECT_PATH=$SCRATCH/plasmaDNA
BED_PATH=$PROJECT_PATH/bedFiles/chrom_split_bed
OUT_PATH=$BED_PATH/demultiplexed
PROGRAM_PATH=$HOME/TGIRT_UMI/dedup_bed
mkdir -p $OUT_PATH

for BED in $BED_PATH/*.bed
do
	SAMLENAME=$(echo $(basename $BED) | cut -d'.' -f1)
	CHROM=$(echo $(basename $BED) | cut -d'.' -f2)
	TMP_DIR=$OUT_PATH/$SAMLENAME
	mkdir -p $TMP_DIR
	echo cat $BED \
		\| tr \'_\' \''\t'\' \
		\| awk \''{print $1,$2,$3,$4,$(NF-1),$(NF)}'\' OFS=\''\t'\' \
		\| sort -k1,1 -k2,2n -k3,3n -T $TMP_DIR \
		\| datamash -g 1,2,3,6 collapse 4 \
		\| tee $OUT_PATH/${SAMLENAME}_collapse.${CHROM}.bed \
		\| $(which python) ${PROGRAM_PATH}/unique_bed.py --infile=- --threshold=1 \
		\> $OUT_PATH/${SAMLENAME}_unique.${CHROM}.bed
done


