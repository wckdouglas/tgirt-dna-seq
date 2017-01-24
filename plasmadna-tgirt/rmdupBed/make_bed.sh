#!/bin/bash

PROJECT_PATH=/stor/work/Lambowitz/cdw2854/plasmaDNA/
BAM_PATH=$PROJECT_PATH/rmdupBamFiles
BED_PATH=$PROJECT_PATH/rmdupBedFiles
CORES=12
mkdir -p $BED_PATH

for BAM in $BAM_PATH/*.bam
do
	SAMPLENAME=$(basename ${BAM%.bam})
	echo sambmba sort \
			--nthreads=$CORES \
			--show-progress \
			--tmpdir=$BED_PATH/$SAMPLENAME \
			--sort-by-name \
			--out=/dev/stdout \
			$BAM \
		\| samtools view -F 1024 -F 2048 -b@ $CORES \
		\| bedtools bamtobed -i - -bedpe -mate1 \
		\| awk \''$1!="."'\' \
		\| bedpe_to_bed.py -i - -o - --min 10 --max 100000\
		\> $BED_PATH/${SAMPLENAME}.bed
done
