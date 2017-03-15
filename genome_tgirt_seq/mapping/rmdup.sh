#!/bin/bash

BAM_PATH=$WORK/cdw2854/genomeDNA/clustered_map

for BAM in $BAM_PATH/*.bam
do
	SAMPLENAME=$(basename ${BAM})
	echo sambamba markdup \
		-p -t 20 -r \
		--tmpdir=$BAM_PATH/rmdup_bam/temp \
		$BAM \
		$BAM_PATH/rmdup_bam/${SAMPLENAME}
done
