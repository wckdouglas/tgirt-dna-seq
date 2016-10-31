#!/bin/bash

PROJECTPATH=$WORK/cdw2854/plasmaDNA
BAMPATH=$PROJECTPATH/bamFiles
RMDUP_BAM=$PROJECTPATH/rmdupBamFiles
SORT_BAM=$PROJECTPATH/sortedBamFiles
CORES=12
mkdir -p $RMDUP_BAM $SORT_BAM

for BAM in $BAMPATH/*bam
do
	SAMPLENAME=$(basename ${BAM%.bam})
	DEDUP_TMP_FILE=$RMDUP_BAM/$SAMPLENAME
	SORT_TMP_FILE=$RMDUP_BAM/$SAMPLENAME
	SORTED_BAM=$SORT_BAM/${SAMPLENAME}.bam
	DEDUP_BAM=$RMDUP_BAM/${SAMPLENAME}.bam
	echo sambamba sort \
			--nthreads=$CORES \
			--tmpdir=$SORT_TMP_FILE \
			--out=$SORTED_BAM \
			--show-progress \
			$BAM \;\
		sambamba markdup \
			--nthreads=$CORES \
			--tmpdir=$DEDUP_TMP_FILE \
			--show-progress \
			$SORTED_BAM \
			$DEDUP_BAM
done
