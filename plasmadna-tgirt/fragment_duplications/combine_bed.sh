#!/bin/bash


PROJECT_PATH=/stor/work/Lambowitz/cdw2854/plasmaDNA
SPLITTED_BED=$PROJECT_PATH/bedFiles/umi_splitted/demultiplexed
COMBINED_BED=$PROJECT_PATH/bedFiles

SAMPLES=$(ls $SPLITTED_BED/*unique.*.bed \
	| rev \
	| cut -d'/' -f1 \
	| rev \
	| cut -d'.' -f1 \
	| uniq  )

for SAMPLENAME in $SAMPLES
do
	echo cat $SPLITTED_BED/${SAMPLENAME}.*.bed \
		\> $COMBINED_BED/${SAMPLENAME}.bed
done
