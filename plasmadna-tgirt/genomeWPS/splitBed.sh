#!/bin/bash

PROJECTPATH=/scratch/02727/cdw2854/plasmaDNA
BEDPATH=$PROJECTPATH/rmdup_bed
SPLITBED=$BEDPATH/splittedBed
mkdir -p $SPLITBED

for BED in $BEDPATH/*.bed
do
	SAMPLENAME=$(basename ${BED%.bed})
	awk \'\$1~/^\[1-9\]+$\|^\[XY\]$/ \{print \$0 \>\> \"$SPLITBED/$SAMPLENAME.\"\$1\".bed\"}\' $BED
done
