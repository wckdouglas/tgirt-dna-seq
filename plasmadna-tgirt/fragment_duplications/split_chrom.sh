#!/bin/bash

PROJECTPATH=/scratch/02727/cdw2854/plasmaDNA
BEDPATH=${PROJECTPATH}/bedFiles
SPLITBED=${BEDPATH}/umi_splitted
mkdir -p ${SPLITBED}
#rm -f ${SPLITBED}/*bed

for BED in $(ls ${BEDPATH}/*.bed | rev | cut -d'/' -f1 | rev |cut -d'-' -f2- | sort| uniq)
do
	SAMPLENAME=${BED%.bed}
	echo cat $BEDPATH/*-${SAMPLENAME}.bed \
		\| awk \'\$1~/^\[0-9\]+$\|^\[XY\]$/ \{print \$0 \>\> \"${SPLITBED}/${SAMPLENAME}.\"\$1\".bed\"}\'
done
