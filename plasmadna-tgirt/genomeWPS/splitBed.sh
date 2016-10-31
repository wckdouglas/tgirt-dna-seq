#!/bin/bash

PROJECTPATH=/scratch/02727/cdw2854/plasmaDNA
BEDPATH=${PROJECTPATH}/bedFiles/merged
SPLITBED=${BEDPATH}/splittedBed
mkdir -p ${SPLITBED}
rm -f ${SPLITBED}/*bed

for BED in ${BEDPATH}/*.bed
do
	SAMPLENAME=$(basename ${BED%.bed})
	echo awk \'\$1~/^\[0-9\]+$\|^\[XY\]$/ \{print \$0 \>\> \"${SPLITBED}/${SAMPLENAME}.\"\$1\".bed\"}\' ${BED}
done
