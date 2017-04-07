#!/bin/bash

PROJECTPATH=/scratch/02727/cdw2854/plasmaDNA
PROJECTPATH=$WORK/cdw2854/plasmaDNA
BEDPATH=${PROJECTPATH}/bedFiles
SPLITBED=${BEDPATH}/splittedBed
mkdir -p ${SPLITBED}
#rm -f ${SPLITBED}/*bed

for BED in ${BEDPATH}/P1022_1113*.bed
do
	SAMPLENAME=$(basename ${BED%.bed})
	echo awk \'\$1~/^\[0-9\]+$\|^\[XY\]$/ \{print \$0 \>\> \"${SPLITBED}/${SAMPLENAME}.\"\$1\".bed\"}\' ${BED}
done
