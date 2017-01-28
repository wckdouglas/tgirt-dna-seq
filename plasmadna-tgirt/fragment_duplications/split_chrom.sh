#!/bin/bash

PROJECTPATH=/scratch/02727/cdw2854/plasmaDNA
PROJECTPATH=$WORK/cdw2854/plasmaDNA
BEDPATH=${PROJECTPATH}/bedFiles
SPLITBED=${BEDPATH}/umi_splitted
mkdir -p ${SPLITBED}
rm -f ${SPLITBED}/*bed

#for BED in $(ls ${BEDPATH}/*umi2id*.bed | rev | cut -d'/' -f1 | rev |cut -d'-' -f2- | sort| uniq)
for BED in $(ls ${BEDPATH}/P11*umi2id*.bed )
do
	SAMPLENAME=$(basename ${BED%.bed})
#	echo cat $BEDPATH/*-${SAMPLENAME}.bed \
	echo cat $BED \
		\| awk \'\$1~/^\[0-9\]+$\|^\[XY\]$/ \{print \$0 \>\> \"${SPLITBED}/${SAMPLENAME}.\"\$1\".bed\"}\'
done
