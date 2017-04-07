#!/bin/bash

PROJECTPATH=$WORK/cdw2854/plasmaDNA
#PROJECTPATH=$WORK/cdw2854/plasmaDNA
BEDPATH=${PROJECTPATH}/bedFiles
SPLITBED=${BEDPATH}/umi_splitted
mkdir -p ${SPLITBED}
#rm -f ${SPLITBED}/*bed

#for BED in $(ls ${BEDPATH}/*umi2id*.bed | rev | cut -d'/' -f1 | rev |cut -d'-' -f2- | sort| uniq)
<<<<<<< HEAD
for BED in $(ls ${BEDPATH}/merged/P1*umi2id*.bed )
=======
for BED in $(ls ${BEDPATH}/P10*umi2id*.bed )
>>>>>>> 57d5eccb6c1ad751fd8b63eae7888be42a69d8cf
do
	SAMPLENAME=$(basename ${BED%.bed})
#	echo cat $BEDPATH/*-${SAMPLENAME}.bed \
	echo cat $BED \
		\| awk \'\$1~/^\[0-9\]+$\|^\[XY\]$/ \{print \$0 \>\> \"${SPLITBED}/${SAMPLENAME}.\"\$1\".bed\"}\'
done
