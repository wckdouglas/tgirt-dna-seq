#!/bin/bash

PROJECTPATH=/scratch/02727/cdw2854/plasmaDNA
<<<<<<< HEAD
BEDPATH=${PROJECTPATH}/bedFiles/umi_splitted/demultiplexed/merged
=======
PROJECTPATH=$WORK/cdw2854/plasmaDNA
BEDPATH=${PROJECTPATH}/bedFiles
>>>>>>> 57d5eccb6c1ad751fd8b63eae7888be42a69d8cf
SPLITBED=${BEDPATH}/splittedBed
mkdir -p ${SPLITBED}
#rm -f ${SPLITBED}/*bed

<<<<<<< HEAD
for BED in ${BEDPATH}/*unique.bed
=======
for BED in ${BEDPATH}/P1022_1113*.bed
>>>>>>> 57d5eccb6c1ad751fd8b63eae7888be42a69d8cf
do
	SAMPLENAME=$(basename ${BED%.bed})
	echo awk \'\$1~/^\[0-9\]+$\|^\[XY\]$/ \{print \$0 \>\> \"${SPLITBED}/${SAMPLENAME}.\"\$1\".bed\"}\' ${BED}
done
