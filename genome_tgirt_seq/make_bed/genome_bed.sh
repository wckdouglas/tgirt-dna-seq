#!/usr/bin

PROJECT_PATH=/scratch/02727/cdw2854/jurkatCells
BAM_PATH=${PROJECT_PATH}/bamFiles
BED_PATH=${PROJECT_PATH}/bedFiles
BED_TO_BEDPE=/work/02727/cdw2854/lonestar/src/bedFileTools/bin/bedpeTobed
mkdir -p ${BED_PATH}

for BAM in ${BAM_PATH}/*bam
do
    SAMPLENAME=$(basename ${BAM%.bam})
    echo samtools view -b -F 2048 $BAM \
	\| bamtools filter -script filterCriteria.json \
        \| samtools sort -n -T ${BED_PATH}/${SAMPLENAME} -O bam \
        \| samtools fixmate -r - - \
        \| bedtools bamtobed -i - -bedpe -mate1 \
        \| ${BED_TO_BEDPE} -i - -m 100000 -l 10 \
	\| sort -k1,1 -k2,2n -k6,6 \
        \> ${BED_PATH}/${SAMPLENAME}.bed
done
