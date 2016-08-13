#!/bin/bash

PROJECT_PATH=/scratch/02727/cdw2854/jurkatCells
BAM_PATH=${PROJECT_PATH}/bamFiles
RESULT_PATH=${PROJECT_PATH}/mismatches
REF_PATH=/corral-repl/utexas/2013lambowitz/Ref/hg19/Sequence/WholeGenomeFasta
REF=${REF_PATH}/genome.fa
SMALL_RNA_BED=${REF_PATH}/smRNA.bed
PILE_UP_BED=${WORK}/lonestar/src/pileup2bed/bin/pileup2bed
mkdir -p ${RESULT_PATH}

for BAM in ${BAM_PATH}/*bam
do
	SAMPLENAME=$(basename ${BAM%.bam})
#	if [ ! -f ${BAM}.bai ]
#	then
#	    samtools index $BAM
#	fi
#	echo python pileupBamToBase.py --bamfile=$BAM --refFasta=$REF\
#		-q 33 -d 10000 -p 12 -s 0 -o  $RESULT_PATH/${SAMPLENAME}.tsv
#	echo python basePosExtraction.py --inBam=$BAM --refFasta=$REF\
#		--qual=33 \> ${RESULT_PATH}/${SAMPLENAME}.tsv
	# bedtools intersect -v -a $BAM -b $SMALL_RNA_BED \
#	echo samtools mpileup -f $REF --max-depth 10000 -q 30 $BAM \
#		\| ${PILE_UP_BED} - 37 0  \
#		\> ${RESULT_PATH}/${SAMPLENAME}.tsv
	echo $(which python) basePosExtraction.py $BAM $CHROM
done
