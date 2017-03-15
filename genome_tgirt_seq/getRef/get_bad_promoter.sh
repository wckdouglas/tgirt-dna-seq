#!/usr/bin/bash

LINK=https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2013-14-5-r51/MediaObjects/13059_2012_3110_MOESM2_ESM.TXT
DESTINATION=$REF/hg19/annotations
curl $LINK \
	| sed 's/^/chr/g' \
	| sed 's/:/\t/g'\
	| sed 's/-/\t/g' \
	| sed 's/,//g' \
	| sed 's/#//g' \
	| sed 's/ratio=//g' \
	| awk '{print $1,$2,$3,$4,"0",".",$5}' OFS='\t' \
	| bedtools sort \
	> $DESTINATION/bad_promoters_100bp_extension.bed
	
#| bedtools slop -b 100 -i - -g $REF/hg19/Sequence/WholeGenomeFasta/genome.fa.fai \
