#!/bin/bash

LINK=ftp://ftp.ensembl.org/pub/release-83/regulation/homo_sapiens/AnnotatedFeatures.gff.gz
REFPATH=/Users/wckdouglas/plasmaDNA/reference
FILENAME=openChrom.gff

curl $LINK \
	| gunzip \
	| grep open_chrom \
	| sed 's/^chr//g' \
	> $REFPATH/$FILENAME
