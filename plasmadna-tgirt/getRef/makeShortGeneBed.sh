#!/bin/bash

BEDPATH=/Users/wckdouglas/plasmaDNA/reference
BEDNAME=genes_count.bed

awk '$7~/(tRNA|snoRNA|YRNA|miRNA|piRNA|snRNA)/' \
	$BEDPATH/$BEDNAME \
| sort -k1,1 -k2,2n \
> $BEDPATH/smallRNA.bed  

