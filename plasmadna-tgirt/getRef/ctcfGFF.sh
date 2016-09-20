#/bin/bash

LINK=ftp://ftp.ensembl.org/pub/release-83/regulation/homo_sapiens//MotifFeatures.gff.gz
REFPATH=${REF}/GRCh38/Bed_for_counts_only

curl ${LINK} | zcat | sed 's/^chr//g' > ${REFPATH}/MotifFeatures.gff
