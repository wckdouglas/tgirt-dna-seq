#/bin/bash

LINK=ftp://ftp.ensembl.org/pub/release-83/regulation/homo_sapiens//MotifFeatures.gff.gz
<<<<<<< HEAD
REFPATH=/scratch/cdw2854/reference/hg38
=======
REFPATH=/corral-repl/utexas/2013lambowitz/Ref/GRCh38/Bed_for_counts_only
>>>>>>> 0b7a2f8b3cca06c6ab23ffd218ca28253fa83e6d

curl $LINK | zcat | sed 's/^chr//g' > $REFPATH/MotifFeatures.gff
