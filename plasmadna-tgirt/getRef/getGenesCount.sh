URL1="http://www.proteinatlas.org/download/rna_tissue.csv.zip"
URL2="http://www.proteinatlas.org/download/rna_celline.csv.zip"
for URL in $URL1 $URL2
do
	FILENAME=$(basename $URL)
	OUTFILE=$WORK/cdw2854/plasmaDNA/genes/$FILENAME
	curl -o $OUTFILE $URL
done

