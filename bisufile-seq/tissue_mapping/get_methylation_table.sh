#!/bin/bash

TARGET_PATH=/stor/work/Lambowitz/ref/hg19/methylation
LINK=http://www.pnas.org/content/suppl/2015/09/15/1508736112.DCSupplemental/pnas.1508736112.sd01.xlsx

FILENAME=$(basename $LINK)
curl -o $TARGET_PATH/$FILENAME $LINK
