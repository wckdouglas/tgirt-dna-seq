#!/bin/bash

python make_regions.py 4
bash qc_no_hiINDEL.sh > command.sh
parallel :::: command.sh
