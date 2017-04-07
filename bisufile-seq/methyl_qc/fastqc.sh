#!/bin/bash

DATAPATH=${DATA}/JA16234/bisulfite

fastqc -t 12 ${DATAPATH}/*.gz
