#!/bin/bash

sh runPipeline.sh | grep -v 'try\|50\|93' > command.sh
python slurmJob.py -c command.sh -t 6:00:00 -j pipeline -N 1 -n 1
