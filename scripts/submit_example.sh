#!/bin/sh
# This is a simple example of a SGE batch script
#
# Run with
#
#   echo "/home/cog/pprins/opt/somatic_pipeline/scripts/submit.sh" | qsub -P SAP42 -cwd

# request Bourne shell as shell for job
#$ -S /bin/bash

# Execute from the current working directory
#$ -cwd

# send output to stdout
#$ -o stdout

echo "Got $NSLOTS slots."

PATH=$SGE_O_PATH:$PATH

# print date and time
date
# Sleep for 20 seconds
hostname
set
pwd
ruby -v
echo $PATH
sleep 2
# print date and time again
date
