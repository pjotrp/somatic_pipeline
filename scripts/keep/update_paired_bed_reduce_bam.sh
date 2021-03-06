#! /bin/bash

BED="$HOME/kinome_design_SS_V2_110811_nochr_annot.bed"

BASEDIR="$( cd "$( dirname "$0" )" && pwd )"
echo $BASEDIR

IFS=$'\t'
while read ref tumor ; do
  echo "tumor=$tumor ref=$ref"
  for x in $tumor $ref ; do 
    echo Reducing $x...
    echo "$BASEDIR/../runners/bedtools/bin/bed_reduce_bam.sh $HOME/BED $x"|~/izip/git/opensource/ruby/once-only/bin/once-only -v --ignore-lock --ignore-queue -d . 
# --pbs "-P SAP42" 
  done
done < somatic_bams.txt



