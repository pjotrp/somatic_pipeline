#! /bin/bash
#
# Reduce BAM files by applying design (BED). Results are written to cwd.
#
# Usage: 
#
#   reduce_bam_using_bed.sh bedfile < list
#
# Create list with something like
#
#   find /data/mapping/cancer -name '*MBC*'|grep bam$

BEDTOOLS=$HOME/opt/bedtools
BEDFILE=$1

BASEDIR="$( cd "$( dirname "$0" )" && pwd )"
echo $BASEDIR

while read bam ; do
  echo Reducing $bam...
  echo "$BEDTOOLS/bin/bed_reduce_bam.sh $BEDFILE $bam"|~/izip/git/opensource/ruby/once-only/bin/once-only -v --ignore-lock --ignore-queue -d .
  done
done



