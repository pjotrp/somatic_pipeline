#! /bin/bash
#
# Reduce BAM files by applying design (BED). Results are written to cwd.
#
# Usage: 
#
#   reduce_bam_using_bed.sh bedfile [once-only opts] < list
#
# Create list with something like
#
#   find /data/mapping/cancer -name '*MBC*'|grep bam$
#
# e.g.
#
#   wgs01:~/data/run5/bam_reduced$ ~/opt/somatic_pipeline/scripts/reduce_bam_using_bed.sh ~/kinome_design_SS_V2_110811_nochr_annot_sorted.bed < /home/cog/pprins/current_mbc_bamlist.txt 
#

BEDTOOLS=$HOME/opt/bedtools
BEDFILE=$1
design="$(basename $BEDFILE .bed)"

BASEDIR="$( cd "$( dirname "$0" )" && pwd )"
echo $BASEDIR

while read bam ; do
  echo Reducing $bam...

  echo "$BEDTOOLS/bin/intersectBed -abam $bam -b $BEDFILE > $(basename $bam .bam).$design.bam"|~/izip/git/opensource/ruby/once-only/bin/once-only -v --pfff -d .
done



