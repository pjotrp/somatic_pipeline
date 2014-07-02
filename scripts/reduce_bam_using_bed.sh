#! /bin/bash
#
# Reduce BAM files by applying design (BED). Results are written to cwd.
#
# Usage: 
#
#   reduce_bam_using_bed.sh bedfile [once-only opts] < list
#
# Create list with something like (removing files with wf in the name)
#
#   find /data/mapping/cancer -name '*MBC*'|grep -v wf|grep bam$
#
# e.g.
#
#   wgs01:~/data/run5/bam_reduced$ ~/opt/somatic_pipeline/scripts/reduce_bam_using_bed.sh ~/kinome_design_SS_V2_110811_nochr_annot_sorted.bed < /home/cog/pprins/current_mbc_bamlist.txt 
#
# Optionally a basedir to the data files can be passed in.

BEDTOOLS=$HOME/opt/bedtools
onceonly="once-only"
BEDFILE=$1
DATAPATH=$2
design="$(basename $BEDFILE .bed)"

BASEDIR="$( cd "$( dirname "$0" )" && pwd )"
echo $BASEDIR

while read bam ; do
  bam="$DATAPATH/$bam"
  echo Reducing $bam...
  outfn=$(basename $bam .bam).$design.bam
  echo "$BEDTOOLS/bin/intersectBed -abam $bam -b $BEDFILE|cat > $outfn"|$onceonly -v --pbs '-q veryshort' --pfff --out $outfn -d .
  # echo "$BEDTOOLS/bin/intersectBed -abam $bam -b $BEDFILE|cat > $outfn"|~/izip/git/opensource/ruby/once-only/bin/once-only -v --pfff --out $outfn -d .
done



