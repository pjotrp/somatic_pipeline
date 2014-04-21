#! /bin/bash
# 
# Simple rmdup loop
#
#   ./rmdup.sh < bamlist

while read bam ; do
  echo ==== Remove duplicates $bam...
  outfn=$(basename $bam .bam)_rmdup.bam
  echo "$sambamba markdup -r $bam $outfn"| $onceonly --pfff -d . -v --out $outfn
  [ $? -ne 0 ] && exit 1
done

