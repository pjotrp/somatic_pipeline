#! /bin/bash
# 
# Simple rmdup loop
#
#   sambamba_rmdup.sh < bamlist
#
# Example
#
#   ls -1 --color=never *.bam | rmdup.sh

sambamba=sambamba
onceonly=$HOME/izip/git/opensource/ruby/once-only/bin/once-only

while read bam ; do
  echo ==== $sambamba remove duplicates $bam...
  outfn=$(basename $bam .bam)_rmdup.bam
  echo "$sambamba markdup -r $bam $outfn"| $onceonly --pfff --pbs '-q veryshort' -d . -v -in $bam --out $outfn
  # echo "$sambamba markdup -r $bam $outfn"| $onceonly --pfff -d . -v -in $bam --out $outfn
  [ $? -ne 0 ] && exit 1
done

