#! /bin/bash
# 
# Simple sambamba index loop
#
#   ./index.sh < bamlist
#
# Example
#
#   ls -1 --color=never *_rmdup.bam | index.sh

sambamba=$HOME/opt/bin/sambamba
onceonly=$HOME/izip/git/opensource/ruby/once-only/bin/once-only

while read bam ; do
  echo ==== $sambamba index $bam...
  outfn=$(basename $bam .bam).bam.bai
  echo "$sambamba index $bam"| $onceonly --pfff -d . -v -in $bam --out $outfn
  [ $? -ne 0 ] && exit 1
done

