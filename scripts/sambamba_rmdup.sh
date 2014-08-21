#! /bin/bash
# 
# Simple rmdup loop
#
#   sambamba_rmdup.sh < bamlist
#
# Example
#
#   ls -1 --color=never *.bam | grep -v _rmdup | /hpc/local/CentOS6/cog_bioinf/CuppenResearch/somatic_pipeline/scripts/sambamba_rmdup.sh

sambamba=sambamba
onceonly=$HOME/izip/git/opensource/ruby/once-only/bin/once-only

while read bam ; do
  echo ==== $sambamba remove duplicates $bam...
  outfn=$(basename $bam .bam)_rmdup.bam
  # echo "$sambamba markdup -r $bam $outfn"| $onceonly --pfff --pbs '-q veryshort' -d . -v -in $bam --out $outfn
  echo "$sambamba markdup -r $bam $outfn"| $onceonly -v --pfff -d . -in $bam --out $outfn
  [ $? -ne 0 ] && exit 1
done

