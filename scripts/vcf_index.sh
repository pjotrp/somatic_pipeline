#! /bin/bash
# 
# Simple vcf sort using tabix
#
# e.g.,
#
#   ls *Somatic*vcf --color=never -1 |vcf_index.sh
#

while read infile ; do
  echo ** VCF $infile
  echo *** Sort...
  sorted=$infile.sorted.vcf
  grep '^\#' $infile > $sorted
  egrep -v '^X|^Y|^#' $infile | sort -n -k1 -k2 >> $sorted
  egrep '^X|^Y' $infile | sort -k1,1d -k2,2n >> $sorted
  echo *** bgzip...
  bgzip -f $sorted
  [ $? -ne 0 ] && exit 1
  echo *** tabix index...
  tabix -p vcf -f $sorted.gz
  [ $? -ne 0 ] && exit 1
done
