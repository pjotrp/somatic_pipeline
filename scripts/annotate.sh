#! /bin/bash
# 
# Simple annotator loop for .vcf and .snp files
#
#   ./annotator.sh < list
#

annotate=/data/common_scripts/SAP42-testing/annotator.pl
onceonly=$HOME/izip/git/opensource/ruby/once-only/bin/once-only

while read infile ; do
  echo ==== $annotate $infile...
  echo "$annotate -in $infile"|$onceonly -in $infile --out $infile.snv1 --out $infile.snv2 --pfff -d . -v
  [ $? -ne 0 ] && exit 1
done
