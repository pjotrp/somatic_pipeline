#! /bin/bash
# 
# Simple annotator loop for .vcf and .snp files
#
#   ./annotator.sh < list
#

for x in *.snp ; do echo $x; ~/izip/git/opensource/ruby/once-only/bin/once-only /data/common_scripts/SAP42-testing/annotator.pl -in $x ; done

annotate=/data/common_scripts/SAP42-testing/annotator.pl
onceonly=$HOME/izip/git/opensource/ruby/once-only/bin/once-only

while read infile ; do
  echo ==== $annotate $infile...
  echo "$annotate -in $infile"|$onceonly --pfff -d . -v
  [ $? -ne 0 ] && exit 1
done
