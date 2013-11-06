#! /bin/bash
#
#  Reindex BAM files
#
#  Usage:
#
#    update_index_bam.sh *.bam
#

for x in $* ; do 
   echo ==== Create index for $x 
   echo "sambamba index $x" |~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -v -d . 
done

