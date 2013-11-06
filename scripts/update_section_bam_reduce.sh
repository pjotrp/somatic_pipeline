#! /bin/bash
#
#  Reduce BAM files by chr number
#
#  Usage:
#
#    update_section_bam_reduce.sh 17 path/*.bam
#
#  The original BAM files can be a full path - the output will be written to the local dir

chr=$1
shift

for x in $* ; do 
   echo $x
   full=${x##*/}
   name=${full%%.*}
   echo ==== Reducing $x 
   echo $name
   # this can also be done with sambamba, provided there is an index
   echo "~/opt/bin/samtools view -bh $x $chr > ${name}_chr$chr.bam" |~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -v -d . 
done

