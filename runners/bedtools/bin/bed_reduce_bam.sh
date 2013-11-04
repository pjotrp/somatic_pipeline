# Reduce BAM file by applying bed
#
# bed_reduce_bam.sh design.bed bamfile 

design=$1
bam=$2
base1=`basename $design`

~/opt/bedtools/bin/intersectBed -abam $bam -b $design > ${bam%.*}_${base1%.*}.bam
