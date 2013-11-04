# Reduce BAM file by applying bed
#
# bed_reduce_bam.sh design.bed bamfile 

design=$1
bam=$2

~/opt/bedtools/bin/intersectBed -abam $bam -b $bed > ${bam%.*}_bed.bam
