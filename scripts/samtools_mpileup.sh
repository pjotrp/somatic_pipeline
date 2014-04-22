#! /bin/bash
# 
# Simple samtools index loop
#
#   ./index.sh < bamlist
#
# Example
#
#   ls -1 --color=never *_rmdup.bam | samtools_mpileup.sh

samtools=$HOME/opt/bin/samtools
onceonly=$HOME/izip/git/opensource/ruby/once-only/bin/once-only

REFSEQ=/data/GENOMES/human_GATK_GRCh37/GRCh37_gatk.fasta
BED="$HOME/kinome_design_SS_V2_110811_nochr_annot.bed"
phred=10  # 1:1000
bed=$BED
refgenome=$REFSEQ


echo "==== Index fasta with samtools ..."
echo "$samtools faidx $refgenome"|$onceonly --pfff -d . -v
[ $? -ne 0 ] && exit 1

while read bam ; do
  echo ==== $samtools mpileup $bam...
  outfn=$(basename $bam .bam).mpileup
  # check -E option
  echo "$samtools mpileup -B -q $phred -f $refgenome -l $bed $bam > $bam.mpileup"|$onceonly --pfff -v -d . --out $outfn
  [ $? -ne 0 ] && exit 1
  exit 1
done
