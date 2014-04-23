#! /bin/bash
# 
# Simple samtools index loop
#
#   ./index.sh < bamlist
#
# Example
#
#   pprins@wgs01:~/data/run5/bam_reduced/mpileup$ ls -1 --color=never ../*rmdup.bam |~/opt/somatic_pipeline/scripts/samtools_mpileup.sh

samtools=$HOME/opt/bin/samtools
onceonly=$HOME/izip/git/opensource/ruby/once-only/bin/once-only

REFSEQ=/data/GENOMES/human_GATK_GRCh37/GRCh37_gatk.fasta
BED="$HOME/kinome_design_SS_V2_110811_nochr_annot_sorted.bed"
phred=10  # 1:10
bed=$BED
refgenome=$REFSEQ

echo "==== Index fasta with samtools ..."
echo "$samtools faidx $refgenome"|$onceonly --pfff -d . -v
[ $? -ne 0 ] && exit 1

while read bam ; do
  echo ==== $samtools mpileup $bam...
  outfn=$(basename $bam .bam).mpileup
  # Optimized for FFPE after rmdup, settings comparable with bcbio-next
  #
  # Using -E instead of BAQ (no -B)
  echo "$samtools mpileup -E -m 3 -q $phred -l $BED -f $refgenome $bam > $outfn"|$onceonly --pfff -v -d . --out $outfn
  [ $? -ne 0 ] && exit 1
done
