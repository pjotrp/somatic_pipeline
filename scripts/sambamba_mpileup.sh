#! /bin/bash
# 
# Simple sambamba index loop
#
#   ./index.sh < bamlist
#
# Example
#
#   pprins@wgs01:~/data/run5/bam_reduced/mpileup$ ls -1 --color=never ../*rmdup.bam |~/opt/somatic_pipeline/scripts/sambamba_mpileup.sh

# ---- Fetch command line and environment
if [ $1 == "--config" ]; then
  config=$2
  shift ; shift
  . $config
fi

sambamba=sambamba-1.1
samtools=samtools
onceonly=$HOME/local/bin/once-only

phred=10  # 1:10
bed="$HOME/MBC/kinome_design_SS_V2_110811_nochr_annot_sorted.bed"
orig_refgenome="/hpc/cog_bioinf/GENOMES/human_GATK_GRCh37/GRCh37_gatk.fasta"
refgenome="GRCh37_gatk.fasta"

set

echo "==== Index fasta with samtools ..."

ln -s $orig_refgenome $refgenome
echo "$samtools faidx $refgenome"|$onceonly --pfff -d . -v
[ $? -ne 0 ] && exit 1

while read bam ; do
  echo ==== $sambamba mpileup $bam...
  outfn=$(basename $bam .bam).mpileup
  # Optimized for FFPE after rmdup, settings comparable with bcbio-next
  #
  # Using -E instead of BAQ (no -B)
  echo "$sambamba mpileup $bam --samtools -E -m 3 -q $phred -l $bed -f $refgenome  > $outfn"|$onceonly --pfff -v -d . --out $outfn
  [ $? -ne 0 ] && exit 1
done
