#! /bin/bash
# 
# Simple samtools index loop
#
#   ./index.sh < bamlist
#
# Example
#
#   pprins@wgs01:~/data/run5/bam_reduced/mpileup$ ls -1 --color=never ../*rmdup.bam |~/opt/somatic_pipeline/scripts/samtools_mpileup.sh

# ---- Fetch command line and environment
if [ $1 == "--config" ]; then
  config=$2
  shift ; shift
  . $config
fi

samtools=samtools-1.1
onceonly=$HOME/local/bin/once-only

phred=10  # 1:10
bed="$HOME/MBC/kinome_design_SS_V2_110811_nochr_annot_sorted.bed"
orig_refgenome="/hpc/cog_bioinf/GENOMES/human_GATK_GRCh37/GRCh37_gatk.fasta"
refgenome="GRCh37_gatk.fasta"

set

echo "==== Index fasta with samtools ..."

cp -vau $orig_refgenome $refgenome
echo "$samtools faidx $refgenome"|$onceonly --pfff -d . -v
[ $? -ne 0 ] && exit 1

while read bam ; do
  echo ==== $samtools mpileup $bam...
  outfn=$(basename $bam .bam).mpileup
  # Optimized for FFPE after rmdup, settings comparable with bcbio-next
  #
  # Using -E instead of BAQ (no -B)
  echo "$samtools mpileup -E -m 3 -q $phred -l $bed -f $refgenome $bam > $outfn"|$onceonly --pfff --pbs '-q veryshort' --pbs-type SGE -v -d . --out $outfn
  [ $? -ne 0 ] && exit 1
done
