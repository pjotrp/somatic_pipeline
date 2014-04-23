#! /bin/bash
#
# Usage
#
#   varscan2.sh [--config env.sh] normalname tumorname normal.mpileup tumor.mpileup
#
# This script is normally run from a controller which creates env.sh. It can also 
# be submittend to PBS with, for example, 'qsub -P SAP42 -cwd'
#
# E.g.
#
#   ~/opt/somatic_pipeline/scripts/run.rb --pbs --config run.json ~/opt/somatic_pipeline/scripts/varscan2.sh all_mbc.txt
#
# Uncomment for testing:
# CHR=17

# ---- Default settings
REFSEQ=/data/GENOMES/human_GATK_GRCh37/GRCh37_gatk.fasta
if [ ! -z $CHR ]; then
  BED="$HOME/full_kinome_CoDeCZ_chr$CHR.bed"
else
  # BED="$HOME/full_kinome_CoDeCZ.bed"
  BED="$HOME/kinome_design_SS_V2_110811_nochr_annot.bed"
fi
SAMTOOLS=$HOME/opt/bin/samtools
SAMBAMBA=$HOME/opt/bin/sambamba
ONCEONLY="$HOME/izip/git/opensource/ruby/once-only/bin/once-only"
CHROMOSOMES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"

# ---- PBS settings
#$ -S /bin/bash
#$ -o stdout
PATH=$SGE_O_PATH:$PATH

# ---- Fetch command line and environment
if [ $1 == "--config" ]; then
  config=$2
  shift ; shift
  . $config
fi
normalname=$1 # unused
tumorname=$2  # unused
normal=$3
tumor=$4

phred=30  # 1:1000
onceonly=$ONCEONLY
refgenome=$REFSEQ
samtools=$SAMTOOLS
sambamba=$SAMBAMBA
bed=$BED

set

mkdir -p varscan2

# --mpileup 1 option (newer undocumented scoring)

outfn=$normal-$tumor.varScan.output

options="../$normal ../$tumor $outfn --min-coverage-normal 5 --min-coverage-tumor 8 --somatic-p-value 0.001"

echo "java -jar $HOME/opt/lib/VarScan.v2.3.6.jar somatic $options"|$onceonly --pfff --in ../$normal --in ../$tumor --skip-glob "$outfn*" -v -d varscan2
[ $? -ne 0 ] && exit 1

echo "java -jar $HOME/opt/lib/VarScan.v2.3.6.jar processSomatic $outfn.snp"|$onceonly -v -d varscan2 --in $outfn.snp
[ $? -ne 0 ] && exit 1

exit 0

# The following runs the alternative readcount tools (older scoring)
# this is no longer used.
#
echo "==== Readcount on tumor $tumor..."
for chr in $CHROMOSOMES ; do
  echo "!!!! chromosome $chr"
  # By chromosome to avoid readcount segfault!
  echo "~/opt/bin/bam-readcount -b $phred -w 5 -f $refgenome  ../$tumor $chr > $tumor.$chr.readcount"|$onceonly --pfff -d varscan2 -v --skip $tumor.$chr.readcount
  [ $? -ne 0 ] && exit 1
  echo "Running fpfilter using $tumor..."
  echo "perl $HOME/opt/bin/fpfilter.pl --output-basename $tumor.$chr $normal-$tumor.varScan.output.snp $tumor.$chr.readcount" | $onceonly --pfff -d varscan2 -v
  [ $? -ne 0 ] && exit 1
done
