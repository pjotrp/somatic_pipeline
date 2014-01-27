#! /bin/bash
#
# Usage
#
#   varscan2.sh [--config env.sh] normaldescr tumordescr normal.bam tumor.bam
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
  BED="$HOME/full_kinome_CoDeCZ.bed"
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
normalname=$1
tumorname=$2
normal=$3
tumor=$4

phred=10  # 1:1000
onceonly=$ONCEONLY
refgenome=$REFSEQ
samtools=$SAMTOOLS
sambamba=$SAMBAMBA
bed=$BED

set

mkdir -p varscan2

echo "normal=$normal tumor=$tumor"

for x in $normal $tumor ; do 
  echo "==== Remove duplicates of $x"
  name="${x%.*}"
  x2=${name}_rmdup.bam
  echo "$sambamba markdup -r $x $x2"| $onceonly --pfff -d . -v --skip $x2
  [ $? -ne 0 ] && exit 1
done
normal=${normal%.*}_rmdup.bam
tumor=${tumor%.*}_rmdup.bam
echo "normal=$normal tumor=$tumor"

for x in $normal $tumor ; do 
  echo "==== Index with Samtools $x ..."
  echo "$sambamba index $x"| $onceonly --pfff -d . -v
  [ $? -ne 0 ] && exit 1
  echo "==== Index fasta with samtools ..."
  echo "$samtools faidx $refgenome"|$onceonly --pfff -d . -v
  [ $? -ne 0 ] && exit 1
  echo "==== Create samtools mpileup of $x"
  # check -E option
  # echo "$samtools mpileup -B -q $phred -f $refgenome -l $bed ../$x > $x.mpileup"|$onceonly --pfff -v -d varscan2 --skip $x.mpileup
  echo "$samtools mpileup -B -q $phred -f $refgenome -l $bed ../$x > $x.mpileup"|$onceonly --pfff -v -d varscan2 --skip $x.mpileup
  [ $? -ne 0 ] && exit 1
done

# --mpileup 1 option (newer undocumented scoring)
# options=$normal.mpileup $tumor.mpileup $normal-$tumor.varScan.output --min-coverage-normal 5 --min-coverage-tumor 8 --somatic-p-value 0.001 --strand-filter --min-var-freq 0.20
options="$normal.mpileup $tumor.mpileup $normal-$tumor.varScan.output --min-coverage-normal 5 --min-coverage-tumor 6 --somatic-p-value 0.01 --min-var-freq 0.15"
echo "java -jar $HOME/opt/lib/VarScan.v2.3.6.jar somatic $options"|$onceonly --pfff -v -d varscan2
[ $? -ne 0 ] && exit 1

#   --min-tumor-freq - Minimum variant allele frequency in tumor [0.10]
#   --max-normal-freq - Maximum variant allele frequency in normal [0.05]
#   --p-value - P-value for high-confidence calling [0.07]
echo "java -jar ~/opt/lib/VarScan.v2.3.6.jar processSomatic $normal-$tumor.varScan.output.snp"|$onceonly -v -d varscan2 --skip $normal-$tumor.varScan.output.snp
[ $? -ne 0 ] && exit 1

exit 0

# The following runs the alternative readcount tools (older scoring)
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
