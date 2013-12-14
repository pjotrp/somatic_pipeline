#! /bin/bash
#
# Usage
#
#   varscan2.sh [--config env.sh] normaldescr tumordescr normal.bam tumor.bam
#
# This script is normally run from a controller which creates env.sh. It can also 
# be submittend to PBS with, for example, 'qsub -P SAP42 -cwd'

# ---- Default settings
REFSEQ=/data/GENOMES/human_GATK_GRCh37/GRCh37_gatk.fasta
BED="$HOME/full_kinome_CoDeCZ_chr17.bed"
SAMTOOLS=$HOME/opt/bin/samtools
SAMBAMBA=$HOME/opt/bin/sambamba
ONCEONLY="$HOME/izip/git/opensource/ruby/once-only/bin/once-only"

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

phred=30  # 1:1000
onceonly=$ONCEONLY
refgenome=$REFSEQ
samtools=$SAMTOOLS
sambamba=$SAMBAMBA
bed=$BED

set

mkdir -p varscan2

echo "tumor=$tumor normal=$normal"

for x in $normal $tumor ; do 
  echo "==== Index with Samtools $x ..."
  echo "$sambamba index $x"| $onceonly --pfff -d . -v
  [ $? -ne 0 ] && exit 1
  echo "==== Index fasta with samtools ..."
  echo "$samtools faidx $refgenome"|$onceonly --pfff -d . -v
  [ $? -ne 0 ] && exit 1
  echo "==== Create samtools mpileup of $x"
  # check -E option
  echo "$samtools mpileup -B -q $phred -f $refgenome -l $bed ../$x > $x.mpileup"|$onceonly --pfff -v -d varscan2 --skip $x.mpileup
  [ $? -ne 0 ] && exit 1
done


# --mpileup 1 option (newer undocumented scoring)
echo "java -jar $HOME/opt/lib/VarScan.v2.3.6.jar somatic $normal.mpileup $tumor.mpileup $normal-$tumor.varScan.output --min-coverage-normal 5 --min-coverage-tumor 8 --somatic-p-value 0.001 --strand-filter --min-var-freq 0.20"|$onceonly --pfff -v -d varscan2
[ $? -ne 0 ] && exit 1

echo "java -jar ~/opt/lib/VarScan.v2.3.6.jar processSomatic $normal-$tumor.varScan.output.snp"|$onceonly -v -d varscan2 --skip $normal-$tumor.varScan.output.snp
[ $? -ne 0 ] && exit 1

# The following runs the alternative readcount tools (older scoring)
#
echo "==== Readcount on tumor $tumor..."
echo "~/opt/bin/bam-readcount -b $phred -w 5 -f $refgenome  ../$tumor 17 > $tumor.readcount"|$onceonly --pfff -d varscan2 -v --skip $tumor.readcount
[ $? -ne 0 ] && exit 1
echo "Running fpfilter using $tumor..."
echo "perl $HOME/opt/bin/fpfilter.pl --output-basename $tumor $normal-$tumor.varScan.output.snp $tumor.readcount" | $onceonly --pfff -d varscan2 -v
[ $? -ne 0 ] && exit 1

