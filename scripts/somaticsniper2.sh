#! /bin/bash
#
# Usage
#
#   somaticsniper.sh [--config env.sh] normaldescr tumordescr normal.bam tumor.bam
#
# This script is normally run from a controller which creates env.sh. It can also 
# be submittend to PBS with, for example, 'qsub -P SAP42 -cwd'

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

phred=30  # 1:1000
onceonly=$ONCEONLY
refgenome=$REFSEQ
samtools=$SAMTOOLS
sambamba=$SAMBAMBA
bed=$BED
somaticsniper=bam-somaticsniper
cachedir=/tmp  # otherwise set to `pwd`

set

mkdir -p somaticsniper

echo "normal=$normal tumor=$tumor"
normalname="${normal%.*}"
tumorname="${tumor%.*}"

df -h

if true ; then 
  for x in $normal $tumor ; do 
    echo "==== Remove duplicates of $x"
    name="${x%.*}"
    x2=$cachedir/${name}_rmdup.bam
    $sambamba markdup -r $x $x2
    # Don't use once-only here with cache
    # echo "$sambamba markdup -r $x $x2"| $onceonly --pfff -d . -v --skip $x2
    [ $? -ne 0 ] && exit 1
  done
  normal=$cachedir/${normal%.*}_rmdup.bam
  tumor=$cachedir/${tumor%.*}_rmdup.bam
  echo "normal=$normal tumor=$tumor"
fi

if true ; then 
  for x in $normal $tumor ; do 
    echo "==== Select design $x"
    name="${x%.*}"
    x2=${name}_bed.bam
    echo "$HOME/opt/bedtools/bin/intersectBed -abam $x -b $bed > $x2"| $onceonly --pfff -d . -v --skip $x2
    # Don't use once-only here with cache
    # $HOME/opt/bedtools/bin/intersectBed -abam $x -b $bed > $x2
    [ $? -ne 0 ] && exit 1
  done
  normal=${normal%.*}_bed.bam
  tumor=${tumor%.*}_bed.bam
  echo "normal=$normal tumor=$tumor"
fi

for x in $normal $tumor ; do 
  echo "==== Index with Samtools $x ..."
  echo "$sambamba index $x"| $onceonly --pfff -d . -v
  [ $? -ne 0 ] && exit 1
  echo "==== Index fasta with samtools ..."
  echo "$samtools faidx $refgenome"|$onceonly --pfff -d . -v
  [ $? -ne 0 ] && exit 1
done

  echo "==== Somatic sniper"
  echo "$somaticsniper -q $phred -Q $phred -J -N 8 -f $refgenome $tumor $normal $normalname-$tumorname.snp"| $onceonly --pfff -d somaticsniper -v --skip $normalname-$tumorname.snp
  [ $? -ne 0 ] && exit 1

# The following runs readcount 
#
echo "==== Readcount on tumor $tumor..."
# CHROMOSOMES="17 18 19 20"
# CHROMOSOMES="17 18"
for chr in $CHROMOSOMES ; do
  echo "!!!! chromosome $chr"
  # By chromosome to avoid readcount segfault!
  echo "~/opt/bin/bam-readcount -b $phred -w 5 -f $refgenome  $tumor $chr > $tumorname.$chr.readcount"|$onceonly --pfff -d somaticsniper -v --skip $tumorname.$chr.readcount --force
  [ $? -ne 0 ] && exit 1
  echo "Running fpfilter using $tumor..."
  echo "perl $HOME/opt/somatic-sniper/src/scripts/fpfilter.pl --output-basename $tumorname.$chr --snp-file $normalname-$tumorname.snp --readcount-file $tumorname.$chr.readcount"|~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -d somaticsniper -v
  [ $? -ne 0 ] && exit 1
  echo "perl $HOME/opt/somatic-sniper/src/scripts/highconfidence.pl --min-mapping-quality $phred --snp-file $tumorname.$chr.fp_pass"|$onceonly -d somaticsniper -v
  # [ $? -ne 0 ] && exit 1  -- throws error on empty fp_pass!
done

