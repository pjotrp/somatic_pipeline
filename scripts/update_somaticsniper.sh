#! /bin/sh

listfn=$1
if [ -z $listfn ] ; then 
  echo "Provide a tab delimited file with bamfile names" 
  exit 1
fi

BASEDIR="$( cd "$( dirname "$0" )" && pwd )"
echo $BASEDIR

onceonly="$HOME/izip/git/opensource/ruby/once-only/bin/once-only"
refgenome=/data/GENOMES/human_GATK_GRCh37/GRCh37_gatk.fasta
bed="$HOME/full_kinome_CoDeCZ_chr17.bed"

mkdir -p somaticsniper

IFS=$'\t'
while read normal tumor ; do
  echo "tumor=$tumor normal=$normal"
  ref=$normal
  name="${tumor%.*}"

  echo "bam-somaticsniper -f $refgenome ../$tumor ../$ref $name.snp"| ~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -d somaticsniper -v

  echo "==== Index with Samtools $tumor ..."
  echo "samtools index $tumor"| ~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -d . -v

  echo "==== Readcount on tumor $tumor..."
  echo "~/opt/bin/bam-readcount -b 15 -w 5 -f /data/GENOMES/human_GATK_GRCh37/GRCh37_gatk.fasta  ../$tumor 17 > $tumor.readcount"| ~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -d somaticsniper -v --skip $tumor.readcount
  echo "Running fpfilter using ref $ref..."
  echo "perl $HOME/opt/somatic-sniper/src/scripts/fpfilter.pl --output-basename $tumor $name.snp --readcount-file $tumor.readcount"|~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -d somaticsniper -v

