#! /bin/sh
#
# Follows the recommendations on http://gmt.genome.wustl.edu/somatic-sniper/1.0.2/documentation.html
#
# Note: the indel check is not yet added

phred=30  # 1:1000

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

  echo "bam-somaticsniper -q $phred -f $refgenome ../$tumor ../$ref $name.snp"| ~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -d somaticsniper -v --skip $name.snp

  echo "==== Index with Samtools $tumor ..."
  echo "samtools index $tumor"| ~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -d . -v

  echo "==== Readcount on tumor $tumor..."
  echo "~/opt/bin/bam-readcount -b $phred -w 5 -f $refgenome  ../$tumor 17 > $tumor.readcount"| ~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -d somaticsniper -v --skip $tumor.readcount
  echo "Running fpfilter using ref $ref..."
  echo "perl $HOME/opt/somatic-sniper/src/scripts/fpfilter.pl --output-basename $tumor --snp-file $name.snp --readcount-file $tumor.readcount"|~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -d somaticsniper -v

  echo "perl $HOME/opt/somatic-sniper/src/scripts/highconfidence.pl --min-mapping-quality $phred --snp-file $tumor.fp_pass"|$onceonly -d somaticsniper -v

done < $listfn

