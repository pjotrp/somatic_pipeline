#! /bin/bash
#
# Usage
#
#   update_varscan2.sh listfile
#
# where listfile contains a list of tab delimited ref-tumor samples
#

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

mkdir -p varscan2

IFS=$'\t'
while read normal tumor ; do
  echo "tumor=$tumor normal=$normal"
  ref=$normal
  for x in $tumor $ref ; do 
    echo "==== Index with Samtools $x ..."
    echo "samtools index $x"| ~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -d . -v
    name="${x%.*}_full_kinome_CoDeCZ_chr17"
    echo "==== Create samtools mpileup of $x (name $name)"
    # check -E option
    echo "~/opt/bin/samtools mpileup -B -q $phred -f $refgenome -l $bed ../$x > $x.mpileup"|~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -v -d varscan2 --skip $x.mpileup
  done

  # --mpileup 1 option
  echo "java -jar ~/opt/lib/VarScan.v2.3.6.jar somatic $normal.mpileup $tumor.mpileup $normal-$tumor.varScan.output --min-coverage-normal 5 --min-coverage-tumor 8 --somatic-p-value 0.001 --strand-filter --min-var-freq 0.20"|~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -v -d varscan2

  echo "java -jar ~/opt/lib/VarScan.v2.3.6.jar processSomatic $normal-$tumor.varScan.output.snp"|~/izip/git/opensource/ruby/once-only/bin/once-only -v -d varscan2 --skip $normal-$tumor.varScan.output.snp

  # The following runs the alternative readcount tools (older?)
  #
  echo "==== Readcount on tumor $tumor (chr17)..."
  echo "~/opt/bin/bam-readcount -b $phred -w 5 -f $refgenome  ../$tumor 17 > $tumor.readcount"| ~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -d varscan2 -v --skip $tumor.readcount
  echo "Running fpfilter using ref $ref..."
  echo "perl $HOME/opt/bin/fpfilter.pl --output-basename $tumor $normal-$tumor.varScan.output.snp $tumor.readcount"|~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -d varscan2 -v

done < $listfn



