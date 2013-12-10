#! /bin/bash
#
# Usage
#
#   update_varscan2.sh listfile
#
# where listfile contains a list of tab delimited ref-tumor samples
#

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
    name="${x%.*}_full_kinome_CoDeCZ_chr17"
    echo "==== Create samtools mpileup of $x (name $name)"
    # check -E option
    echo "~/opt/bin/samtools mpileup -B -q 10 -f $refgenome -l $bed ../$x > $x.mpileup"|~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -v -d varscan2
  done

  # --mpileup 1 option
  echo "java -jar ~/opt/lib/VarScan.v2.3.6.jar somatic $normal.mpileup $tumor.mpileup $normal-$tumor.varScan.output --min-coverage-normal 5 --min-coverage-tumor 8 --somatic-p-value 0.001 --strand-filter --min-var-freq 0.20"|~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -v -d varscan2

  echo "java -jar ~/opt/lib/VarScan.v2.3.6.jar processSomatic $normal-$tumor.varScan.output.snp"|~/izip/git/opensource/ruby/once-only/bin/once-only -v -d varscan2

  # The following runs the alternative readcount tools (older?)
  #
  # echo "==== Index with Samtools $tumor2.bam ..."
  # echo "samtools index $tumor2.bam"| ~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -d . -v
  # echo "==== Readcount on tumor $tumor2..."
  # echo "~/opt/bin/bam-readcount -w 5 -f /data/GENOMES/human_GATK_GRCh37/GRCh37_gatk.fasta  ../$tumor2.bam 17 > $tumor2.readcount"| ~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -d varscan2 -v
  # echo "Running fpfilter using ref $ref2..."
  # echo "perl $HOME/opt/bin/fpfilter.pl --output-basename $tumor2 $ref2-varscan.snp $tumor2.readcount"|~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -d varscan2 -v 

done < $listfn



