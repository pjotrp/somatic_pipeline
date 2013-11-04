#! /bin/bash

BASEDIR="$( cd "$( dirname "$0" )" && pwd )"
echo $BASEDIR

onceonly="~/izip/git/opensource/ruby/once-only/bin/once-only"
mkdir -p varscan2

IFS=$'\t'
while read ref tumor ; do
  echo "tumor=$tumor ref=$ref"
  for x in $tumor $ref ; do 
    name="${x%.*}_full_kinome_CoDeCZ_chr17"
    echo Samtools mpileup of $name...
    echo "~/opt/bin/samtools mpileup -B -f /data/GENOMES/human_GATK_GRCh37/GRCh37_gatk.fasta ../$name.bam > $name.mpileup"|~/izip/git/opensource/ruby/once-only/bin/once-only -v -d varscan2
  done
  ref="${ref%.*}_full_kinome_CoDeCZ_chr17"
  tumor="${tumor%.*}_full_kinome_CoDeCZ_chr17"
  echo "Running varscan2..."
  echo "java -jar ~/opt/lib/VarScan.v2.3.6.jar somatic $ref.mpileup $tumor.mpileup $ref-varscan --min-coverage 20 --somatic-p-value 0.001" | ~/izip/git/opensource/ruby/once-only/bin/once-only -d varscan2 -v 

  @@ Check tumor & ref:
  echo "Readcount..."
  echo "~/opt/bin/bam-readcount  -f /data/GENOMES/human_GATK_GRCh37/GRCh37_gatk.fasta  ../$tumor.bam 17 > $tumor.readcount"| ~/izip/git/opensource/ruby/once-only/bin/once-only -d varscan 2 -v
  echo "Running fpfilter..."
  echo "perl $HOME/opt/bin/fpfilter.pl $ref-varscan.snp $tumor.readcount"|~/izip/git/opensource/ruby/once-only/bin/once-only -d varscan2 -v 

  exit 1
done < somatic_bams.txt



