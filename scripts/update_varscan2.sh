#! /bin/bash

BASEDIR="$( cd "$( dirname "$0" )" && pwd )"
echo $BASEDIR

mkdir -p varscan2

IFS=$'\t'
while read ref tumor ; do
  echo "tumor=$tumor ref=$ref"
  for x in $tumor $ref ; do 
    # echo "$BASEDIR/../runners/bedtools/bin/bed_reduce_bam.sh $HOME/full_kinome_CoDeCZ_chr17.bed $x"|~/izip/git/opensource/ruby/once-only/bin/once-only -v --ignore-lock -d . 
  # --pbs "-P SAP42" 
    name="${x%.*}_full_kinome_CoDeCZ_chr17"
    echo Varscan on $name...
    echo "~/opt/bin/samtools mpileup -B -f /data/GENOMES/human_GATK_GRCh37/GRCh37_gatk.fasta ../$name.bam > $name.mpileup"|~/izip/git/opensource/ruby/once-only/bin/once-only -v -d varscan2
  done
done < somatic_bams.txt



