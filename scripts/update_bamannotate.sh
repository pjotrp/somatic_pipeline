#! /bin/sh

listfn=$1
if [ -z $listfn ] ; then 
  echo "Provide a tab delimited file with bamfile names" 
  exit 1
fi

bamannotate=/data/common_scripts/stef/getBamAnnotate.pl
inheritance=/data/common_scripts/stef/getInheritance.pl

dir=bamannotate
mkdir -p $dir

BASEDIR="$( cd "$( dirname "$0" )" && pwd )"
echo $BASEDIR

onceonly="$HOME/izip/git/opensource/ruby/once-only/bin/once-only"
refgenome=/data/GENOMES/human_GATK_GRCh37/GRCh37_gatk.fasta
bed="$HOME/full_kinome_CoDeCZ_chr17.bed"

IFS=$'\t'
while read ref tumor ; do
  echo "tumor=$tumor ref=$ref"
  echo ---- Tumor $tumor
  echo ---- Reference $ref
  tumorname="${tumor%.*}"
  refname="${ref%.*}"
  tumorout=${tumorname##*/}
  refout=${refname##*/}

  tumorbam=`find /data/mapping/cancer -name $tumor`
  tumorbase="${tumorbam%.*}"
  refbam=`find /data/mapping/cancer -name $ref`
  refbase="${refbam%.*}"
  echo "*** $tumorbam ***"
  echo "*** $refbam ***"
  echo "*** $tumorout ***"
  echo "*** $refout ***"
  echo "$bamannotate -out $tumorout -snv $refbase*.snvs.snv1 -snv $tumorbase*.snvs.snv1 -ref /data/GENOMES/human_GATK_GRCh37/GRCh37_gatk.fasta -bam $refbam -bam $tumorbam"|  ~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -d $dir -v 
  echo "$inheritance -d -in $tumorout.bamAnn -out $tumorout -uaf_ids 3 -aff_ids 4 -min_cov 8"|~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -d $dir -v 
done < $listfn

