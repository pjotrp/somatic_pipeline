#! /bin/sh

bamannotate=/data/common_scripts/stef/getBamAnnotate.pl
inheritance=/data/common_scripts/stef/getInheritance.pl

dir=bamannotate
mkdir -p $dir

IFS=$'\t'
while read ref tumor ; do
  echo "tumor=$tumor ref=$ref"
  echo ---- Tumor $tumor
  echo ---- Reference $ref
  tumorname="${tumor%.*}"
  refname="${ref%.*}"
  tumorout=${tumorname##*/}
  refout=${refname##*/}

    echo "*** $tumorout ***"
    echo "*** $refout ***"
    echo "$bamannotate -out $tumorout -snv $refname*.snvs.snv1 -snv $tumorname*.snvs.snv1 -ref /data/GENOMES/human_GATK_GRCh37/GRCh37_gatk.fasta -bam $ref -bam $tumor"|  ~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -d $dir -v --ignore-lock
    echo "$inheritance -d -in $tumorout.bamAnn -out $tumorout -uaf_ids 3 -aff_ids 4 -min_cov 4"|~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -d $dir -v --ignore-lock
done < bamfiles_2.txt

