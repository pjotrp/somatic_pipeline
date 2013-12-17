#! /bin/bash
#
# Usage
#
#   bamannotate.sh [--config env.sh] normaldescr tumordescr normal.bam tumor.bam
#
# This script is normally run from a controller which creates env.sh. It can also 
# be submittend to PBS with, for example, 'qsub -P SAP42 -cwd'

# Uncomment for testing:
# CHR=17

# ---- Default settings
normalSEQ=/data/GENOMES/human_GATK_GRCh37/GRCh37_gatk.fasta
if [ ! -z $CHR ]; then
  BED="$HOME/full_kinome_CoDeCZ_chr$CHR.bed"
else
  BED="$HOME/full_kinome_CoDeCZ.bed"
fi
SAMTOOLS=$HOME/opt/bin/samtools
SAMBAMBA=$HOME/opt/bin/sambamba
BAMANNOTATE=/data/common_scripts/stef/getBamAnnotate.pl
INHERITANCE=/data/common_scripts/stef/getInheritance.pl
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
normalfull=$5
tumorfull=$6

phred=30  # 1:1000
onceonly=$ONCEONLY
normalgenome=$normalSEQ
samtools=$SAMTOOLS
sambamba=$SAMBAMBA
bed=$BED
bamannotate=$BAMANNOTATE
inheritance=$INHERITANCE
cachedir=`pwd`  # otherwise set to `pwd`

set

dir=./bamannotate
mkdir -p ./bamannotate

echo "normal=$normal tumor=$tumor"
normalbase=`dirname $normalfull`
tumorbase=`dirname $tumorfull`
echo "*** normal base $normalbase ***"
echo "*** tumor base $tumorbase ***"
tumorout=${tumorname##*/}
normalout=${normalname##*/}
echo "*** normal out $normalout ***"
echo "*** tumor out $tumorout ***"
  
echo "$bamannotate -out $tumorout -snv $normalbase/*.snvs.snv1 -snv $tumorbase/*.snvs.snv1 -normal /data/GENOMES/human_GATK_GRCh37/GRCh37_gatk.fasta -bam $cachedir/$normal -bam $cachedir/$tumor"|  ~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -d $dir -v 
echo "$inheritance -d -in $tumorout.bamAnn -out $tumorout -uaf_ids 3 -aff_ids 4 -min_cov 8"|~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -d $dir -v 


