#! /bin/bash
#
# Usage
#
#   freebayes.sh [--config env.sh] normalname tumorname normal.mpileup tumor.mpileup
#
# This script is normally run from a controller which creates env.sh. It can also 
# be submittend to PBS with, for example, 'qsub -P SAP42 -cwd'
#
# E.g.
#
#   ~/opt/somatic_pipeline/scripts/run.rb --pbs --config run.json ~/opt/somatic_pipeline/scripts/freebayes.sh all_mbc.txt
#

# ---- Default settings
REFSEQ=/data/GENOMES/human_GATK_GRCh37/GRCh37_gatk.fasta
BED="$HOME/kinome_design_SS_V2_110811_nochr_annot_sorted.bed"

SAMTOOLS=$HOME/opt/bin/samtools
SAMBAMBA=$HOME/opt/bin/sambamba
ONCEONLY="$HOME/izip/git/opensource/ruby/once-only/bin/once-only"
FREEBAYES=$HOME/opt/freebayes/bin/freebayes

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
normalname=$1 # unused
tumorname=$2  # unused
normal=$3
tumor=$4

phred=30  # 1:1000
onceonly=$ONCEONLY
refgenome=$REFSEQ
samtools=$SAMTOOLS
sambamba=$SAMBAMBA
freebayes=$FREEBAYES
bed=$BED

set

mkdir -p freebayes

outfn=$normal-$tumor.freebayes.output

options="-f $refgenome -C 3 -t $bed --pooled-discrete --genotype-qualities --min-coverage 5"
echo "$freebayes $options ../$tumor ../$normal $outfn.vcf > "|$onceonly --pfff --in ../$normal --in ../$tumor --skip-glob $outfn.vcf -v -d freebayes
[ $? -ne 0 ] && exit 1

