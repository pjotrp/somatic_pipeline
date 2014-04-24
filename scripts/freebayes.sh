#! /bin/bash
#
# Usage
#
#   freebayes.sh [--config env.sh] normalname tumorname normal.bam tumor.bam
#
# This script is normally run from a controller which creates env.sh. It can also 
# be submittend to PBS with, for example, 'qsub -P SAP42 -cwd'
#
# E.g.
#
#  ~/opt/somatic_pipeline/scripts/make_paired_tumor_normal_list.rb /home/cog/pprins/data/run5/bam_reduced/*p.bam > ../../../paired_tumor_normal_bamlist.txt
#  ~/opt/somatic_pipeline/scripts/run.rb --config ../../run.json ~/opt/somatic_pipeline/scripts/freebayes.sh ../../paired_tumor_normal_bamlist.txt
#
#
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

for bam in $normal $tumor ; do
  echo ==== $sambamba index $bam...
  outfn1=$(basename $bam .bam).bam.bai
  echo "$sambamba index $bam"| $onceonly --pfff -d . -v -in $bam --out $outfn1
  [ $? -ne 0 ] && exit 1
done


outfn=$normal-$tumor.freebayes.output
options="-f $refgenome -C 3 -t $bed --pooled-discrete --genotype-qualities --min-coverage 5"
echo "$freebayes $options ../$normal ../$tumor > $outfn.vcf "|$onceonly --pfff --in ../$normal --in ../$tumor --out $outfn.vcf -v -d freebayes
[ $? -ne 0 ] && exit 1

samples=`~/izip/git/opensource/ruby/bioruby-vcf/bin/bio-vcf -q --eval-once 'header.samples.join(" ")' < freebayes/$outfn.vcf`

echo "##### "$samples
echo "$HOME/opt/vcflib/bin/vcfsamplediff VT $samples $outfn.vcf > $outfn.diff.vcf"|$onceonly --pfff --in $outfn.vcf --out $outfn.diff.vcf -v -d freebayes
[ $? -ne 0 ] && exit 1

grep -i VT=germline freebayes/$outfn.diff.vcf > freebayes/$outfn.Germline.vcf
grep -i VT=Somatic freebayes/$outfn.diff.vcf > freebayes/$outfn.Somatic.vcf
grep -i VT=LOH freebayes/$outfn.diff.vcf > freebayes/$outfn.LOH.vcf

exit 0

