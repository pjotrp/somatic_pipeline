#! /bin/bash
#
# Usage
#
#   freebayes.sh [--config env.sh] normalname tumorname normal.bam tumor.bam
#
# This script is normally run from a controller which creates env.sh. It can also 
# be submitted to PBS with, for example, 'qsub -P SAP42 -cwd'
#
# E.g.
#
# ~/opt/somatic_pipeline/scripts/make_paired_tumor_normal_list.rb `pwd`/bam_reduced/*bam > paired_tumor_normal_bamlist.txt
#  ~/opt/somatic_pipeline/scripts/run.rb --config run.json ~/opt/somatic_pipeline/scripts/freebayes.sh paired_tumor_normal_bamlist.txt
#

# ---- PBS settings
#$ -S /bin/bash
#$ -o stdout
# PATH=$SGE_O_PATH:$PATH

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

set

mkdir -p freebayes

for bam in $normal $tumor ; do
  echo ==== $sambamba index $bam...
  outfn1=$(basename $bam .bam).bam.bai
  echo "$sambamba index $bam"| $onceonly --pfff -d . -v -in $bam --out $outfn1
  [ $? -ne 0 ] && exit 1
done

outfn=$normal-$tumor.freebayes.output
options="-f $refgenome -C 3 -t $bed --pooled-discrete --genotype-qualities --min-coverage 5 --no-indels --no-mnps --no-complex"
echo "$freebayes $options ../$normal ../$tumor > $outfn.vcf "|$onceonly --pfff --in ../$normal --in ../$tumor --out $outfn.vcf -v -d freebayes

samples=`$HOME/izip/git/opensource/ruby/bioruby-vcf/bin/bio-vcf -q --skip-header --eval-once 'header.samples.join(" ")' < freebayes/$outfn.vcf`
# [ $? -ne 0 ] && exit 1  <-- note freebayes always generates an exit error

echo "##### "$samples
echo "$HOME/opt/vcflib/bin/vcfsamplediff VT $samples $outfn.vcf > $outfn.diff.vcf"|$onceonly --pfff --in $outfn.vcf --out $outfn.diff.vcf -v -d freebayes
[ $? -ne 0 ] && exit 1

head -57 freebayes/$outfn.diff.vcf > freebayes/$outfn.Germline.vcf
head -57 freebayes/$outfn.diff.vcf > freebayes/$outfn.Somatic.vcf
grep -i VT=germline freebayes/$outfn.diff.vcf >> freebayes/$outfn.Germline.vcf
grep -i VT=Somatic freebayes/$outfn.diff.vcf >> freebayes/$outfn.Somatic.vcf

exit 0

