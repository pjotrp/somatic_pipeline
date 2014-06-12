#! /bin/bash
#
# Usage
#
#   varscan2.sh [--config env.sh] normalname tumorname normal.mpileup tumor.mpileup
#
# This script is normally run from a controller which creates env.sh. It can also 
# be submittend to PBS with, for example, 'qsub -P SAP42 -cwd'
#
# E.g.
#
#  ./make_paired_tumor_normal_list.rb ~/data/run5/bam_reduced/mpileup/*.mpileup > ~/data/run5/paired_tumor_normal_list.txt
#
#
#
#   ~/opt/somatic_pipeline/scripts/run.rb --pbs --config run.json ~/opt/somatic_pipeline/scripts/varscan2.sh all_mbc.txt
#
# If you want additional VCF output set varscan_vcf=1 in run.json.

# ---- Default settings
CHROMOSOMES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"

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

set

mkdir -p varscan2

# --mpileup 1 option (newer undocumented scoring)

outfn=$normal-$tumor.varScan.output

options="../$normal ../$tumor $outfn --min-coverage-normal 5 --min-coverage-tumor 8 --somatic-p-value 0.001"

# Make sure the inputs are pileups!
echo "java -jar $HOME/opt/lib/VarScan.v2.3.6.jar somatic $options"|$onceonly --pfff --in ../$normal --in ../$tumor --skip-glob "$outfn*" -v -d varscan2
[ $? -ne 0 ] && exit 1

echo "java -jar $HOME/opt/lib/VarScan.v2.3.6.jar processSomatic $outfn.snp"|$onceonly -v -d varscan2 --in $outfn.snp
[ $? -ne 0 ] && exit 1

if [ ! -z $VARSCAN_VCF ]; then
  # Create (optional) VCF output
  echo "java -jar $HOME/opt/lib/VarScan.v2.3.6.jar somatic $options --output-vcf "|$onceonly --pfff --in ../$normal --in ../$tumor --skip-glob "$outfn*" -v -d varscan2
  [ $? -ne 0 ] && exit 1
fi

exit 0

# The following runs the alternative readcount tools (older scoring)
# this is no longer used.
#
echo "==== Readcount on tumor $tumor..."
for chr in $CHROMOSOMES ; do
  echo "!!!! chromosome $chr"
  # By chromosome to avoid readcount segfault!
  echo "~/opt/bin/bam-readcount -b $phred -w 5 -f $refgenome  ../$tumor $chr > $tumor.$chr.readcount"|$onceonly --pfff -d varscan2 -v --skip $tumor.$chr.readcount
  [ $? -ne 0 ] && exit 1
  echo "Running fpfilter using $tumor..."
  echo "perl $HOME/opt/bin/fpfilter.pl --output-basename $tumor.$chr $normal-$tumor.varScan.output.snp $tumor.$chr.readcount" | $onceonly --pfff -d varscan2 -v
  [ $? -ne 0 ] && exit 1
done
