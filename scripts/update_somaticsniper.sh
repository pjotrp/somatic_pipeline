#! /bin/sh

mkdir -p somaticsniper

tumor=merged_MBC023T_F3_20130828_chr17_full_kinome_CoDeCZ_chr17.bam
ref=merged_MBC023R_F3_20130828_chr17_full_kinome_CoDeCZ_chr17.bam
# tumor=merged_MBC008T_F3_20130528_chr17_full_kinome_CoDeCZ_chr17.bam
# ref=merged_MBC008R_F3_20130528_chr17_full_kinome_CoDeCZ_chr17.bam
echo ---- Tumor $tumor
echo ---- Reference $ref
  name="${tumor%.*}"

  echo "bam-somaticsniper -f /data/GENOMES/human_GATK_GRCh37/GRCh37_gatk.fasta ../$tumor ../$ref $name.snp"| ~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -d somaticsniper -v

  echo "==== Index with Samtools $tumor ..."
  echo "samtools index $tumor"| ~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -d . -v

  echo "==== Readcount on tumor $tumor..."
  echo "~/opt/bin/bam-readcount -b 15 -w 5 -f /data/GENOMES/human_GATK_GRCh37/GRCh37_gatk.fasta  ../$tumor 17 > $tumor.readcount"| ~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -d somaticsniper -v
  echo "Running fpfilter using ref $ref..."
  echo "perl $HOME/opt/somatic-sniper/src/scripts/fpfilter.pl --output-basename $tumor $name.snp --readcount-file $tumor.readcount"|~/izip/git/opensource/ruby/once-only/bin/once-only --pfff -d somaticsniper -v $*

