#! /bin/sh

tumor=../merged_MBC023T_F3_20130828_chr17_full_kinome_CoDeCZ_chr17.bam
ref=../merged_MBC023R_F3_20130828_chr17_full_kinome_CoDeCZ_chr17.bam

tumor=../rmdup_merged_MBC051T_F3_20130828_chr17_full_kinome_CoDeCZ_chr17.bam
ref=../rmdup_merged_MBC051R_F3_20130828_chr17_full_kinome_CoDeCZ_chr17.bam

# tumor=../merged_MBC008T_F3_20130528_chr17_full_kinome_CoDeCZ_chr17.bam
# ref=../merged_MBC008R_F3_20130528_chr17_full_kinome_CoDeCZ_chr17.bam
echo ---- Tumor $tumor
echo ---- Reference $ref

name=${tumor#*/}
java6 -Xmx2g -Djava.io.tmpdir=$HOME/tmp -jar ~/opt/mutect/muTect-1.1.4.jar \
  --analysis_type MuTect \
  --reference_sequence /data/GENOMES/human_GATK_GRCh37/GRCh37_gatk.fasta \
  --input_file:normal $ref \
  --input_file:tumor $tumor \
  --out $name.call_stats.txt \
  --coverage_file $name.coverage.wig.txt 

# --dbsnp ~/Broad/dbsnp_132_b37.leftAligned.vcf \
# --intervals 17:78000000-79000000 \
# --cosmic hg19_cosmic_v54_120711.vcf


