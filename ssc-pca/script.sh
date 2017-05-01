#!/bin/bash

### input - 1kg ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
### input - SSC - ssc.chr21.vcf.gz
### change population labels in plot_pca.py

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
bcftools view ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --samples-file ref.samples.txt --output-type z --output-file 1000g.ref.vcf.gz

bcftools index 1000g.ref.vcf.gz
bcftools merge -m all 1000g.ref.vcf.gz ssc.chr21.vcf.gz --output-type z --output chr21.combined.vcf.gz
plink --vcf chr21.combined.vcf.gz --geno 0.3 --make-bed --out chr21.combined

plink --bfile chr21.combined --pca 10 --out pca_10
python plot_pca.py
