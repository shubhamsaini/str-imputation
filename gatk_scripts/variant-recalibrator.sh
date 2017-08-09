#!/bin/bash

#SBATCH -A ddp268
#SBATCH -p shared
#SBATCH --get-user-env
#SBATCH -t 1-20:00 # Walltime Days-Hours-Minutes
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH -o /home/s1saini/chr15.post.out
#SBATCH -e /home/s1saini/chr15.post.err

c=15

ls -d1v /oasis/scratch/comet/s1saini/temp_project/chr$c/*.gz > /home/s1saini/chr$c.cat.list


java -cp /home/s1saini/gatk.jar org.broadinstitute.gatk.tools.CatVariants \
    -R /oasis/scratch/comet/s1saini/temp_project/ref/human_g1k_v37.fasta \
    -V /home/s1saini/chr$c.cat.list \
    -out /oasis/scratch/comet/s1saini/temp_project/chr$c/chr$c.catvar.vcf.gz \
    -assumeSorted
    
java -jar /home/s1saini/gatk.jar \
    -R /oasis/scratch/comet/s1saini/temp_project/ref/human_g1k_v37.fasta \
    -T VariantAnnotator \
    -V /oasis/scratch/comet/s1saini/temp_project/chr$c/chr$c.catvar.vcf.gz \
    -o /oasis/scratch/comet/s1saini/temp_project/chr$c/chr$c.catvar.dbsnp.vcf.gz \
    --dbsnp /oasis/scratch/comet/s1saini/temp_project/ref/dbsnp_138.b37.vcf


java -jar /home/s1saini/gatk.jar \
   -T VariantRecalibrator \
   -R /oasis/scratch/comet/s1saini/temp_project/ref/human_g1k_v37.fasta \
   -input /oasis/scratch/comet/s1saini/temp_project/chr$c/chr$c.catvar.dbsnp.vcf.gz \
   -recalFile /oasis/scratch/comet/s1saini/temp_project/chr$c/output.recal \
   -tranchesFile /oasis/scratch/comet/s1saini/temp_project/chr$c/output.tranches \
   -nt 4 \
   -resource:hapmap,known=false,training=true,truth=true,prior=15.0 /oasis/scratch/comet/s1saini/temp_project/ref/hapmap_3.3.b37.vcf \
   -resource:omni,known=false,training=true,truth=true,prior=12.0 /oasis/scratch/comet/s1saini/temp_project/ref/1000G_omni2.5.b37.vcf \
   -resource:1000G,known=false,training=true,truth=false,prior=10.0 /oasis/scratch/comet/s1saini/temp_project/ref/1000G_phase1.snps.high_confidence.b37.vcf \
   -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /oasis/scratch/comet/s1saini/temp_project/ref/dbsnp_138.b37.vcf \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -an InbreedingCoeff \
   -mode SNP
   
   
java -jar /home/s1saini/gatk.jar \
   -T ApplyRecalibration \
   -R /oasis/scratch/comet/s1saini/temp_project/ref/human_g1k_v37.fasta \
   -input /oasis/scratch/comet/s1saini/temp_project/chr$c/chr$c.catvar.dbsnp.vcf.gz \
   -tranchesFile /oasis/scratch/comet/s1saini/temp_project/chr$c/output.tranches \
   -recalFile /oasis/scratch/comet/s1saini/temp_project/chr$c/output.recal \
   -o /oasis/scratch/comet/s1saini/temp_project/chr$c/chr$c.output.recalibrated.filtered.vcf.gz \
    --ts_filter_level 99.5 \
   -mode SNP \
   -ef
   

bcftools view /oasis/scratch/comet/s1saini/temp_project/chr$c/chr$c.output.recalibrated.filtered.vcf.gz -f PASS -m2 -M2 -v snps -O u | bcftools annotate -I +'%CHROM:%POS:%REF:%ALT' -O z -o /oasis/scratch/comet/s1saini/temp_project/chr$c/chr$c.recalibrated.pass.dbsnp.annotated.vcf.gz

bcftools index /oasis/scratch/comet/s1saini/temp_project/chr$c/chr$c.recalibrated.pass.dbsnp.annotated.vcf.gz

mkdir /oasis/scratch/comet/s1saini/temp_project/chr$c/shapeit/

bcftools query -f '%POS\n' /oasis/scratch/comet/s1saini/temp_project/chr$c/chr$c.recalibrated.pass.dbsnp.annotated.vcf.gz > /oasis/scratch/comet/s1saini/temp_project/chr$c/shapeit/chr$c.POS.txt

sed -n '0~80000p' /oasis/scratch/comet/s1saini/temp_project/chr$c/shapeit/chr$c.POS.txt > /oasis/scratch/comet/s1saini/temp_project/chr$c/shapeit/chr$c.POS.sed.txt

(echo "1"; cat /oasis/scratch/comet/s1saini/temp_project/chr$c/shapeit/chr$c.POS.sed.txt; tail -1 /oasis/scratch/comet/s1saini/temp_project/chr$c/shapeit/chr$c.POS.txt) > /oasis/scratch/comet/s1saini/temp_project/chr$c/shapeit/chr$c.POS.v2.txt

readarray -t POS < /oasis/scratch/comet/s1saini/temp_project/chr$c/shapeit/chr$c.POS.v2.txt

nlines=`wc -l < /oasis/scratch/comet/s1saini/temp_project/chr$c/shapeit/chr$c.POS.v2.txt`
let "nlines=$nlines-2"

for i in `seq 0 $nlines`;
do

bcftools view /oasis/scratch/comet/s1saini/temp_project/chr$c/chr$c.recalibrated.pass.dbsnp.annotated.vcf.gz -r $c:${POS[$i]}-${POS[$i+1]} --output-type z --output-file /oasis/scratch/comet/s1saini/temp_project/chr$c/shapeit/chr$c.reg$i.vcf.gz



bcftools query -f '%REF\t%ID\n' /oasis/scratch/comet/s1saini/temp_project/chr$c/shapeit/chr$c.reg$i.vcf.gz > /oasis/scratch/comet/s1saini/temp_project/chr$c/shapeit/chr$c.reg$i.a1allele.txt

plink --vcf /oasis/scratch/comet/s1saini/temp_project/chr$c/shapeit/chr$c.reg$i.vcf.gz --keep-allele-order --a1-allele /oasis/scratch/comet/s1saini/temp_project/chr$c/shapeit/chr$c.reg$i.a1allele.txt 1 2 --make-bed --out /oasis/scratch/comet/s1saini/temp_project/chr$c/shapeit/chr$c.reg$i

python /home/s1saini/createFam.py /home/s1saini/famDetail.txt /oasis/scratch/comet/s1saini/temp_project/chr$c/shapeit/chr$c.reg$i.fam /oasis/scratch/comet/s1saini/temp_project/chr$c/shapeit/chr$c.reg$i.fam

plink --bfile /oasis/scratch/comet/s1saini/temp_project/chr$c/shapeit/chr$c.reg$i --keep-allele-order --me 1 1 --set-me-missing --make-bed --out /oasis/scratch/comet/s1saini/temp_project/chr$c/shapeit/chr$c.reg$i\_me

done