#!/bin/bash

### usage: ./runme.sh phasedstr.vcf.gz snp.vcf.gz snp.str.vcf.gz chrom

phasedstr=$1
snp=$2
snpstr=$3
chrom=$4

rm ID.SS*.txt
rm imputed.str.SS*.vcf.gz*
rm SS*.diff.txt

# extract parents
bcftools view ${phasedstr} --samples-file ssc_parents_id.txt -O z -o hipstr.chr${chrom}.phased.parents.vcf.gz --threads 4 --force-samples
bcftools index hipstr.chr${chrom}.phased.parents.vcf.gz
bcftools view ${snp} --samples-file ssc_parents_id.txt -O z -o shapeit.chr${chrom}.parents.vcf.gz --threads 4 --force-samples
bcftools index shapeit.chr${chrom}.parents.vcf.gz

phasedstr=hipstr.chr${chrom}.phased.parents.vcf.gz
snp=shapeit.chr${chrom}.parents.vcf.gz

./concat.sh -s ${phasedstr} -p ${snp} -o ${snpstr}

bcftools query -l ${snpstr} > sample.txt
cat sample.txt | xargs -i -P8 ./L1O.sh -s {} -p pedigree.fam -l ${snpstr} -t ${phasedstr} -b beagle.08Jun17.d8b.jar -q ${snp}

ls *diff*txt | sed 's/.diff.txt//' > samplesID.txt

bcftools query -f '%ID\t%CHROM:%POS\n' ${phasedstr} | sort -f > POS_ID.txt

rm l1o.results.csv
python analyse_L1O.py .

mv l1o.results.csv l1o.results.chr${chrom}.csv

cat sample.txt | awk '{print "imputed.str."$1".vcf.gz"}' > files.list
cat files.list | xargs -i -P8 bcftools index {}
bcftools merge -m id -l files.list -O z -o l1o.imputed.chr${chrom}.str.vcf.gz
bcftools index l1o.imputed.chr${chrom}.str.vcf.gz
