# usage: ./gen-snp-str-vcf.sh shapeit.output.vcf.gz final-str-snp.vcf.gz

bcftools query -l $1 > sampleID.txt
bcftools view str-snp-synced.vcf --samples-file sampleID.txt --output-type z --output-file str-snp-synced-reorder.vcf.gz
bcftools view $1 --regions-file regions.txt --output-type z --output-file shapeit.str.regions.vcf.gz
vcf-sort -c str-snp-synced-reorder.vcf.gz > str-snp-synced-reorder-v2.vcf
cat str-snp-synced-reorder-v2.vcf | sed 's/^chr//g' > str-snp-synced-reorder-v3.vcf
bcftools index shapeit.str.regions.vcf.gz
bgzip str-snp-synced-reorder-v3.vcf
bcftools index str-snp-synced-reorder-v3.vcf.gz 
bcftools concat -a shapeit.str.regions.vcf.gz str-snp-synced-reorder-v3.vcf.gz --output-type z --output $2

