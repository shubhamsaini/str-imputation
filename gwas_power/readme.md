## Finding SNP with highest LD with STRs
Usage

    ./driver.sh chromosome str.vcf.gz snp.vcf.gz numThreads window

Example:

    ./driver.sh 21 hipstr_calls_chr21.vcf.gz ssc.shapeit.chr21.vcf.gz 4 100000

Note:

 - This script does not take genotype phase into consideration
 - It is assumed that samples are in the same order in both STR and SNP VCFs.
	 - If not, uncomment the three lines to reorder the samples in the SNP VCF