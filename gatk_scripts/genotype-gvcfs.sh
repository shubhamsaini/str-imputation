#!/bin/bash
#SBATCH -A ddp268
#SBATCH -p shared
#SBATCH --get-user-env
#SBATCH -t 1-00:00 # Walltime Days-Hours-Minutes
#SBATCH --job-name=chr5part6
#SBATCH -o /home/s1saini/log/chr5part6_out
#SBATCH -e /home/s1saini/log/chr5part6_err
#SBATCH --mem=16G
java -jar /home/s1saini/gatk.jar -T GenotypeGVCFs -V /home/s1saini/chr5.list -R /oasis/scratch/comet/s1saini/temp_project/ref/human_g1k_v37.fasta -newQual -L 5:25000001-30000000 -o /oasis/scratch/comet/s1saini/temp_project/chr5.part6.vcf.gz