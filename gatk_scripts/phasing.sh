#!/bin/bash
#SBATCH -A ddp268
#SBATCH -p shared
#SBATCH --get-user-env
#SBATCH -t 0-12:00 # Walltime Days-Hours-Minutes
#SBATCH --job-name=chr3reg12
#SBATCH -o /home/s1saini/shapeitlog/chr3reg12.out
#SBATCH -e /home/s1saini/shapeitlog/chr3reg12.err
#SBATCH --mem=16G
#SBATCH -c 16
/home/s1saini/shapeit/bin/shapeit \
-B /oasis/scratch/comet/s1saini/temp_project/chr3/shapeit/chr3.reg12\_me \
-M /oasis/projects/nsf/ddp268/s1saini/genetic_map_b37/genetic_map_chr3_combined_b37.txt \
--duohmm -W 5 \
--thread 16 \
-O /oasis/scratch/comet/s1saini/temp_project/chr3/shapeit/shapeit.chr3.reg12 
/home/s1saini/shapeit/bin/shapeit \
-convert \
--input-haps /oasis/scratch/comet/s1saini/temp_project/chr3/shapeit/shapeit.chr3.reg12 \
--output-vcf /oasis/scratch/comet/s1saini/temp_project/chr3/shapeit/shapeit.chr3.reg12.vcf 