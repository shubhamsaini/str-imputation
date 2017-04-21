#!/bin/bash

### usage: 
### ./runbeagle.sh -s sampleID -p pedigree.txt -f samples.txt -l final.str.snp.vcf.gz -t hipstrcalls.vcf.gz -b path/to/beagle.exe
### Eg: ./runbeagle.sh -s SSC00092 -p pedigree.txt -f sampleID.txt -l final.str.snp.vcf.gz -t str-vcf-paired-v2.vcf.gz -b beagle.27Jul16.86a.jar 
###

while getopts s:p:f:t:b:l: option
do
        case "${option}"
        in
                s) sss=${OPTARG};;
                p) ped=${OPTARG};;
                f) samFile=$OPTARG;;
                t) str=$OPTARG;;
                b) beagle=$OPTARG;;
                l) procStr=$OPTARG;;
        esac
done

export SAMPLEID=$sss
echo $SAMPLEID
export FATHERID=`cat $ped | grep $SAMPLEID | awk '{print $2}' | head -1`
export MOTHERID=`cat $ped | grep $SAMPLEID | awk '{print $3}' | head -1`
export CHILD1=`cat $ped | grep $FATHERID | awk '{print $1}' | head -1`
export CHILD2=`cat $ped | grep $FATHERID | awk '{print $1}' | tail -n 1`

rm sampleExclude.txt
echo $FATHERID > sampleExclude.txt
echo $MOTHERID >> sampleExclude.txt
echo $CHILD1 >> sampleExclude.txt
echo $CHILD2 >> sampleExclude.txt

cat $samFile | grep -v $FATHERID | grep -v $MOTHERID | grep -v $CHILD1 | grep -v $CHILD2 > sampleRef.txt

bcftools view $procStr --samples-file sampleRef.txt --output-type z --output-file ref.vcf.gz
bcftools view $procStr --samples $SAMPLEID -f PASS --output-type z --output-file exclude.vcf.gz

zcat < ref.vcf.gz | awk '!x[$0]++' > ref.vcf
rm ref.vcf.gz
bgzip ref.vcf

java  -Xmx8g -jar $beagle gt=ref.vcf.gz out=ref.imputed
java  -Xmx8g -jar $beagle gt=exclude.vcf.gz ref=ref.imputed.vcf.gz out=imputed
bcftools index imputed.vcf.gz 
bcftools view imputed.vcf.gz -m 3 > imputed.str

bcftools query -f '%CHROM:%POS\t[%GT\t]\n' imputed.str | gsort -V | sed 's/^chr//g' > imputeResult.txt
bcftools query --samples $SAMPLEID -f '%CHROM:%POS\t[%GT\t]\n' $str | gsort -V | sed 's/^chr//g' | grep -v '17:' > groundTruth.txt

sdiff imputeResult.txt groundTruth.txt 
echo "done"




