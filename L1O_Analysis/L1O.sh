#!/bin/bash

### usage: 
### ./L1O.sh -s sampleID -p pedigree.fam -l final.str.snp.vcf.gz -t hipstrcalls.vcf.gz -b path/to/beagle
### Eg: ./L1O.sh -s SSC00092 -p pedigree.fam -l final.str.snp.vcf.gz -t str-vcf-paired-v2.vcf.gz -b beagle.27Jul16.86a.jar 
###

while getopts s:p:f:t:b:l: option
do
        case "${option}"
        in
                s) sss=${OPTARG};;
                p) ped=${OPTARG};;
                t) str=$OPTARG;;
                b) beagle=$OPTARG;;
                l) procStr=$OPTARG;;
        esac
done

export SAMPLEID=$sss
echo $SAMPLEID
export famID=`cat $ped | grep $SAMPLEID | head -1 | awk '{print $1}'`
export CHILD1=`cat $ped | awk -v var="$famID" '$1==var {print $0}' | awk '$3!=0 {print $0}' | head -1 | awk '{print $2}'`
export CHILD2=`cat $ped | awk -v var="$famID" '$1==var {print $0}' | awk '$3!=0 {print $0}' | tail -n 1 | awk '{print $2}'`
export FATHERID=`cat $ped | awk -v var="$famID" '$1==var {print $0}' | awk '$3!=0 {print $0}' | head -1 | awk '{print $3}'`
export MOTHERID=`cat $ped | awk -v var="$famID" '$1==var {print $0}' | awk '$3!=0 {print $0}' | head -1 | awk '{print $4}'`

cat $ped | awk '$3==0 {print $2}' | grep -v $FATHERID | grep -v $MOTHERID | grep -v $CHILD1 | grep -v $CHILD2 > sampleRef.txt

bcftools query -f '%ID\n' $str > ID.txt
bcftools view $procStr --samples-file sampleRef.txt --output-type z --output-file ref.vcf.gz --force-samples
bcftools view $procStr --samples $SAMPLEID --exclude ID=@ID.txt --output-type z --output-file exclude.vcf.gz --force-samples

bcftools norm -d any ref.vcf.gz -O z -o ref.rem.vcf.gz
rm ref.vcf.gz
mv ref.rem.vcf.gz ref.vcf.gz
bcftools index -f ref.vcf.gz

java  -Xmx8g -jar $beagle gt=exclude.vcf.gz ref=ref.vcf.gz out=imputed
bcftools index imputed.vcf.gz 
bcftools view imputed.vcf.gz --include ID=@ID.txt > imputed.str

python write_bases.py imputed.str $SAMPLEID $SAMPLEID.imputeResult.txt
python write_bases.py $str $SAMPLEID $SAMPLEID.groundTruth.txt

echo "done"

rm exclude.vcf.gz ref.vcf.gz sampleExclude.txt sampleRef.txt ref.vcf.gz.csi imputed.vcf.gz imputed.vcf.gz.csi imputed.str
