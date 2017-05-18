#!/bin/sh

ls *.bam > files.list
cat files.list | xargs -I% -n1 samtools sort % %
mkdir sorted
mv *.bam.bam sorted/
rm *.bam
mv sorted/*.bam .
cat files.list | xargs -I% -n1 mv %.bam %
cat files.list | xargs -I% -n1 samtools index %