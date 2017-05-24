#!/bin/sh

ls *.bam > files.list
cat files.list | xargs -I% -P10 -n1 samtools sort % -o %
cat files.list | xargs -I% -P10 -n1 samtools index %