#!/bin/bash

### Script to analyse the L1O analysis results
### Outputs a file *.csv

### Edit the python file before using this

ls *ground*txt | sed 's/.groundTruth.txt//' > samplesID.txt
cat samplesID.txt | xargs -I% sh -c "sed 's/\./NA\/NA/'  %.groundTruth.txt  | awk -F\"/\" '\$1=\$1' OFS=\"\t\" >  %.groundTruth.sum"
cat samplesID.txt | xargs -I% sh -c "sed 's/\./NA\/NA/'  %.imputeResult.txt  | awk -F\"|\" '\$1=\$1' OFS=\"\t\" >  %.imputeResult.sum"

