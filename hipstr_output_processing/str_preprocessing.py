

## python str_processing.py hipstr.output.vcf.gz pedigree.fam shapeit.phased.snps.vcf.gz 40000


from cyvcf2 import Variant, VCF, Writer
import pandas as pd
import numpy as np
from collections import defaultdict
import sys
import os

header = ['fam','child','father','mother','sex','pheno']
pedigree = pd.read_csv(sys.argv[2],names=header,delim_whitespace=True)
pedigree = pedigree[pedigree.sex==0]
pedigree = pedigree[['child','father','mother']]

vcf = VCF(sys.argv[1])
sampleID = dict(zip(range(len(vcf.samples)),vcf.samples))

for i in range(len(sampleID)):
    pedigree[pedigree==sampleID[i]] = i
    
vcf = VCF(sys.argv[1])
w = Writer("temp.partialphased.vcf", vcf)
for v in vcf():
    gtlist = v.genotypes
    for idx in range(pedigree.values.shape[0]):
        childID = pedigree.values[idx,0]
        fatherID = pedigree.values[idx,1]
        motherID = pedigree.values[idx,2]
        #if motherID == 54:
            #print v.CHROM, v.POS, gtlist[motherID]
        #print childID,fatherID,motherID
        if (v.genotypes[childID][0] in v.genotypes[fatherID][0:2]) and (v.genotypes[childID][0] not in v.genotypes[motherID][0:2]):
            print ""
            #print "F|M"
        elif (v.genotypes[childID][0] in v.genotypes[motherID][0:2]) and (v.genotypes[childID][0] not in v.genotypes[fatherID][0:2]):
            if len(gtlist[childID])>2:
                print ""
                #print "M|F"
                allA = gtlist[childID][0]
                gtlist[childID][0] = gtlist[childID][1]
                gtlist[childID][1] = allA
        elif (v.genotypes[childID][1] in v.genotypes[motherID][0:2]) and (v.genotypes[childID][1] not in v.genotypes[fatherID][0:2]):
            print ""
            #print "F|M"
        elif (v.genotypes[childID][1] in v.genotypes[fatherID][0:2]) and (v.genotypes[childID][1] not in v.genotypes[motherID][0:2]):
            if len(gtlist[childID])>2:
                print ""
                #print "M|F"
                allA = gtlist[childID][0]
                gtlist[childID][0] = gtlist[childID][1]
                gtlist[childID][1] = allA
        else:
            if v.format('PQ')[childID][0] < 0.8:
                #gtlist[childID][len(gtlist[childID])-1] = False
                print ""
                #print "unphased"
            else:
                print ""
                #print "unchanged"
        
        if len(gtlist[childID]) < 3:
            gtlist[childID][-1] = False
        if len(gtlist[fatherID]) < 3:
            gtlist[fatherID][-1] = False
        if len(gtlist[motherID]) < 3:
            gtlist[motherID][-1] = False
    v.genotypes = gtlist
    w.write_record(v)   
w.close()

header = ['fam','child','father','mother','sex','pheno']
pedigree = pd.read_csv(sys.argv[2],names=header,delim_whitespace=True)
pedigree = pedigree[pedigree.sex==0]
pedigree = pedigree[['child','father','mother']]

parentToChild = defaultdict(list)

fathers = list(np.unique(pedigree.transpose().values[1]))
mothers = list(np.unique(pedigree.transpose().values[2]))
for i in fathers:
    parentToChild[i] = parentToChild[i] + list(pedigree[pedigree.father==i]['child'].values)
for i in mothers:
    parentToChild[i] = parentToChild[i] + list(pedigree[pedigree.mother==i]['child'].values)
    
    
buff = int(sys.argv[4])
strvcf = VCF('temp.partialphased.vcf')
w = Writer("str-snp-synced.vcf", strvcf)

sampleToIndex = dict(zip(strvcf.samples,range(len(strvcf.samples))))
for v in strvcf():
    pos = v.POS
    chrom = int(filter(str.isdigit, str(v.CHROM)))
    start = pos - (buff/2)
    end = pos + (buff/2)
    region = "".join([str(chrom),":",str(start),"-",str(end)])
    gtlist = v.genotypes
    for i in parentToChild:
        child1ID = sampleToIndex[parentToChild[i][0]]
        child2ID = sampleToIndex[parentToChild[i][1]]
        parentID = sampleToIndex[i]
        child1STR = v.genotypes[child1ID]
        child2STR = v.genotypes[child2ID]
        parentSTR = v.genotypes[parentID]
        if child1STR[-1]:
            childID = child1ID
            childSample = parentToChild[i][0]
        elif child2STR[-1]:
            childID = child2ID
            childSample = parentToChild[i][1]
        else:
            continue
        isFather = i in fathers
        childHap1 = list()
        parentHap1 = list()
        childHap2 = list()
        parentHap2 = list()
        snpvcf = VCF(sys.argv[3],samples=[childSample,i])
        for record in snpvcf(region):
            genotypeData = record.genotypes
            if len(genotypeData[0]) > 2 and len(genotypeData[1]) > 2:
                childHap1.append(genotypeData[0][0])
                childHap2.append(genotypeData[0][1])
                parentHap1.append(genotypeData[1][0])
                parentHap2.append(genotypeData[1][1])
        if isFather:
            if np.sum(np.array(childHap1) == np.array(parentHap1)) > np.sum(np.array(childHap1) == np.array(parentHap2)):
                if v.genotypes[parentID][1] == v.genotypes[childID][0]:
                    allA = gtlist[parentID][0]
                    gtlist[parentID][0] = gtlist[parentID][1]
                    gtlist[parentID][1] = allA
                    print "case 1a"
                    #print pos
                    break
            else:
                if v.genotypes[parentID][0] == v.genotypes[childID][0]:
                    allA = gtlist[parentID][0]
                    gtlist[parentID][0] = gtlist[parentID][1]
                    gtlist[parentID][1] = allA
                    print "case 1b"
        else:
            if np.sum(np.array(childHap2) == np.array(parentHap2)) > np.sum(np.array(childHap2) == np.array(parentHap1)):
                if v.genotypes[parentID][0] == v.genotypes[childID][1]:
                    allA = gtlist[parentID][0]
                    gtlist[parentID][0] = gtlist[parentID][1]
                    gtlist[parentID][1] = allA
                    print "case 2a"
            else:
                if v.genotypes[parentID][1] == v.genotypes[childID][1]:
                    allA = gtlist[parentID][0]
                    gtlist[parentID][0] = gtlist[parentID][1]
                    gtlist[parentID][1] = allA
                    print "case 2b"
    v.genotypes = gtlist
    w.write_record(v)
    
w.close()

delCmd = "rm regions.txt"
os.popen(delCmd).read()

strvcf = VCF('temp.partialphased.vcf')
for v in strvcf():
    pos = v.POS
    chrom = int(filter(str.isdigit, str(v.CHROM)))
    start = pos - (buff/2)
    end = pos + (buff/2)
    data = str(chrom) + "\t" + str(start) + "\t" + str(end)
    with open("regions.txt", "a") as myfile:
        myfile.write(data+"\n")