import pandas as pd
import numpy as np
from scipy.stats.stats import pearsonr   
import sys

dir = sys.argv[1]

sampleHead = ["sampleID"]
samples = pd.read_csv(dir+"/"+"samplesID.txt",names=sampleHead)

header = ['str','gral1','gral2','imal1','imal2']
mergeData = pd.read_csv(dir+"/"+samples.values[0][0]+".diff.txt",names=header,delim_whitespace=True)

for i in samples.values[1:]:
    data = pd.read_csv(dir+"/"+i[0]+".diff.txt",names=header,delim_whitespace=True)
    mergeData = mergeData.append(data)


sumMerged = pd.DataFrame({'str':mergeData['str'], 'im':(mergeData['imal1'].values + mergeData['imal2'].values), 'gr':(mergeData['gral1'].values + mergeData['gral2'].values)})

### remove calls with allele counts < 3
counts = sumMerged.groupby(['str','gr']).count().reset_index()
counts.columns = ['str','gr','counts']
sumMerged = sumMerged.merge(counts, how="inner", on=["str","gr"])
sumMerged = sumMerged[sumMerged.counts>=3]

tmp1 = sumMerged.groupby('str')[['gr','im']].corr()
corr = (tmp1[~tmp1['gr'].eq(1)].reset_index(1, drop=True)['gr'].rename('r').reset_index())
corr = corr.dropna()

### find number of samples
numSamples = sumMerged.groupby(['str']).count().reset_index()[['str','counts']]
numSamples.columns = ['str','numSamples']

### calculate p-Value
def pVal(group, col1, col2):
    c1 = group[col1]
    c2 = group[col2]
    return pearsonr(c1, c2)[1]


pVals = sumMerged.groupby('str').apply(pVal, 'gr', 'im').reset_index()
pVals.columns = ['str','pVal']
corr = corr.merge(pVals, how="inner", on="str")

droppedNa = mergeData.dropna(axis=0)
concord = list()
for i in droppedNa.values:
    listA = set(i[1:3].astype(int))
    listB = set(i[3:5].astype(int))
    concord.append( (2-(max(len(listA-listB) , len(listB-listA))))/2.0 )
    
concordance = pd.DataFrame({'str':droppedNa['str'], 'concord':concord})
#concordance['concord'].value_counts()
concordance = concordance.groupby('str').mean().reset_index()

bpHeader= ['str','bpdiff']
numAllele = list()
bpdiff = pd.read_csv(dir+"/"+"BPDIFF.txt",names=bpHeader,delim_whitespace=True)
for i in bpdiff.values:
    n = len(i[1].split(","))
    numAllele.append(n)

bpdiff = pd.DataFrame({'str':bpdiff['str'], 'numAllele':numAllele})

prHeader= ['str','motif_len']
period = pd.read_csv(dir+"/"+"PERIOD.txt",names=prHeader,delim_whitespace=True)

finalData = corr.merge(concordance, how="inner", on="str")
finalData = finalData.merge(bpdiff,how="inner",on="str")
finalData = finalData.merge(period,how="inner",on="str")
finalData = finalData.merge(numSamples,how="inner",on="str")

finalData.to_csv("l1o.results.csv")
