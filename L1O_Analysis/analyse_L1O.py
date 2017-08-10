import pandas as pd
import numpy as np
from scipy.stats.stats import pearsonr   


sampleHead = ["sampleID"]
samples = pd.read_csv("results/LOO_0810/samplesID.txt",names=sampleHead)

headerIm = ['str','imal1','imal2']
headerGr = ['str','gral1','gral2']
imputed = pd.read_csv("results/LOO_0810//"+samples.values[0][0]+".imputeResult.sum",names=headerIm,delim_whitespace=True)
ground = pd.read_csv("results/LOO_0810//"+samples.values[0][0]+".groundTruth.sum",names=headerGr,delim_whitespace=True)
mergeData = imputed.merge(ground, how="inner", on="str")

for i in samples.values[1:]:
    imputed = pd.read_csv("results/LOO_0810//"+i[0]+".imputeResult.sum",names=headerIm,delim_whitespace=True)
    ground = pd.read_csv("results/LOO_0810//"+i[0]+".groundTruth.sum",names=headerGr,delim_whitespace=True)
    data = imputed.merge(ground, how="inner", on="str")
    mergeData = mergeData.append(data)
    
sumMerged = pd.DataFrame({'str':mergeData['str'], 'im':(mergeData['imal1'].values + mergeData['imal2'].values), 'gr':(mergeData['gral1'].values + mergeData['gral2'].values)})
corr = sumMerged.groupby('str')[['im','gr']].corr().ix[0::2,'gr']
corr = corr.reset_index()[['str','gr']]
corr.columns = ['str','r2']

droppedNa = mergeData.dropna(axis=0)
concord = list()
for i in droppedNa.values:
    listA = set(i[1:3].astype(int))
    listB = set(i[3:5].astype(int))
    concord.append( (2-(max(len(listA-listB) , len(listB-listA))))/2.0 )
    
concordance = pd.DataFrame({'str':droppedNa['str'], 'concord':concord})
concordance = concordance.groupby('str').mean().reset_index()

bpHeader= ['str','bpdiff']
numAllele = list()
bpdiff = pd.read_csv("results/BPDIFF.txt",names=bpHeader,delim_whitespace=True)
for i in bpdiff.values:
    n = len(i[1].split(","))
    numAllele.append(n)

bpdiff = pd.DataFrame({'str':bpdiff['str'], 'numAllele':numAllele})

prHeader= ['str','motif_len']
period = pd.read_csv("results/PERIOD.txt",names=prHeader,delim_whitespace=True)

finalData = corr.merge(concordance, how="inner", on="str")
finalData = finalData.merge(bpdiff,how="inner",on="str")
finalData = finalData.merge(period,how="inner",on="str")


finalData.to_csv("chr21.stats.csv")