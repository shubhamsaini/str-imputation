from cyvcf2 import Variant, VCF, Writer
import pandas as pd
import numpy as np
from collections import defaultdict
import sys
import os


vcf = VCF(sys.argv[1])
w = Writer(sys.argv[1]+"unphased", vcf)
for v in vcf():
    gtlist = v.genotypes
    for idx in range(len(vcf.samples)):
        gtlist[idx][len(gtlist[idx])-1] = False
    v.genotypes = gtlist
    w.write_record(v)   
w.close()