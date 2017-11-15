######
### usage: python corr.py "STR_868537" str.txt snp.txt snp.id.txt output.txt
######


import numpy as np
import pandas as pd
import sys

str_id = sys.argv[1]

str = pd.read_csv(sys.argv[2], names=["GT"])
str = np.array([len(i) for i in str.GT.values])
snp = np.loadtxt(sys.argv[3])

snp_sum = snp.reshape(-1, int(snp.shape[0]/2), 2).sum(axis=2).transpose()
str_sum = str.reshape(-1, int(str.shape[0]/2), 2).sum(axis=2).transpose()

output = list()
for i in range(snp_sum.shape[1]):
    if np.std(str_sum.flatten()) == 0:
        output.append(-9)
    elif np.std(snp_sum[:,i]) == 0:
        output.append(-99)
    else:
        output.append(np.corrcoef(str_sum.flatten(),snp_sum[:,i])[0,1])

header = ["ID"]
snp_id = pd.read_csv(sys.argv[4], names=header)
max_id = snp_id.ID[np.argmax(output)]
max_r = np.max(output)

write_to_file = " ".join([str_id, max_id, max_r.astype("str"), "\n"])

with open(sys.argv[5], "a") as myfile:
    myfile.write(write_to_file)
