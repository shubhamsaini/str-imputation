from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import pylab
import matplotlib

header = ["sam","fam"]+["pca"+str(i) for i in range(10)]
pca_eigvec = pd.read_csv("pca_10.eigenvec", sep=" ",names=header)

label = [0]*70 + [1]*50 + [2]*50 + [3]*50 + [4]*50 + [5]*160

colors = ['brown','black','green','blue','purple', 'red']
pop = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS', 'SSC']

eigvec = pca_eigvec.values[:,2:12]

fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.scatter(eigvec[:,0], eigvec[:,1], c=color, cmap=pylab.cm.cool)
# ax.set_xlabel("PC 1", size=15)
# ax.set_ylabel("PC 2", size=15)
# ax.set_xticklabels(ax.get_xticks(), size=12)
# ax.set_yticklabels(ax.get_yticks(), size=12)
# ax.spines["top"].set_visible(False);
# ax.spines["right"].set_visible(False);
# ax.get_xaxis().tick_bottom();
# ax.get_yaxis().tick_left();


plt.scatter(eigvec[:,0], eigvec[:,1], c=label, cmap=matplotlib.colors.ListedColormap(colors))

cb = plt.colorbar()
loc = np.arange(0,max(label),max(label)/float(len(colors)))
cb.set_ticks(loc)
cb.set_ticklabels(pop)
fig.savefig("pca.pdf")
