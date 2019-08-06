from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib import rcParams as rc
#rc.update({'font.size': 12})
#import seaborn as sns
from uproot import open as upopen
import numpy as np

file = upopen("../root_files/file_after_weights_3inner.root")
#trees = file.allkeys(filterclass=lambda x: issubclass(x, up.tree.TTreeMethods))
tree = file["data"]
geneta = tree.array('geneta')
showerid = [tree.array('shower_number')[:,0],tree.array('shower_number')[:,1],tree.array('shower_number')[:,2]]

fig, ax = plt.subplots(1, 3, figsize=(15,7), 
                       sharex=False, sharey=False, squeeze=False)
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
for ireg in range(3):
    frac = []
    counts, bins = np.histogram(geneta, bins=100)
    for ib in range(len(bins[:-1])):
        tmp = [ sid for geta,sid in zip(geneta,showerid[ireg]) if geta>bins[ib] and geta<bins[ib+1] ]
        frac.append( 1. - (np.count_nonzero(np.array(tmp)) / float(len(tmp))) )
    print(frac)
    ax[0,ireg].hist(bins[:-1], bins, weights=frac)
plt.xlabel('$|\eta|$')
plt.ylabel('Fraction of complete showers')
plt.savefig('figs/complete_fraction.png')
plt.savefig('/eos/user/b/bfontana/www/ResolutionStudies/complete_fraction.png')
