from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib import rcParams as rc
#rc.update({'font.size': 12})
#import seaborn as sns
from uproot import open as upopen
import numpy as np

file = upopen("../root_files/fileweights_3inner.root")
keys = file.keys()

nreg = 3
fig, ax = plt.subplots(2, nreg, figsize=(15,14), 
                       sharex=False, sharey=False, squeeze=False)
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
#energy distributions
for ireg in range(nreg):
    counts, bins = file['en{}_layer_sign;1'.format(ireg+1)].numpy()
    ax[0,ireg].hist(bins[:-1], bins, weights=counts, histtype='step', label='signal')
    for iw in range(3):
        counts, bins = file['en{}_layers_bckg{};1'.format(ireg+1,iw+1)].numpy()
        ax[0,ireg].hist(bins[:-1], bins, weights=counts, histtype='step', label='bckg'+str(iw+1))
    ax[0,ireg].legend()
#weights
for ireg in range(nreg):
    for iw in range(3):
        counts, bins = file['weight{}_sr{};1'.format(iw+1,ireg+1)].numpy()
        ax[1,ireg].hist(bins[:-1], bins, weights=counts, histtype='step', label='weight'+str(iw+1))
        ax[1,ireg].legend()
plt.xlabel('Layers')
plt.savefig('figs/weights.png')
plt.savefig('/eos/user/b/bfontana/www/ResolutionStudies/weights.png')
