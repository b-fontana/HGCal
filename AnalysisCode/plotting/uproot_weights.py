import sys
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib import rcParams as rc
#rc.update({'font.size': 12})
#import seaborn as sns
from uproot import open as upopen
import numpy as np

mask = sys.argv[1]
samples = sys.argv[2]
file = upopen("root_files/fileweights_"+str(mask)+samples+".root")
keys = file.keys()

nreg = 3
figsize = (15.,12.)
figscale = figsize[0] / figsize[1]
figratios = [1,.5]
fig, ax = plt.subplots(len(figratios), nreg, figsize=figsize, 
                       sharex=True, gridspec_kw={'height_ratios': figratios}, 
                       sharey=False, squeeze=False)
fig.add_subplot(111, frameon=False)
plt.subplots_adjust(wspace=0.2, hspace=0.)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
colors = ['blue', 'orange','green','brown']
#energy distributions
for ireg in range(nreg):
    counts, bins = file['en{}_layer_sign;1'.format(ireg+1)].numpy()
    ax[0,ireg].hist(bins[:-1], bins, weights=counts, histtype='step', label='signal', color=colors[0])
    for iw in range(3):
        counts, bins = file['en{}_layers_bckg{};1'.format(ireg+1,iw+1)].numpy()
        ax[0,ireg].hist(bins[:-1], bins, weights=counts, histtype='step', label='bckg'+str(iw+1), color=colors[iw+1])
    ax[0,ireg].legend()
#weights
for ireg in range(nreg):
    for iw in range(3):
        counts, bins = file['weight{}_sr{};1'.format(iw+1,ireg+1)].numpy()
        ax[1,ireg].hist(bins[:-1], bins, weights=counts, histtype='step', label='weight'+str(iw+1), color=colors[iw+1])
        ax[1,ireg].legend(loc='upper right')
radius = dict({'1':'1.3cm', '2':'2.6cm', '3':'5.3cm'})
for ireg in range(nreg):
    plt.text(.05, 1.015, 'Integration cylinder radius: '+radius[str(ireg+1)], transform=ax[0,ireg].transAxes)
    ax[0,ireg].set_yticks([0.025, 0.05, 0.075, 0.1, 0.125])
ax[1,1].set_xlabel('Layers', fontsize=str(12*figscale), labelpad=1.7*figscale)
ax[0,0].set_ylabel('Fraction of deposited energy', fontsize=str(12*figscale), labelpad=1.7*figscale)
ax[1,0].set_ylabel('Weights', fontsize=str(12*figscale), labelpad=1.7*figscale)
plt.savefig('/eos/user/b/bfontana/www/PartialWafers/'+samples+'/mask'+str(mask)+'/weights_'+str(mask)+samples+'.png')
