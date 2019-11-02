import sys
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib import rcParams as rc
#rc.update({'font.size': 12})
#import seaborn as sns
from uproot import open as upopen
import numpy as np

nreg = 3
mask = sys.argv[1]
samples = sys.argv[2]
file = upopen('root_files/final_'+str(mask)+samples+'_ed.root')
df = file["data"].pandas.df(['abs_geneta', 'sid'], flatten=True).unstack()
geneta = df.abs_geneta[2]
showerid = df.sid

fig, ax = plt.subplots(1, nreg, figsize=(15,7), sharex=False, sharey=True, squeeze=False)
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
for ireg in range(nreg):
    frac = [[] for _ in range(4)]
    counts, bins = np.histogram(geneta, bins=100)
    for ib in range(len(bins[:-1])):
        tmp = np.array([ sid for geta,sid in zip(geneta,showerid[ireg]) if geta>bins[ib] and geta<bins[ib+1] ])
        frac[0].append( 1. - float(np.count_nonzero(tmp)) / tmp.size)
        frac[1].append( float(tmp[tmp==1].size) / tmp.size )
        frac[2].append( float(tmp[tmp==2].size) / tmp.size )
        frac[3].append( float(tmp[tmp==3].size) / tmp.size )
    labels = ['signal', 'bckg 1', 'bckg 2', 'bckg 3']
    for iw in range(4):
        ax[0,ireg].hist(bins[:-1], bins, weights=frac[iw], histtype='step', label=labels[iw])
    if ireg==0: ax[0,ireg].legend(loc='lower left')

ax[0,1].set_xlabel('$|\eta|$', fontsize=str(12), labelpad=1.7)
ax[0,0].set_ylabel('Fraction of showers in each category', fontsize=str(12), labelpad=1.7)
radius = dict({'1':'1.3cm', '2':'2.6cm', '3':'5.3cm'})
for ireg in range(nreg):
    plt.text(.05, 1.015, 'Integration cylinder radius: '+radius[str(ireg+1)], transform=ax[0,ireg].transAxes)
    if samples=='outer': ax[0,ireg].set_xticks([1.47, 1.52, 1.57, 1.62])
plt.subplots_adjust(wspace=0.0, hspace=0.)
plt.savefig('/eos/user/b/bfontana/www/PartialWafers/'+samples+'/mask'+str(mask)+'/complete_fraction_'+str(mask)+samples+'.png')
