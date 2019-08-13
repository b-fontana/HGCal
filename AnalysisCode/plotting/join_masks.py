import sys
import numpy as np
from matplotlib import pyplot as plt
nreg = 3
method = sys.argv[1]
samples = sys.argv[2]
width, height = 3, 2
assert(width==nreg)

fig, ax = plt.subplots(height, width, figsize=(15,7), 
                       sharex=False, sharey=False, squeeze=False)
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', 
                top=False, bottom=False, left=False, right=False)
plt.grid(False)
for iax1 in range(height):
    for iax2 in range(width):
            for imask in range(3,7):
                with np.load('numpy_files/arrays_reg'+str(iax2+1)+'_'+str(imask)+samples+'_'+method+'.npz') as f:
                    if iax1==0:
                        ax[iax1,iax2].errorbar(f['xcenters'], f['bias'], yerr=f['sigmas'],
                                               linestyle='', markersize=3, capsize=3,
                                               label='mask'+str(imask))
                        ax[iax1,iax2].set_ylim([-.5,.5])
                    elif iax1==1:
                        ax[iax1,iax2].errorbar(f['xcenters'], f['indep'], 
                                               yerr=f['indep_error'],
                                               linestyle='', markersize=3, capsize=3,
                                               label='mask'+str(imask))
                        ax[iax1,iax2].set_ylim([0.,6.])
                    ax[iax1,iax2].legend()
plt.xlabel(r'|$\eta$|')
plt.savefig('figs/join_' + samples + '_fineeta.png')
plt.savefig( ('/eos/user/b/bfontana/www/ResolutionStudies' 
              + '/join_' + samples + '_' + method + '.png') )
