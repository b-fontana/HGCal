import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import params
sns.set(style='ticks', font_scale=3)
sns.despine()

def diag_ticks(ax1,ax2):
    d = .005  # how big to make the diagonal lines in axes coordinates
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

method, samples = sys.argv[1], sys.argv[2]
df = []
bias, biaserr, ratio_bad, nsamples = ([] for _ in range(4))
etavalues = params.etavalues(samples)
etalabels = ['['+str(np.round(x,3))+';'+str(np.round(y,3))+'[' for x,y in zip(etavalues[:-1],np.roll(etavalues,-1)[:-1])]
for imask in range(3,7):
    with np.load('numpy_files/arrays_reg3_'+str(imask)+samples+'_'+method+'.npz') as f:
        df.append(pd.DataFrame(f['response_eta']).transpose())
        df[-1].columns = ['response', 'eta']
        df[-1] = df[-1][ (df[-1].response>=-1) & (df[-1].response<2.5) ]
        df[-1]['eta_cat'] = pd.cut(df[-1].eta, bins=etavalues, right=False, labels=[i for i in range(len(etavalues)-1)])
        df[-1]['mask'] = imask
        bias.append(f['bias'])        
        biaserr.append(f['biaserr'])        
        ratio_bad.append(f['ratio_bad'])
        nsamples.append(f['nsamples'])

#join the datasets of the fours geometries (masks)
df_tot = pd.concat([df[0],df[1],df[2],df[3]], axis=0, sort=False)

fig, ax = plt.subplots(4, 1, figsize=(60,25), sharex=True,
                       gridspec_kw={'height_ratios': [7,.3,.7,1]})

ax[0].plot([-0.5, 40], [0, 0], linewidth=2, color='grey', linestyle='dashed')
ax[2].plot([-0.5, 40], [0, 0], linewidth=2, color='grey', linestyle='dashed')
if samples == 'outer':
    ax[2].plot([-0.5, 40], [0.05, 0.05], linewidth=2, color='grey', linestyle='dashed')

#boxplot
flierprops = dict(markerfacecolor='0.75', markersize=8,
                  linestyle='none')
axsns = sns.boxplot(x='eta_cat', y='response', hue='mask', data=df_tot, whis=1.5, flierprops=flierprops, ax=ax[0])
for patch in axsns.artists:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r,g,b,.5))
plt.setp(ax[0].get_legend().get_texts(), fontsize='40')
plt.setp(ax[0].get_legend().get_title(), fontsize='40')

lin = (-0.3, len(etavalues)-2.3, len(etavalues)-1)
colors = ['blue', 'orange', 'green', 'red']
for imask in range(4):
    ax[0].set_xticklabels([])
    ax[1].set_xticklabels(etalabels)
    incr = imask*0.2
    imask2 = imask+3
    ax[0].errorbar(np.linspace(lin[0]+incr, lin[1]+incr, lin[2]), 
                   bias[imask], yerr=biaserr[imask], 
                   color=colors[imask], linestyle='', 
                   marker='o', markersize=12, capsize=8, label='fit '+str(imask2))
    ax[0].legend()
    ax[1].plot(np.linspace(lin[0]+incr,lin[1]+incr,lin[2]), 
             ratio_bad[imask], 
             color=colors[imask], linestyle='', 
             marker='s', markersize=14)
    ax[2].plot(np.linspace(lin[0]+incr,lin[1]+incr,lin[2]), 
             ratio_bad[imask], 
             color=colors[imask], linestyle='', 
             marker='s', markersize=14)

    ax[3].plot(np.linspace(lin[0]+incr,lin[1]+incr,lin[2]), 
               nsamples[imask], 
               color=colors[imask], linestyle='', 
               marker='s', markersize=14)

ax[0].set_ylabel(r'Response: $(E_{reco}-E_{gen})/E_{gen}$',
                 fontsize=38)
ax[0].set_ylim([-1.1,2.5])

ax[1].set_ylabel(r'$\frac{N_{\Delta E/E<-0.95}}{N_{total}}$',
                 fontsize=42, labelpad=50)
#ax[1].set_ylim([-0.05,0.2] if samples=='inner' else [-0.15,0.85])
ax[1].set_yticks([1.] if samples=='inner' else [0.6])
ax[1].grid()
ax[1].set_ylim(0.95,1.05) if samples == 'inner' else  ax[1].set_ylim(0.45,0.75)
ax[2].set_ylim(-0.2,0.4) if samples == 'inner' else  ax[2].set_ylim(-0.05,0.2)
ax[2].set_yticks([0.,0.2] if samples=='inner' else [0.,0.15])
ax[1].spines['bottom'].set_visible(False)
ax[2].spines['top'].set_visible(False)
ax[1].xaxis.tick_top()
ax[1].tick_params(labeltop='off')  # don't put tick labels at the top
ax[2].xaxis.tick_bottom()
diag_ticks(ax[1],ax[2])

ax[3].set_ylabel(r'$N_{total}$',
                 fontsize=35)
ax[3].set_xlabel(r'$|\eta|$', fontsize=38)
ax[3].grid()
#ax[3].set_yscale('log')
ax[3].set_yticks([350,3500] if samples=='inner' else [200,20000])

plt.subplots_adjust(wspace=0., hspace=0.)
plt.savefig('final_'+samples+'.png')
plt.savefig('/eos/user/b/bfontana/www/ResolutionStudies/final_'+samples+'.png')
