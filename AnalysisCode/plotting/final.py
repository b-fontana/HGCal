import sys, os
import copy
import csv
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt
import seaborn as sns
sns.set(style='ticks', font_scale=3)
import params
from dump_data import dump_tcl_resolution_formula

extra_str = ''
figsize = (60.,30.)
figscale = figsize[0] / figsize[1]
method, samples, region = sys.argv[1], sys.argv[2], sys.argv[3]
df = []
bias, biaserr, bad75, bad50, nsamples, quart_dist = ([] for _ in range(6))
etavalues = params.etavalues(samples)
etameans = [(x+y)/2 for x,y in zip(etavalues[:-1],np.roll(etavalues,-1)[:-1])]
etalabels = ['['+str(np.round(x,3))+';'+str(np.round(y,3))+'[' for x,y in zip(etavalues[:-1],np.roll(etavalues,-1)[:-1])]
radiuslabels = ['\n ' + str(round(322.103*np.tan(2*np.arctan(np.exp(-etameans[i]))),2)) + ' cm' for i in range(len(etameans))]
etalabels = [etalabels[i]+radiuslabels[i] for i in range(len(etalabels))]
method_str = 'corr_fineeta' if method == 'fineeta' else 'corr_ed' if method == 'ed' else 'nocorr'
for imask in range(3,7):
    with np.load('numpy_files/arrays_reg'+str(region)+'_'+str(imask)+samples+'_'+method_str+'.npz') as f:
        df.append(pd.DataFrame(f['response_eta']).transpose())
        df[-1].columns = ['response', 'eta']
        df[-1] = df[-1][ (df[-1].response>=-1) & (df[-1].response<2.5) ]
        df[-1]['eta_cat'] = pd.cut(df[-1].eta, bins=etavalues, right=False, labels=[i for i in range(len(etavalues)-1)])
        df[-1]['mask'] = imask
        bias.append(f['bias']) 
        biaserr.append(f['biaserr'])        
        nsamples.append(df[-1].groupby('eta_cat')['response'].count().to_numpy(dtype=np.float32))
        badratio = lambda x: df[-1][df[-1].response<x].groupby('eta_cat')['response'].count().to_numpy()
        bad75.append(badratio(-.75) / nsamples[-1])
        bad50.append(badratio(-.50) / nsamples[-1])
        quantiles = df[-1].groupby('eta_cat')['response'].quantile([.25,.75]).to_numpy()
        quart_dist.append(np.abs(quantiles[0::2]-quantiles[1::2]))

#join the datasets of the fours geometries (masks)
df_tot = pd.concat([df[0],df[1],df[2],df[3]], axis=0, sort=False)

figratios = [1,.2]
fig, ax = plt.subplots(len(figratios), 1, figsize=figsize, sharex=True,
                       gridspec_kw={'height_ratios': figratios})
axes = dict({'main':ax[0], 'ntot':ax[1]})
axes['main'].plot([-0.5, 40], [0, 0], linewidth=2, color='grey', linestyle='dashed')

#boxplot
flierprops = dict(markerfacecolor='0.75', markersize=5*figscale,
                  linestyle='none')
axsns = sns.boxplot(x='eta_cat', y='response', hue='mask', data=df_tot, whis=1.5, flierprops=flierprops, ax=axes['main'])
for patch in axsns.artists:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r,g,b,.5))
plt.setp(axes['main'].get_legend().get_texts(), fontsize=str(6*figscale))

lin = (-0.3, len(etavalues)-2.3 if samples=='inner' else len(etavalues)-2.3, len(etavalues)-1)
colors = ['blue', 'orange', 'green', 'red']
handles = []
for imask in range(4):
    axes['ntot'].set_xticklabels(etalabels)
    incr = imask*0.2
    imask2 = imask+3
    axes['main'].errorbar(np.linspace(lin[0]+incr, lin[1]+incr, lin[2]), 
                          bias[imask], yerr=biaserr[imask], 
                          color=colors[imask], linestyle='', 
                          marker='o', markersize=7*figscale, capsize=5*figscale, label='fit '+str(imask2))

    options = dict(color=colors[imask], linestyle='', 
                   markersize=7*figscale, markeredgecolor='black', markeredgewidth=3)
    axes['ntot'].plot(np.linspace(lin[0]+incr,lin[1]+incr,lin[2]), 
                      nsamples[imask], marker='s', **options)

axes['main'].errorbar(np.NaN, np.NaN, marker='D', color='grey', label='outliers', ms=13) #add legend entry for outliers
axes['main'].legend(loc='upper right' if samples=='outer' else 'upper left', title='Masks')
axes['main'].set_ylabel(r'Response: $(E_{reco}-E_{gen})/E_{gen}$',
                        fontsize=15*figscale)
axes['main'].set_ylim([-1.1,2.])
axes['ntot'].set_ylabel(r'$N_{total}$', fontsize=15*figscale)
axes['ntot'].set_xlabel(r'$|\eta|$, radius of layer #1', fontsize=15*figscale, labelpad=18*figscale)
axes['ntot'].grid(linewidth=3)
axes['ntot'].set_yticks([350,2000,3000] if samples=='inner' else [200,10000,20000])
axes['ntot'].set_ylim([0,3550] if samples=='inner' else [0,30000])

plt.subplots_adjust(wspace=0., hspace=0.)

texth = 1.015
plt.text(.0, texth, 'CMS Preliminary', ha='left', transform=axes['main'].transAxes)
radius = dict({'1':'1.3cm', '2':'2.6cm', '3':'5.3cm'})
plt.text(.815, texth, 'Integration cylinder radius: '+radius[region], transform=axes['main'].transAxes)
methodmap = dict({'ed':'Shower reconstruction', 'fineeta':'Brute force calibration', 'nocorr':'No Correction'})
plt.text(.34, texth, 'Software correction method: '+methodmap[method], transform=axes['main'].transAxes)
plt.savefig('/eos/user/b/bfontana/www/PartialWafers/'+samples+'/final_reg'+str(region)+'_'+samples+'_'+method+extra_str+'.png')

########################################################################
########################################################################
figratios = [.6,1]
fig, ax = plt.subplots(len(figratios), 1, figsize=figsize, sharex=True,
                       gridspec_kw={'height_ratios': figratios})
axes = dict({'res':ax[0], 'bad':ax[1]})

for imask in range(4):
    axes['bad'].set_xticks([x for x in range(len(etavalues)-1)])
    axes['bad'].set_xticklabels(etalabels)
    m = 12.5 if samples=='inner' else 13.5
    for i in [x for x in np.linspace(1.5,m,len(etavalues)-2)-1]:
        axes['bad'].plot([i, i], [0, 0.55] if samples=='inner' else [0, 0.95], linewidth=2, color='grey', linestyle=(0, (5,10)))
        axes['res'].plot([i, i], [0.05, 1.35] if samples=='inner' else [0.05, 1.56], linewidth=2, color='grey', linestyle=(0, (5,10)))
    incr = imask*0.2
    imask2 = imask+3
    options = dict(color=colors[imask], linestyle='', 
                   markersize=10*figscale, markeredgecolor='black', markeredgewidth=3)
    if(imask==0):
        h, = axes['bad'].plot(np.linspace(lin[0]+incr,lin[1]+incr,lin[2]), 
                              bad75[imask], label='$k=-0.75$', marker='s', **options)
        handles.append(copy.copy(h))
        h, = axes['bad'].plot(np.linspace(lin[0]+incr,lin[1]+incr,lin[2]), 
                              bad50[imask], label='$k=-0.50$', marker='v', **options)
        handles.append(copy.copy(h))
    else:
        axes['bad'].plot(np.linspace(lin[0]+incr,lin[1]+incr,lin[2]), 
                         bad75[imask], marker='s', **options)
        axes['bad'].plot(np.linspace(lin[0]+incr,lin[1]+incr,lin[2]), 
                         bad50[imask], marker='v', **options)

    axes['res'].plot(np.linspace(lin[0]+incr,lin[1]+incr,lin[2]), 
                     quart_dist[imask], marker='s', label=str(imask2), **options)

#manually change the color of the legend
for h in handles:
    h.set_color('black')
axes['bad'].legend(handles=handles, fontsize=14*figscale)
axes['bad'].set_ylabel(r'$\frac{N_{\Delta E/E<k}}{N_{total}}$',
                 fontsize=25*figscale, labelpad=50)
axes['bad'].set_yticks([0.,.1,.2,.3,.4,.5] if samples=='inner' else [0.,0.3,0.6,0.9])
axes['bad'].set_ylim([-0.03,0.6] if samples=='inner' else [-0.03,0.95])
axes['bad'].grid(linewidth=3)
axes['bad'].set_xlabel(r'$|\eta|$, radius of layer #1', fontsize=18*figscale, labelpad=18*figscale)

axes['res'].set_ylim([-0.1,1.5] if samples=='inner' else [-0.1,1.6])
axes['res'].set_ylabel(r'Resolution ($|3^{\mathrm{rd}}\mathrm{Q}-1^{\mathrm{st}}\mathrm{Q}|$)', fontsize=18*figscale, labelpad=50)
axes['res'].grid(linewidth=3)
axes['res'].legend(loc='upper right' if samples=='outer' else 'upper left')

plt.subplots_adjust(wspace=0., hspace=0.)

axes['res'].grid(axis='x')
axes['bad'].grid(axis='x')
texth = 1.04
plt.text(.0, texth, 'CMS Preliminary', ha='left', transform=axes['res'].transAxes)
radius = dict({'1':'1.3cm', '2':'2.6cm', '3':'5.3cm'})
plt.text(.815, texth, 'Integration cylinder radius: '+radius[region], transform=axes['res'].transAxes)
methodmap = dict({'ed':'Shower reconstruction', 'fineeta':'Brute force calibration', 'nocorr':'No Correction'})
plt.text(.34, texth, 'Software correction method: '+methodmap[method], transform=axes['res'].transAxes)
plt.savefig('/eos/user/b/bfontana/www/PartialWafers/'+samples+'/final2_reg'+str(region)+'_'+samples+'_'+method+extra_str+'.png')
