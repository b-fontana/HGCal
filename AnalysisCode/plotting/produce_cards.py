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

method, region = sys.argv[1], sys.argv[2]
df = []
bias, biaserr, bad75, bad50, nsamples, quart_dist = ([] for _ in range(6))

#Inner region
etavalues = params.etavalues('inner')
method_str = 'corr_fineeta' if method == 'fineeta' else 'corr_ed' if method == 'ed' else 'nocorr'
for imask in range(3,7):
    with np.load('numpy_files/arrays_reg'+str(region)+'_'+str(imask)+'inner_'+method_str+'.npz') as f:
        df.append(pd.DataFrame(f['response_eta']).transpose())
        df[-1].columns = ['response', 'eta']
        df[-1] = df[-1][ (df[-1].response>=-1) & (df[-1].response<2.5) ]
        df[-1]['eta_cat'] = pd.cut(df[-1].eta, bins=etavalues, right=False, labels=[i for i in range(len(etavalues)-1)])
        df[-1]['mask'] = imask
        bias.append(f['bias']) 
        quantiles = df[-1].groupby('eta_cat')['response'].quantile([.25,.75]).to_numpy()
        quart_dist.append(np.abs(quantiles[0::2]-quantiles[1::2]))

card_etavals = etavalues

#Outer region
etavalues = params.etavalues('outer')
method_str = 'corr_fineeta' if method == 'fineeta' else 'corr_ed' if method == 'ed' else 'nocorr'
for imask in range(3,7):
    with np.load('numpy_files/arrays_reg'+str(region)+'_'+str(imask)+'outer_'+method_str+'.npz') as f:
        df.append(pd.DataFrame(f['response_eta']).transpose())
        df[-1].columns = ['response', 'eta']
        df[-1] = df[-1][ (df[-1].response>=-1) & (df[-1].response<2.5) ]
        df[-1]['eta_cat'] = pd.cut(df[-1].eta, bins=etavalues, right=False, labels=[i for i in range(len(etavalues)-1)])
        df[-1]['mask'] = imask
        
        #the NaN's are required due to the interface between the outer and inner region (the bin between 1.65 and 2.7 should not exist)
        nan_array = np.empty((1))
        nan_array.fill(np.nan)
        bias[imask-3] = np.concatenate( (f['bias'], nan_array, bias[imask-3]), axis=0 )
        quantiles = df[-1].groupby('eta_cat')['response'].quantile([.25,.75]).to_numpy()
        quart_dist[imask-3] = np.concatenate( (np.abs(quantiles[0::2]-quantiles[1::2]), nan_array, quart_dist[imask-3]), axis=0 )

card_etavals = np.concatenate((etavalues, card_etavals), axis=0)
dump_path = '/afs/cern.ch/user/b/bfontana/CMSSW_10_6_0/src/UserCode/DelphesNtuplizer/delphes/cards/'
if not os.path.isdir(dump_path):
    raise ValueError('Wrong path.')
#dump_tcl_resolution_formula(dump_path+'photonEnergyResolution', card_etavals, quart_dist, bias)
nan_array = np.empty((4,1))
nan_array.fill(np.nan)
zeros = np.concatenate( (np.zeros((4,14)), nan_array, np.zeros((4,13))), axis=1 )
dump_tcl_resolution_formula(dump_path+'photonEnergyResolution', card_etavals, quart_dist, zeros)

"""
def save_to_csv(name, bins, res, biases):
    assert(len(res)==4 and len(biases)==4)
    for r in res:
        assert(len(bins)==len(r)+1)
    for bias in biases:
        assert(len(bins)==len(bias)+1)
    with open(name, 'w') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',')
        csv_writer.writerow(['Edge1', 'Edge2', 
                             'ResolutionMask3', 'ResolutionMask4', 'ResolutionMask5', 'ResolutionMask6',
                             'BiasMask3', 'BiasMask4', 'BiasMask5', 'BiasMask6'])
        for i, (edge1,edge2) in enumerate(zip(bins[:-1],bins[1:])):
            csv_writer.writerow([edge1, edge2, 
          round(res[0][i]/res[0][i],3), round(res[1][i]/res[0][i],3), round(res[2][i]/res[0][i],3),  round(res[3][i]/res[0][i],3),
          round(biases[0][i],3), round(biases[1][i],3), round(biases[2][i],3),  round(biases[3][i],3)])
"""
