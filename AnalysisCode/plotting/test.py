import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
sns.set(style='ticks', font_scale=.9, rc={'font.sans-serif': 'Bitstream Vera Sans'})
sns.despine()

method, samples = sys.argv[1], sys.argv[2]
df = []
bias = []
etavalues = np.linspace(1.45,1.55,12)
#etavalues = np.linspace(1.45,1.55,4)
etalabels = ['['+str(np.round(x,3))+';'+str(np.round(y,3))+'[' for x,y in zip(etavalues[:-1],np.roll(etavalues,-1)[:-1])]
for imask in range(3,7):
    with np.load('numpy_files/arrays_reg3_'+str(imask)+samples+'_'+method+'.npz') as f:
        df.append(pd.DataFrame(f['response_eta']).transpose())
        df[-1].columns = ['response', 'eta']
        df[-1] = df[-1][ (df[-1].response>=-1) & (df[-1].response<3) ]
        df[-1]['eta_cat'] = pd.cut(df[-1].eta, bins=len(etavalues)-1, right=False, labels=[i for i in range(len(etavalues)-1)])
        df[-1]['mask'] = imask
        bias.append(f['bias'])

#join the datasets of the fours geometries (masks)
df_tot = pd.concat([df[0],df[1],df[2],df[3]], axis=0, sort=False)

#plot
plt.figure(figsize=(25,12))
plt.plot([-0.5, 40], [0, 0], linewidth=2, color='grey', linestyle='dashed')
ax = sns.boxplot(x='eta_cat', y='response', hue='mask', data=df_tot, whis=1.5)
ax.set_xticklabels(etalabels)
plt.plot(np.linspace(-0.3,10.7,12), bias[0], color='blue', linestyle='', marker='X', markersize=12)
plt.plot(np.linspace(-0.1,10.9,12), bias[1], color='orange', linestyle='', marker='X', markersize=12)
plt.plot(np.linspace(0.1,11.1,12), bias[2], color='green', linestyle='', marker='X', markersize=12)
plt.plot(np.linspace(0.3,11.3,12), bias[3], color='red', linestyle='', marker='X', markersize=12)
ax.figure.savefig('pic.png')
ax.figure.savefig('/eos/user/b/bfontana/www/ResolutionStudies/pic.png')
