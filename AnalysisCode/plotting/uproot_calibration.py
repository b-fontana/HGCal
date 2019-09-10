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
file = upopen("root_files/filefriend_"+str(mask)+samples+".root")
#trees = file.allkeys(filterclass=lambda x: issubclass(x, up.tree.TTreeMethods))
tree = file["tfriend1"]
f1 = [tree.array('f1')[:,0], tree.array('f1')[:,1], tree.array('f1')[:,2]]
#f2 = [tree.array('f2')[:,0], tree.array('f2')[:,1], tree.array('f2')[:,2]]

fig, ax = plt.subplots(1, 3, figsize=(15,7), 
                       sharex=False, sharey=False, squeeze=False)
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
for ireg in range(3):
    ax[0,ireg].hist(f1[ireg], bins=100, label='f1')
    #ax[0,ireg].hist(f2[ireg], label='f2')
plt.ylabel('Counts')
plt.legend()
plt.savefig('figs/calibration_'+str(mask)+samples+'.png')
plt.savefig('/eos/user/b/bfontana/www/ResolutionStudies/mask'+str(mask)+'/calibration_'+str(mask)+samples+'.png')
