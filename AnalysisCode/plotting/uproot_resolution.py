import sys
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib import rcParams as rc
#rc.update({'font.size': 12})
#import seaborn as sns
import uproot as up

mask = sys.argv[1]
samples = sys.argv[2]
file = up.open("root_files/file_after_weights_"+str(mask)+samples+".root")
#trees = file.allkeys(filterclass=lambda x: issubclass(x, up.tree.TTreeMethods))
tree = file["data"]
geneta = tree.array('geneta')
deltaE = [tree.array('deltaE')[:,0], tree.array('deltaE')[:,1], 
          tree.array('deltaE')[:,2], tree.array('deltaE_corr')[:,0], 
          tree.array('deltaE_corr')[:,1], tree.array('deltaE_corr')[:,2]]

###############################################################################
###############################################################################
width, height = 3, 2
fig, ax = plt.subplots(height, width, figsize=(15,7), 
                       sharex=False, sharey=False, squeeze=False)
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
for iax1 in range(height):
    for iax2 in range(width):        
        h = ax[iax1,iax2].hist2d(x=geneta, y=deltaE[iax2 + width*iax1],
                                 bins=(10,100), range=( (2.69,3.04), (-1.2,1.)),
                                 cmap='YlOrRd', norm=LogNorm())
plt.colorbar(h[3], ax=ax.ravel().tolist())
plt.xlabel('$|\eta|$')
plt.gca().xaxis.set_label_coords(0.4,-0.05)
plt.ylabel('$\Delta E / E_{gen}$')
plt.savefig('resolution_vs_eta.png')

###############################################################################
###############################################################################
width, height = 3, 2
fig, ax = plt.subplots(height, width, figsize=(15,7), 
                       sharex=False, sharey=False, squeeze=False)
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
for iax1 in range(height):
    for iax2 in range(width):        
        h = ax[iax1,iax2].hist(x=deltaE[iax2 + width*iax1],
                               bins=(100), range=(-1.2,1.))
plt.xlabel('$\Delta E / E_{gen}$')
plt.savefig('figs/resolution_'+str(mask)+samples+'.png')
plt.savefig('/eos/user/b/bfontana/www/ResolutionStudies/resolution_'+str(mask)+samples+'.png')
