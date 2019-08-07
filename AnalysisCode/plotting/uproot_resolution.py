import sys
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
#from matplotlib import rcParams as rc
#rc.update({'font.size': 12})
#import seaborn as sns
import uproot as up

mask = sys.argv[1]
samples = sys.argv[2]
file = up.open("root_files/file_after_weights_"+str(mask)+samples+".root")
#trees = file.allkeys(filterclass=lambda x: issubclass(x, up.tree.TTreeMethods))
tree = file["data"]
geneta = tree.array('geneta')
deltaE_array = tree.array('deltaE')
deltaEcorr_array = tree.array('deltaE_corr')
deltaE = [deltaE_array[:,0], deltaE_array[:,1], deltaE_array[:,2], 
          deltaEcorr_array[:,0], deltaEcorr_array[:,1], deltaEcorr_array[:,2]]
showerid_array = tree.array('shower_number')
showerid = [showerid_array[:,0], showerid_array[:,1], showerid_array[:,2]]
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
plt.ylabel(r'$(E_{reco}-E_{gen})/E_{gen}$')
plt.savefig('resolution_vs_eta_'+str(mask)+samples+'.png')
plt.savefig('/eos/user/b/bfontana/www/ResolutionStudies/resolution_vs_eta_'+str(mask)+samples+'.png')

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
plt.xlabel(r'$(E_{reco}-E_{gen})/E_{gen}$')
plt.ylabel('Counts')
plt.savefig('figs/resolution_'+str(mask)+samples+'.png')
plt.savefig('/eos/user/b/bfontana/www/ResolutionStudies/resolution_'+str(mask)+samples+'.png')

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
        if iax1==0:
            x1 = [elem for i,elem in enumerate(deltaE[iax2]) if showerid[iax2][i]==0]
            ax[iax1,iax2].hist(x=x1, bins=(200), range=(-1.1,1.1), label='complete before correction')
            x2 = [elem for i,elem in enumerate(deltaE[iax2 + 3]) if showerid[iax2][i]==0]
            ax[iax1,iax2].hist(x=x2, bins=(200), range=(-1.1,1.1), label='complete after correction')
        elif iax1==1:
            x1 = [elem for i,elem in enumerate(deltaE[iax2]) if showerid[iax2][i]>0]
            ax[iax1,iax2].hist(x=x1, bins=(200), range=(-1.1,1.1), label='incomplete before correction')
            x2 = [elem for i,elem in enumerate(deltaE[iax2 + 3]) if showerid[iax2][i]>0]
            ax[iax1,iax2].hist(x=x2, bins=(200), range=(-1.1,1.1), label='incomplete after correction')
        ax[iax1,iax2].legend()
plt.ylabel('Counts')
plt.xlabel(r'$(E_{reco}-E_{gen})/E_{gen}$')
plt.savefig('figs/res_split_'+str(mask)+samples+'.png')
plt.savefig('/eos/user/b/bfontana/www/ResolutionStudies/res_split_'+str(mask)+samples+'.png')
