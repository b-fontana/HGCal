from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib import rcParams as rc
#rc.update({'font.size': 12})
#import seaborn as sns
import uproot as up

file = up.open("../CCode/file_after_weights_3inner.root")
#trees = file.allkeys(filterclass=lambda x: issubclass(x, up.tree.TTreeMethods))
tree = file["data"]
geneta = tree.array('geneta')
showerid = [tree.array('shower_number')[:,0],tree.array('shower_number')[:,1],tree.array('shower_number')[:,2]]

###############################################################################
###############################################################################
width, height = 3, 1
fig, ax = plt.subplots(height, width, figsize=(15,7), 
                       sharex=False, sharey=False, squeeze=False)
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
for iax1 in range(height):
    for iax2 in range(width):        
        h = ax[iax1,iax2].hist(x=showerid[iax2 + width*iax1],
                               bins=(100), range=(-1.2,1.))
plt.xlabel('$\Delta E / E_{gen}$')
plt.savefig('resolution.png')




#g = sns.JointGrid(x=deltaE_corr[:,0], y=geneta)
#g = g.plot_joint(sns.distplot)
#g.savefig('pic.png')
