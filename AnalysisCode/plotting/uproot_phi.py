import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import uproot as up

def reduce_radian_range(v, r):
    """converts all values in 'v' to the -r to r range"""
    assert r>0
    f = np.floor(abs(v)/r)
    prod = f*r
    return [v[i]-prod[i] if v[i]>0 else v[i]+prod[i] for i in range(len(v))]

nreg = 3
mask = sys.argv[1]
samples = sys.argv[2]
bound = np.pi/6
file = up.open("root_files/final_"+str(mask)+samples+"_ed.root")

tree = file["data"]
geneta = tree.array('geneta')
genphi = reduce_radian_range(tree.array('genphi'), bound)
assert len(geneta) == len(genphi)
deltaE_array = tree.array('deltaE')
deltaEcorr_array = tree.array('deltaE_corr')
deltaE = [deltaE_array[:,0], deltaE_array[:,1], deltaE_array[:,2]]
deltaE_corr = [deltaEcorr_array[:,0], deltaEcorr_array[:,1], deltaEcorr_array[:,2]]

if samples == 'inner':
    etacut = 2.9
    genphi = [genphi[i] for i in range(len(geneta)) if geneta[i] > etacut]
    for ie in range(len(deltaE)):
        deltaE[ie] = [deltaE[ie][i] for i in range(len(geneta)) if geneta[i] > etacut]
    for ie in range(len(deltaE)):
        deltaE_corr[ie] = [deltaE_corr[ie][i] for i in range(len(geneta)) if geneta[i] > etacut]
elif samples == 'outer':
    etacut = 1.55
    genphi = [genphi[i] for i in range(len(geneta)) if geneta[i] < etacut]
    for ie in range(len(deltaE)):
        deltaE[ie] = [deltaE[ie][i] for i in range(len(geneta)) if geneta[i] < etacut]
    for ie in range(len(deltaE)):
        deltaE_corr[ie] = [deltaE_corr[ie][i] for i in range(len(geneta)) if geneta[i] < etacut]
        
###############################################################################
###############################################################################
width, height = 3, 1
fig, ax = plt.subplots(height, width, figsize=(15,7), 
                       sharex=False, sharey=False, squeeze=False)
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
bins = (40,25)
xrangetuple = (-bound,bound)
yrangetuple = (-.2,.2) if samples=='inner' else (-.4,.2)
for iax1 in range(height):
    for iax2 in range(width):
        h = ax[0,iax2].hist2d(x=genphi, y=deltaE_corr[iax2], bins=bins,
                              range=(xrangetuple,yrangetuple),
                              cmap='YlOrRd', norm=LogNorm())
        ax[iax1,iax2].legend()
fig.colorbar(h[3], ax=ax.ravel().tolist())
plt.xlabel(r'$\phi$')
plt.ylabel(r'$(E_{reco}-E_{gen})/E_{gen}$')
plt.savefig('figs/phi_'+str(mask)+samples+'_ed.png')
plt.savefig('/eos/user/b/bfontana/www/ResolutionStudies/mask'+str(mask)+'/phi_'+str(mask)+samples+'_ed.png')
