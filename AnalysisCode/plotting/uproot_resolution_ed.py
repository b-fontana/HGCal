import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
#from matplotlib import rcParams as rc
#rc.update({'font.size': 12})
#import seaborn as sns
from scipy.optimize import curve_fit, OptimizeWarning
import warnings
warnings.simplefilter('error', OptimizeWarning)
from scipy import exp
import uproot as up

nreg = 3
mask = sys.argv[1]
samples = sys.argv[2]
file = up.open("root_files/final_"+str(mask)+samples+"_ed.root")
#trees = file.allkeys(filterclass=lambda x: issubclass(x, up.tree.TTreeMethods))
tree = file["data"]
geneta = tree.array('geneta')
deltaE_array = tree.array('deltaE')
deltaEcorr_array = tree.array('deltaE_corr')
deltaE = [deltaE_array[:,0], deltaE_array[:,1], deltaE_array[:,2]]
deltaE_corr = [deltaEcorr_array[:,0], deltaEcorr_array[:,1], deltaEcorr_array[:,2]]

###############################################################################
###############################################################################
width, height = 3, 3
fitcut = 1.5
fig_proj, ax_proj = [[] for _ in range(nreg)], [[] for _ in range(nreg)]
for i in range(nreg):
    fig_proj[i], ax_proj[i] = plt.subplots(4, 5, figsize=(25,15))
fig, ax = plt.subplots(height, width, figsize=(15,12), 
                       sharex=False, sharey=False, squeeze=False)
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
xrangetuple = (2.69,3.04) if samples == 'inner' else (1.44, 1.66)
yrangetuple = (-1.05,1.05) if samples == 'inner' else (-.99, 3.)
gaus = lambda x,a,m,sigma: a*exp(-(x-m)**2/(2*sigma**2))
for iax2 in range(width):        
    bins = (12,100) if samples == 'inner' else (12,200)
    h, xedges, yedges, _ = ax[0,iax2].hist2d(x=geneta, y=deltaE_corr[iax2], bins=bins,
                                             range=(xrangetuple, yrangetuple),
                                             cmap='YlOrRd', norm=LogNorm())
    cols_with_zero = 0
    for ibin in range(bins[0]):
        if not np.any(h[ibin,:]):
            cols_with_zero += 1
    bias = np.zeros((bins[0])      - cols_with_zero)
    sigmas = np.zeros((bins[0])    - cols_with_zero)
    rms = np.zeros((bins[0])       - cols_with_zero)
    median = np.zeros((bins[0])    - cols_with_zero)
    stddev = np.zeros((bins[0])    - cols_with_zero)
    nsamples = np.zeros((bins[0])  - cols_with_zero)

    xy = np.array([np.array(deltaE_corr[iax2]), np.array(geneta)])
    it = 0
    for ix1,ix2 in zip(xedges[:-1],np.roll(xedges,-1)[:-1]):
        if ix1 < xedges[cols_with_zero]: 
            continue
        tmp = xy[0][(xy[1]>=ix1) & (xy[1]<ix2)]
        rms[it] = np.sqrt(sum(tmp**2)/len(tmp))
        median[it] = np.median(tmp)
        stddev[it] = np.std(tmp)
        nsamples[it] = len(tmp)
        it += 1
    xcenters = (xedges[:-1]+np.roll(xedges,-1)[:-1])/2
    ycenters = (yedges[:-1]+np.roll(yedges,-1)[:-1])/2
    ycenters_cut = ycenters[abs(ycenters)<fitcut]
    for ibin in range(bins[0]-cols_with_zero):
        h_cut = np.array(h[ibin+cols_with_zero,:])[abs(ycenters)<fitcut]
        mean = h_cut.dot(ycenters_cut)/sum(h_cut)
        sigma = h_cut.dot(ycenters_cut**2)/sum(h_cut) - mean**2
        ax_proj[iax2][int(np.floor(ibin/5)),int(ibin%5)].plot(ycenters_cut, h_cut, label=str(ibin)+', reg'+str(iax2+1))
        try:
            popt, pcov = curve_fit(gaus, ycenters_cut, h_cut, p0=[150,mean,sigma/2.],
                                   bounds=([0.,-np.inf,0.], np.inf))
            if abs(popt[1]) > .8 or popt[2]>2.:
                raise RuntimeError()
            elif abs(popt[1]) > .4 or popt[2]>.7 or nsamples[ibin]<10 or len(ycenters[ycenters<-0.95]) > 0.9*len(ycenters):
                #the last conditions select situations where more than 90% of the data points lie close to resolution=-1
                bias[ibin] = median[ibin]
                sigmas[ibin] = 1.253*stddev[ibin]/np.sqrt(nsamples[ibin])
            else:
                bias[ibin] = popt[1]
                sigmas[ibin] = popt[2]
                ax_proj[iax2][int(np.floor(ibin/5)),int(ibin%5)].plot(ycenters_cut, gaus(ycenters_cut,*popt), color='red') 
        except (RuntimeError, OptimizeWarning):
            try:
                popt, pcov = curve_fit(gaus, ycenters_cut, h_cut, p0=[150,mean,sigma/4.], maxfev=10000,
                                       bounds=([0.,-np.inf,0.], np.inf))
                if abs(popt[1]) > .4 or popt[2]>.7 or nsamples[ibin]<10 or len(ycenters[ycenters<-0.95]) > 0.9*len(ycenters):
                    bias[ibin] = median[ibin]
                    sigmas[ibin] = 1.253*stddev[ibin]/np.sqrt(nsamples[ibin])
                else:
                    bias[ibin] = popt[1]
                    sigmas[ibin] = popt[2]
                    ax_proj[iax2][int(np.floor(ibin/5)),int(ibin%5)].plot(ycenters_cut, gaus(ycenters_cut,*popt), color='orange') 
            except (RuntimeError, OptimizeWarning):
                bias[ibin] = median[ibin]
                sigmas[ibin] = 1.253*stddev[ibin]/np.sqrt(nsamples[ibin])
    indep = [rms[i]/(1+bias[i]) if bias[i]!=-1 else -99 for i in range(len(bias))]
    indep_error = sigmas
    np.savez('numpy_files/arrays_reg'+str(iax2+1)+'_'+str(mask)+str(samples)+'_ed', 
             xcenters=xcenters[cols_with_zero:], bias=bias, sigmas=sigmas,
             indep=indep, indep_error=indep_error)
    ax[1,iax2].errorbar(xcenters[cols_with_zero:], bias, yerr=sigmas, 
                        color='green', ecolor='green', 
                        label='bias', linestyle='', markersize=3, capsize=3)
    ax[1,iax2].legend()
    #ax[2,iax2].set_ylim([-0.5,0.5])
    ax[2,iax2].errorbar(xcenters[cols_with_zero:], indep, yerr=indep_error, 
                        color='blue', ecolor='blue', 
                        label='rms/(1+bias)', linestyle='', markersize=3, capsize=3)
    ax[2,iax2].set_ylim([-10.,10.])
    ax[2,iax2].legend()

fig_proj[0].savefig('/eos/user/b/bfontana/www/ResolutionStudies/proj_ed_reg1_'+str(mask)+samples+'.png')
fig_proj[1].savefig('/eos/user/b/bfontana/www/ResolutionStudies/proj_ed_reg2_'+str(mask)+samples+'.png')
fig_proj[2].savefig('/eos/user/b/bfontana/www/ResolutionStudies/proj_ed_reg3_'+str(mask)+samples+'.png')
fig_proj[0].savefig('figs/proj_ed_reg1_'+str(mask)+samples+'.png')
fig_proj[1].savefig('figs/proj_ed_reg2_'+str(mask)+samples+'.png')
fig_proj[2].savefig('figs/proj_ed_reg3_'+str(mask)+samples+'.png')
#plt.colorbar(h[3], ax=ax.ravel().tolist())
plt.xlabel('$|\eta|$')
plt.gca().xaxis.set_label_coords(0.4,-0.05)
plt.ylabel(r'$(E_{reco}-E_{gen})/E_{gen}$')
plt.savefig('figs/resolution_vs_eta_'+str(mask)+samples+'_ed.png')
plt.savefig('/eos/user/b/bfontana/www/ResolutionStudies/mask'+str(mask)+'/resolution_vs_eta_'+str(mask)+samples+'_ed.png')

###############################################################################
###############################################################################
width, height = 3, 2
fig, ax = plt.subplots(height, width, figsize=(15,7), 
                       sharex=False, sharey=False, squeeze=False)
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
xrangetuple = (-1.05,1.05) if samples == 'inner' else (-.99,1.05)

for iax1 in range(height):
    for iax2 in range(width):
        if iax1 == 0:
            ax[iax1,iax2].hist(x=deltaE[iax2], bins=(100), range=xrangetuple, 
                               histtype='step', color='blue', 
                               label='total resolution before ed correction')
        elif iax1 == 1:
            ax[iax1,iax2].hist(x=deltaE_corr[iax2], bins=(100), range=xrangetuple, 
                               histtype='step', color='blue',
                               label='total resolution after ed correction')
        ax[iax1,iax2].legend()
plt.xlabel(r'$(E_{reco}-E_{gen})/E_{gen}$')
plt.ylabel('Counts')
plt.savefig('figs/resolution_'+str(mask)+samples+'_ed.png')
plt.savefig('/eos/user/b/bfontana/www/ResolutionStudies/mask'+str(mask)+'/resolution_'+str(mask)+samples+'_ed.png')
