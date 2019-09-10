import sys
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, OptimizeWarning
import warnings
warnings.simplefilter('error', OptimizeWarning)
import uproot as up
import params

nreg = 3
mask = sys.argv[1]
samples = sys.argv[2]
method = sys.argv[3]
correction = sys.argv[4]
file = up.open('root_files/final_'+str(mask)+samples+'_'+method+'.root')
#trees = file.allkeys(filterclass=lambda x: issubclass(x, up.tree.TTreeMethods))
tree = file["data"]
geneta = tree.array('geneta')
deltaE_array, deltaEcorr_array = tree.array('deltaE'), tree.array('deltaE_corr')

###############################################################################
###############################################################################
width, height = 3, 3
fitcut = 100.5
fig_proj, ax_proj = [[] for _ in range(nreg)], [[] for _ in range(nreg)]
for i in range(nreg):
    fig_proj[i], ax_proj[i] = plt.subplots(4, 5, figsize=(25,15))
fig, ax = plt.subplots(height, width, figsize=(15,12), 
                       squeeze=False)
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
xrangetuple = (2.69,3.04) if samples == 'inner' else (1.44, 1.66)
yrangetuple = (-1.05,1.05) if samples == 'inner' else (-.97, 2.)
gaus = lambda x,a,m,sigma: a*np.exp(-(x-m)**2/(2*sigma**2))
for iax2 in range(width):        
    deltaE_ = deltaEcorr_array[:,iax2] if correction else deltaE_array[:,iax2]
    xbins = params.etavalues(samples)
    ybins = params.resvalues()
    options = dict(x=geneta, bins=(xbins,ybins), range=(xrangetuple, yrangetuple),
                   cmap='YlOrRd', norm=LogNorm())
    h, xedges, yedges, _ = ax[0,iax2].hist2d(y=deltaE_, **options)
    bias = np.zeros((len(xbins)-1))
    biaserr = np.zeros((len(xbins)-1))
    rms = np.zeros((len(xbins)-1))
    median = np.zeros((len(xbins)-1))
    stddev = np.zeros((len(xbins)-1))
    nsamples = np.zeros((len(xbins)-1))

    xy = np.array([np.array(deltaE_), np.array(geneta)])
    it = 0
    
    for ix1,ix2 in zip(xedges[:-1],np.roll(xedges,-1)[:-1]):
        tmp = xy[0][(xy[1]>=ix1) & (xy[1]<ix2)]
        rms[it] = np.sqrt(sum(tmp**2)/len(tmp))
        median[it] = np.median(tmp)
        stddev[it] = np.std(tmp)
        nsamples[it] = len(tmp)
        it += 1
    xcenters = (xedges[:-1]+np.roll(xedges,-1)[:-1])/2
    ycenters = (yedges[:-1]+np.roll(yedges,-1)[:-1])/2
    ycenters_cut = ycenters[abs(ycenters)<fitcut]
    for ibin in range(len(xbins)-1):
        h_cut = np.array(h[ibin,:])[abs(ycenters)<fitcut]
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
                biaserr[ibin] = 1.253*stddev[ibin]/np.sqrt(nsamples[ibin])
            else:
                bias[ibin] = popt[1]
                biaserr[ibin] = np.sqrt(np.diag(pcov)[1])
                ax_proj[iax2][int(np.floor(ibin/5)),int(ibin%5)].plot(ycenters_cut, gaus(ycenters_cut,*popt), color='red') 
        except (RuntimeError, OptimizeWarning):
            try:
                popt, pcov = curve_fit(gaus, ycenters_cut, h_cut, p0=[150,mean,sigma/4.], maxfev=10000,
                                       bounds=([0.,-np.inf,0.], np.inf))
                if abs(popt[1]) > .4 or popt[2]>.7 or nsamples[ibin]<10 or len(ycenters[ycenters<-0.95]) > 0.9*len(ycenters):
                    bias[ibin] = median[ibin]
                    biaserr[ibin] = 1.253*stddev[ibin]/np.sqrt(nsamples[ibin])
                else:
                    bias[ibin] = popt[1]
                    biaserr[ibin] = np.sqrt(np.diag(pcov)[1])
                    ax_proj[iax2][int(np.floor(ibin/5)),int(ibin%5)].plot(ycenters_cut, gaus(ycenters_cut,*popt), color='orange') 
            except (RuntimeError, OptimizeWarning):
                bias[ibin] = median[ibin]
                biaserr[ibin] = 1.253*stddev[ibin]/np.sqrt(nsamples[ibin])
    indep = [rms[i]/(1+bias[i]) if bias[i]!=-1 else -99 for i in range(len(bias))]
    indep_error = biaserr
    correction_str = '_corr' if correction else '_nocorr'
    save_str = ( 'numpy_files/arrays_reg' + str(iax2+1) + '_'
                 + str(mask) + str(samples) + correction_str + '_' + method )
    np.savez(save_str, 
             xcenters=xcenters, bias=bias, 
             indep=indep, indep_error=indep_error, 
             biaserr=biaserr, response_eta=xy)
    ax[1,iax2].errorbar(xcenters, bias, yerr=biaserr, color='green', ecolor='green', 
                        label='bias', linestyle='', markersize=3, capsize=3)
    ax[1,iax2].legend()
    ax[2,iax2].set_ylim([-0.5,0.5])
    ax[2,iax2].errorbar(xcenters, indep, yerr=indep_error, color='blue', ecolor='blue', 
                        label='rms/(1+bias)', linestyle='', markersize=3, capsize=3)
    ax[2,iax2].legend()

for i in range(nreg):
    fig_proj[i].savefig('/eos/user/b/bfontana/www/ResolutionStudies/proj_'+method+'_reg'+str(i+1)+'_'+str(mask)+samples+'.png')
    fig_proj[i].savefig('figs/proj_'+method+'_reg'+str(i+1)+'_'+str(mask)+samples+'.png')

plt.xlabel('$|\eta|$')
plt.gca().xaxis.set_label_coords(0.4,-0.05)
plt.ylabel(r'$(E_{reco}-E_{gen})/E_{gen}$')
plt.savefig('figs/resolution_vs_eta_'+str(mask)+samples+'_'+method+'.png')
plt.savefig('/eos/user/b/bfontana/www/ResolutionStudies/mask'+str(mask)+'/resolution_vs_eta_'+str(mask)+samples+'_'+method+'.png')

###############################################################################
###############################################################################
width, height = 3, 1 if method == 'fineeta' else 2
fig, ax = plt.subplots(height, width, figsize=(15,7), 
                       sharex=False, sharey=False, squeeze=False)
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
xrangetuple = (-1.05,1.05) if samples == 'inner' else (-.97,2.)
label = 'total resolution after fine calibration'
options = dict(bins=(100), range=xrangetuple, histtype='step', color='blue')
for iax1 in range(height):
    for iax2 in range(width):
        if method == 'ed':
            ax[iax1,iax2].hist(x=deltaE_array[:,iax2], 
                               label='uncorrected', **options)
        ax[iax1,iax2].hist(x=deltaEcorr_array[:,iax2], 
                           label='corrected', **options)
        ax[iax1,iax2].legend()
plt.xlabel(r'$(E_{reco}-E_{gen})/E_{gen}$')
plt.ylabel('Counts')
plt.savefig('figs/resolution_'+str(mask)+samples+'_'+method+'.png')
plt.savefig('/eos/user/b/bfontana/www/ResolutionStudies/mask'+str(mask)+'/resolution_'+str(mask)+samples+'_'+method+'.png')
