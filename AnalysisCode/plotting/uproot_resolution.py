import sys
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, OptimizeWarning
import warnings
warnings.simplefilter('error', OptimizeWarning)
import uproot as up
import params

extra_str = '_new'
nreg = 3
mask = sys.argv[1]
samples = sys.argv[2]
method = sys.argv[3]
correction = False if int(sys.argv[4])==0 else True
correction_str = '_corr_'+method if correction else '_nocorr'
file = up.open('root_files/final_'+str(mask)+samples+'_'+method+'.root')
#trees = file.allkeys(filterclass=lambda x: issubclass(x, up.tree.TTreeMethods))
tree = file["data"]
geneta = tree.array('abs_geneta')
deltaE_array, deltaEcorr_array = tree.array('deltaE'), tree.array('deltaE_corr')

###############################################################################
###############################################################################
width, height = 3, 3
fitcut = 2.
fig_proj, ax_proj = [[] for _ in range(nreg)], [[] for _ in range(nreg)]
for i in range(nreg):
    fig_proj[i], ax_proj[i] = plt.subplots(3, 5, figsize=(25,15),
                                           sharex=True, sharey=False, squeeze=False)
fig, ax = plt.subplots(1, width+1, figsize=(13,7), squeeze=False)
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
    h, xedges, yedges, cbar = ax[0,iax2].hist2d(y=deltaE_, **options)
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
        if(len(tmp)==0):
            rms[it] = 0.
            median[it] = -1.
        else:
            rms[it] = np.sqrt(sum(tmp**2)/len(tmp))
            if(len(tmp[tmp>-.9])==0):
                median[it] = np.median(tmp)
            else:
                median[it] = np.median(tmp[tmp>-.9])
        stddev[it] = np.std(tmp)
        nsamples[it] = len(tmp)
        it += 1

    xcenters = (xedges[:-1]+np.roll(xedges,-1)[:-1])/2
    ycenters = (yedges[:-1]+np.roll(yedges,-1)[:-1])/2
    ycenters_cut = ycenters[abs(ycenters)<fitcut]
    for ibin in range(len(xbins)-1):
        h_cut = np.array(h[ibin,:])[abs(ycenters)<fitcut]
        if(h_cut.size==0):
            raise ValueError('This bin is empty!')
        print('value', sum(h_cut))
        mean = h_cut.dot(ycenters_cut)/sum(h_cut)
        sigma = h_cut.dot(ycenters_cut**2)/sum(h_cut) - mean**2
        ax_proj[iax2][int(np.floor(ibin/5)),int(ibin%5)].plot(ycenters_cut, h_cut, label=str(ibin)+', reg'+str(iax2+1))

        if nsamples[ibin]<10 or len(ycenters[ycenters<-0.9]) > 0.9*len(ycenters):
            if(nsamples[ibin]<2):
                raise ValueError('Not enough samples!')
            #the last conditions select situations where more than 90% of the data points lie close to resolution=-1
            bias[ibin] = median[ibin]
            biaserr[ibin] = 1.253*stddev[ibin]/np.sqrt(nsamples[ibin])
        else:
            if len(ycenters[ycenters<-0.8]) > 0.4*len(ycenters):
                h_cut = h_cut[h_cut>-0.8]
            else:
                h_cut = h_cut[h_cut>-0.9]
            remove_values = [-.75, -.7, -.65, -.6]
            ntries, ntrieslimit = 0, 20
            nremovals = 0
            qualdiff = .1
            previous_quality = np.inf
            bad_fit = True
            while bad_fit:
                try:
                    popt, pcov = curve_fit(gaus, ycenters_cut, h_cut, 
                                           p0=[150,mean,sigma],
                                           bounds=([0.,-abs(1.5*mean)-0.01,0.], [np.inf,abs(1.5*mean)+0.01,2*sigma+0.005]), 
                                           maxfev=5000)
                except (RuntimeError, ValueError):
                    bias[ibin] = median[ibin]
                    biaserr[ibin] = 1.253*stddev[ibin]/np.sqrt(nsamples[ibin])
                    bad_fit = False
                    print("end: RuntimeError")
                    continue
                ntries += 1
                if ntries > ntrieslimit:
                    bias[ibin] = median[ibin]
                    biaserr[ibin] = 1.253*stddev[ibin]/np.sqrt(nsamples[ibin])
                    bad_fit = False
                    print("end: maximum number of tries exceeded")
                    continue

                if abs(popt[1]) > .4 or popt[2]>1.:
                    bad_fit = True
                    print("abs(popt[1]) > .4 or popt[2]>1.")
                    continue
                else:
                    curr_quality = sum((gaus(ycenters_cut, *popt) - h_cut)**2)
                    print(curr_quality)
                    if previous_quality - curr_quality > qualdiff:
                        previous_quality = curr_quality
                        mean = popt[1]
                        sigma = popt[2]
                        bad_fit = True
                        print("1: previous_quality - curr_quality > qualdiff")
                        continue
                    else:
                        nremovals += 1
                        h_cut = h_cut[h_cut>remove_values[nremovals]]
                        popt, pcov = curve_fit(gaus, ycenters_cut, h_cut, 
                                               p0=[150,mean,sigma], 
                                               bounds=([0.,-abs(2*mean)-0.01,0.], [np.inf, abs(2*mean)+0.01, 2*sigma+0.005]), 
                                               maxfev=5000)
                        curr_quality = sum((gaus(ycenters_cut, *popt) - h_cut)**2)
                        if previous_quality - curr_quality > qualdiff and nremovals<len(remove_values):
                            previous_quality = curr_quality
                            mean = popt[1]
                            sigma = popt[2]
                            bad_fit = True
                            print("2: previous_quality - curr_quality > qualdiff")
                            continue

                    bad_fit = False
                    print("end")
                    bias[ibin] = popt[1]
                    biaserr[ibin] = np.sqrt(np.diag(pcov)[1])
                    curr_axis = ax_proj[iax2][int(np.floor(ibin/5)),int(ibin%5)]
                    curr_axis.text(.01, 1.02, '#Iterations / Maximum: '+str(ntries)+' / '+str(ntrieslimit), transform=curr_axis.transAxes)
                    curr_axis.plot(ycenters_cut, gaus(ycenters_cut,*popt), color='red')

    indep = [rms[i]/(1+bias[i]) if bias[i]!=-1 else -99 for i in range(len(bias))]
    indep_error = biaserr
    save_str = ( 'numpy_files/arrays_reg' + str(iax2+1) + '_'
                 + str(mask) + str(samples) + correction_str )
    np.savez(save_str, xcenters=xcenters, bias=bias, indep=indep, indep_error=indep_error,
             biaserr=biaserr, response_eta=xy)
    """
    ax[1,iax2].errorbar(xcenters, bias, yerr=biaserr, color='green', ecolor='green',
                        label='bias', linestyle='', markersize=3, capsize=3)
    ax[1,iax2].legend()
    ax[2,iax2].errorbar(xcenters, indep, yerr=indep_error, color='blue', ecolor='blue',
                        label='rms/(1+bias)', linestyle='', markersize=3, capsize=3)
    ax[2,iax2].legend()
    """
fig.colorbar(cbar, ax=ax.ravel().tolist())
fig.subplots_adjust(wspace=0.2, hspace=0.1)
fig.savefig('/eos/user/b/bfontana/www/PartialWafers/'+samples+'/mask'+str(mask)+'/test'+correction_str+'_reg'+str(i+1)+'_'+str(mask)+samples+extra_str+'.png')
for i in range(nreg):
    fig_proj[i].subplots_adjust(wspace=0.2, hspace=0.1)
    fig_proj[i].savefig('/eos/user/b/bfontana/www/PartialWafers/'+samples+'/mask'+str(mask)+'/proj'+correction_str+'_reg'+str(i+1)+'_'+str(mask)+samples+extra_str+'.png')

"""
plt.xlabel('$|\eta|$')
plt.gca().xaxis.set_label_coords(0.4,-0.05)
plt.ylabel(r'$(E_{reco}-E_{gen})/E_{gen}$')
plt.savefig('/eos/user/b/bfontana/www/PartialWafers/'+samples+'/mask'+str(mask)+'/resolution_vs_eta_'+str(mask)+samples+'_'+method+'.png')
"""
###############################################################################
###############################################################################
width, height = nreg, 2
fig, ax = plt.subplots(height, width, figsize=(15,7),
                       sharex=True, sharey=False, squeeze=False)
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
xrangetuple = (-1.05,2.50) if samples == 'inner' else (-.97,2.)
label = 'total resolution after fine calibration'
options = dict(bins=(100), range=xrangetuple, histtype='step')
for iax1 in range(width):
    frac1 = float(deltaEcorr_array[:,iax1][deltaEcorr_array[:,iax1]<-0.5].size) / deltaE_array[:,iax1].size
    plt.text(.02, .89, r'$\frac{N_{\mathrm{R}<-0.5}}{N_{\mathrm{Total}}}=$'+str(round(frac1*100,2))+'%',
             transform=ax[0,iax1].transAxes, fontsize=12)
    ax[0,iax1].hist(x=deltaEcorr_array[:,iax1], label='corrected', color='olive', **options)
    ax[0,iax1].legend()

    frac2 = float(deltaE_array[:,iax1][deltaE_array[:,iax1]<-0.5].size) / deltaE_array[:,iax1].size
    plt.text(.02, .89, r'$\frac{N_{\mathrm{R}<-0.5}}{N_{\mathrm{Total}}}=$'+str(round(frac2*100,2))+'%',
             transform=ax[1,iax1].transAxes, fontsize=12)
    ax[1,iax1].hist(x=deltaE_array[:,iax1],
                    label='uncorrected', color='blue', **options)
    ax[1,iax1].legend()

radius = dict({'1':'1.3cm', '2':'2.6cm', '3':'5.3cm'})
for ireg in range(nreg):
    plt.text(.02, 1.015, 'Integration cylinder radius: '+radius[str(ireg+1)], transform=ax[0,ireg].transAxes)
plt.xlabel(r'Response (R): $(E_{reco}-E_{gen})/E_{gen}$', fontsize=12)
plt.ylabel('Counts', fontsize=12)
plt.subplots_adjust(wspace=0.2, hspace=0.)
plt.savefig('/eos/user/b/bfontana/www/PartialWafers/'+samples+'/mask'+str(mask)+'/resolution_'+str(mask)+samples+'_'+method+extra_str+'.png')
