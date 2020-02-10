import sys
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.optimize import curve_fit, OptimizeWarning
import warnings
warnings.simplefilter('error', OptimizeWarning)
import uproot as up
import params

extra_str = ''
nreg = 3

###############################################################################
#######plot 2 distributions and perform iterative fitting procedure############
###############################################################################
def detector_note_function():
    mask = sys.argv[1]
    samples = sys.argv[2]
    method = sys.argv[3]
    correction = False if int(sys.argv[4])==0 else True
    correction_str = '_corr_'+method if correction else '_nocorr'
    file = up.open('root_files/final_'+str(mask)+samples+'_'+method+'.root')
    #trees = file.allkeys(filterclass=lambda x: issubclass(x, up.tree.TTreeMethods))
    tree = file['data']
    geneta = tree.array('abs_geneta')
    deltaE_array, deltaEcorr_array = tree.array('deltaE'), tree.array('deltaE_corr')

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
#######plot the superimposed resolutions for the first two signal regions######
###############################################################################
def poster_plot():
    masks = (str(3), str(6)) 
    samples = 'inner'
    method = 'fineeta'
    files = (up.open('root_files/final_'+masks[0]+samples+'_'+method+'.root'),
             up.open('root_files/final_'+masks[1]+samples+'_'+method+'.root'))
    trees = (files[0]['data'], files[1]['data'])
    genetas = (trees[0].array('abs_geneta'), trees[1].array('abs_geneta'))
    sr = 2 #signal region with 5.3 cm of radius
    deltaE_1, deltaEcorr_1 = (trees[0].array('deltaE'))[:,sr], (trees[0].array('deltaE_corr'))[:,sr]
    deltaE_2, deltaEcorr_2 = (trees[1].array('deltaE'))[:,sr], (trees[1].array('deltaE_corr'))[:,sr]

    textx = .71
    texty = .57, .72
    textsize = 13
    figsize = (18,8)
    width, height = 2, 2
    xrangetuple = (-.98,.5)
    etas = (2.92, 2.93)
    alpha = 0.4
    fig, ax = plt.subplots(height, width, figsize=figsize, sharex=True, sharey=False, squeeze=False)
    fig.add_subplot(111, frameon=False)
    plt.subplots_adjust(wspace=0.05, hspace=0.)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.grid(False)
    options = dict(bins=(100), range=xrangetuple, histtype='step', density=True)
    colors = ('darkblue', 'mediumturquoise', 'darkred', 'tomato')

    #full eta range
    frac_tot_corr_1 = float(deltaEcorr_1[deltaEcorr_1<-0.5].size) / deltaEcorr_1.size
    frac_tot_uncorr_1 = float(deltaE_1[deltaE_1<-0.5].size) / deltaE_1.size
    plt.text(textx, texty[0], r'$\frac{N_{\mathrm{corrected}<-0.5}}{N_{\mathrm{Total}}}=$'+str(round(frac_tot_corr_1*100,2))+'%', transform=ax[0,0].transAxes, fontsize=textsize)
    plt.text(textx, texty[1], r'$\frac{N_{\mathrm{uncorrected}<-0.5}}{N_{\mathrm{Total}}}=$'+str(round(frac_tot_uncorr_1*100,2))+'%', transform=ax[0,0].transAxes, fontsize=textsize)
    frac_tot_corr_2 = float(deltaEcorr_2[deltaEcorr_2<-0.5].size) / deltaEcorr_2.size
    frac_tot_uncorr_2 = float(deltaE_2[deltaE_2<-0.5].size) / deltaE_2.size
    plt.text(textx, texty[0], r'$\frac{N_{\mathrm{corrected}<-0.5}}{N_{\mathrm{Total}}}=$'+str(round(frac_tot_corr_2*100,2))+'%', transform=ax[0,1].transAxes, fontsize=textsize)
    plt.text(textx, texty[1], r'$\frac{N_{\mathrm{uncorrected}<-0.5}}{N_{\mathrm{Total}}}=$'+str(round(frac_tot_uncorr_2*100,2))+'%', transform=ax[0,1].transAxes, fontsize=textsize)
    ax[0,0].hist(x=deltaE_1, label='Mask '+masks[0]+' uncorrected', color=colors[0], **options)
    ax[0,0].hist(x=deltaEcorr_1, label='Mask '+masks[0]+' corrected', color=colors[1], **options)
    ax[0,0].legend()
    ax[0,1].hist(x=deltaE_2, label='Mask '+masks[1]+' uncorrected', color=colors[2], **options)
    ax[0,1].hist(x=deltaEcorr_2, label='Mask '+masks[1]+' corrected', color=colors[3], **options)
    ax[0,1].legend()
    plt.text(0.02, .88, r'Full $|\eta|$ region', bbox=dict(facecolor='blue', alpha=alpha), transform=ax[0,0].transAxes, fontsize=textsize)
    plt.text(0.02, .88, r'$'+str(etas[0])+' < |\eta| < '+str(etas[1])+'$', bbox=dict(facecolor='blue', alpha=alpha), transform=ax[1,0].transAxes, fontsize=textsize)
    plt.text(0.01, 1.015, r'Correction: Fine calibration', transform=ax[0,0].transAxes, fontsize=textsize)
    plt.text(0.52, 1.015, r'Signal region: 5.3cm; Inner region', transform=ax[0,1].transAxes, fontsize=textsize)
    
    #specific eta range
    uncorr_2d_1 = np.array([np.array(deltaE_1), np.array(genetas[0])])
    corr_2d_1 = np.array([np.array(deltaEcorr_1), np.array(genetas[1])])
    uncorr_2d_1 = uncorr_2d_1[0][ (uncorr_2d_1[1]>etas[0]) & (uncorr_2d_1[1]<etas[1]) ]
    corr_2d_1 = corr_2d_1[0][ (corr_2d_1[1]>etas[0]) & (corr_2d_1[1]<etas[1]) ]

    uncorr_2d_2 = np.array([np.array(deltaE_2), np.array(genetas[0])])
    corr_2d_2 = np.array([np.array(deltaEcorr_2), np.array(genetas[1])])
    uncorr_2d_2 = uncorr_2d_2[0][ (uncorr_2d_2[1]>etas[0]) & (uncorr_2d_2[1]<etas[1]) ]
    corr_2d_2 = corr_2d_2[0][ (corr_2d_2[1]>etas[0]) & (corr_2d_2[1]<etas[1]) ]

    frac_part_uncorr_1 = float(uncorr_2d_1[uncorr_2d_1<-0.5].size) / uncorr_2d_1.size
    frac_part_corr_1 = float(corr_2d_1[corr_2d_1<-0.5].size) / corr_2d_1.size
    plt.text(textx, texty[0], r'$\frac{N_{\mathrm{corrected}<-0.5}}{N_{\mathrm{Total}}}=$'+str(round(frac_part_uncorr_1*100,2))+'%', transform=ax[1,0].transAxes, fontsize=textsize)
    plt.text(textx, texty[1], r'$\frac{N_{\mathrm{corrected}<-0.5}}{N_{\mathrm{Total}}}=$'+str(round(frac_part_corr_1*100,2))+'%', transform=ax[1,0].transAxes, fontsize=textsize)
    frac_part_uncorr_2 = float(uncorr_2d_2[uncorr_2d_2<-0.5].size) / uncorr_2d_2.size
    frac_part_corr_2 = float(corr_2d_2[corr_2d_2<-0.5].size) / corr_2d_2.size
    plt.text(textx, texty[0], r'$\frac{N_{\mathrm{uncorrected}<-0.5}}{N_{\mathrm{Total}}}=$'+str(round(frac_part_uncorr_2*100,2))+'%', transform=ax[1,1].transAxes, fontsize=textsize)
    plt.text(textx, texty[1], r'$\frac{N_{\mathrm{corrected}<-0.5}}{N_{\mathrm{Total}}}=$'+str(round(frac_part_corr_2*100,2))+'%', transform=ax[1,1].transAxes, fontsize=textsize)
    ax[1,0].hist(x=uncorr_2d_1, label='Mask '+masks[0]+' uncorrected', color=colors[0], **options)
    ax[1,0].hist(x=corr_2d_1, label='Mask '+masks[0]+' corrected', color=colors[1], **options)
    ax[1,0].legend()
    ax[1,1].hist(x=uncorr_2d_2, label='Mask '+masks[1]+' uncorrected', color=colors[2], **options)
    ax[1,1].hist(x=corr_2d_2, label='Mask '+masks[1]+' corrected', color=colors[3], **options)
    ax[1,1].legend()
    plt.text(0.75, .35, r'Full $|\eta|$ region', bbox=dict(facecolor='red', alpha=alpha), transform=ax[0,1].transAxes, fontsize=textsize)
    plt.text(0.75, .35, r'$'+str(etas[0])+' < |\eta| < '+str(etas[1])+'$', bbox=dict(facecolor='red', alpha=alpha), transform=ax[1,1].transAxes, fontsize=textsize)

    ax[0,0].yaxis.set_ticklabels([])
    ax[0,1].yaxis.set_ticklabels([])
    ax[1,0].yaxis.set_ticklabels([])
    ax[1,1].yaxis.set_ticklabels([])

    quart25 = np.quantile(corr_2d_1, .25)
    quart50 = np.quantile(corr_2d_1, .50)
    quart75 = np.quantile(corr_2d_1, .75)
    quartdist = 1.5*abs(quart75 - quart25)
    kwargs ={'width': 0.004, 'head_width': 0.05, 'head_length': 0.55, 'head_starts_at_zero': False, 'color':'mediumturquoise'}
    """
    ax[1,0].arrow(quart25, 8.5, 0, -4, **kwargs)
    ax[1,0].arrow(quart75, 8.5, 0, -4, **kwargs)
    """
    patchkwargs = {'linewidth': 2., 'alpha': 1.}
    linekwargs = {'linewidth': 1, 'alpha': 1.}
    rect1 = patches.Rectangle(xy=(quart25,0.5), width=abs(quart50-quart25), height=5, edgecolor='royalblue', facecolor='None', **patchkwargs)
    rect2 = patches.Rectangle((quart50,0.5),abs(quart75-quart50),5,edgecolor='royalblue',facecolor='None',**patchkwargs)
    lineleft = patches.Rectangle(xy=(quart25-quartdist,2.75), width=quartdist, height=0.06, edgecolor='royalblue', facecolor='royalblue', **linekwargs)
    lineright = patches.Rectangle(xy=(quart75,2.75), width=quartdist, height=0.06, edgecolor='royalblue', facecolor='royalblue', **linekwargs)
    wiskleft = patches.Rectangle(xy=(quart25-quartdist,2.25), width=0.007, height=1, edgecolor='royalblue', facecolor='royalblue', **linekwargs)
    wiskright = patches.Rectangle(xy=(quart75+quartdist,2.25), width=0.007, height=1, edgecolor='royalblue', facecolor='royalblue', **linekwargs)
    ax[1,0].add_patch(rect1)
    ax[1,0].add_patch(rect2)
    ax[1,0].add_patch(lineleft)
    ax[1,0].add_patch(lineright)
    ax[1,0].add_patch(wiskleft)
    ax[1,0].add_patch(wiskright)

    plt.xlabel(r'Response: $(E_{reco}-E_{gen})/E_{gen}$', fontsize=12)
    plt.ylabel('Normalized counts', fontsize=textsize)
    plt.savefig('/eos/user/b/bfontana/www/PartialWafers/'+samples+'/resolution_overlayed_'+samples+'_'+method+extra_str+'.png')


###############################################################################
#######main####################################################################
###############################################################################
poster_plot()
