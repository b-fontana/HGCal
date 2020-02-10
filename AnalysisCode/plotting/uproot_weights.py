import sys
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams['font.size'] = 20
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib import rcParams as rc
#rc.update({'font.size': 12})
#import seaborn as sns
from uproot import open as upopen
import numpy as np

nreg = 3

##############################################################################
########Full plot for the Detector's Note#####################################
#########Run with 'python uproot_weights 3 outer' ############################
##############################################################################
def full_weights_plot():
    mask = sys.argv[1]
    samples = sys.argv[2]
    file = upopen("root_files/fileweights_"+str(mask)+samples+".root")
    keys = file.keys()

    figsize = (15.,12.)
    figscale = figsize[0] / figsize[1]
    figratios = [1,.5]
    fig, ax = plt.subplots(len(figratios), nreg, figsize=figsize, 
                           sharex=True, gridspec_kw={'height_ratios': figratios}, 
                           sharey=False, squeeze=False)
    fig.add_subplot(111, frameon=False)
    plt.subplots_adjust(wspace=0.2, hspace=0.)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.grid(False)
    colors = ['blue', 'orange','green','brown']
    #energy distributions
    for ireg in range(nreg):
        counts, bins = file['en{}_layer_sign;1'.format(ireg+1)].numpy()
        ax[0,ireg].hist(bins[:-1], bins, weights=counts, histtype='step', label='signal', color=colors[0])
        for iw in range(3):
            counts, bins = file['en{}_layers_bckg{};1'.format(ireg+1,iw+1)].numpy()
            ax[0,ireg].hist(bins[:-1], bins, weights=counts, histtype='step', label='bckg'+str(iw+1), color=colors[iw+1])
        ax[0,ireg].legend()
    #weights
    for ireg in range(nreg):
        for iw in range(3):
            counts, bins = file['weight{}_sr{};1'.format(iw+1,ireg+1)].numpy()
            ax[1,ireg].hist(bins[:-1], bins, weights=counts, histtype='step', label='weight'+str(iw+1), color=colors[iw+1])
            ax[1,ireg].legend(loc='upper right')
    radius = dict({'1':'1.3cm', '2':'2.6cm', '3':'5.3cm'})
    for ireg in range(nreg):
        plt.text(.05, 1.015, 'Integration cylinder radius: '+radius[str(ireg+1)], transform=ax[0,ireg].transAxes)
        ax[0,ireg].set_yticks([0.025, 0.05, 0.075, 0.1, 0.125])
    ax[1,1].set_xlabel('Layers', fontsize=str(12*figscale), labelpad=1.7*figscale)
    ax[0,0].set_ylabel('Fraction of deposited energy', fontsize=str(12*figscale), labelpad=1.7*figscale)
    ax[1,0].set_ylabel('Weights', fontsize=str(12*figscale), labelpad=1.7*figscale)
    plt.savefig('/eos/user/b/bfontana/www/PartialWafers/'+samples+'/mask'+str(mask)+'/weights_'+str(mask)+samples+'.png')

##############################################################################
########Partial plot for the Posters@LHCC#####################################
########Run with 'python uproot_weights' #####################################
##############################################################################
def poster_weight_plot():
    mask = str(3)
    samples = ('inner', 'outer')
    files = (upopen("root_files/fileweights_"+mask+samples[0]+".root"),
             upopen("root_files/fileweights_"+mask+samples[1]+".root"))
    keys = (files[0].keys(), files[0].keys())

    figsize = (11.,20.)
    figscale = figsize[0] / figsize[1]
    fig, ax = plt.subplots(2, 1, figsize=figsize, sharex=True, sharey=False, squeeze=False)
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.grid(False)
    legprops = {'size': 20}
    colors = ['darkblue', 'mediumturquoise', 'darkred', 'tomato']

    #inner region energy distributions
    counts, bins = ( files[0]['en3_layer_sign;1'] ).numpy()
    ax[0,0].hist(bins[:-1], bins, weights=counts, histtype='stepfilled', label='Inner region, fully contained', color=colors[0])
    counts, bins = ( files[0]['en3_layers_bckg3;1'] ).numpy()
    ax[0,0].hist(bins[:-1], bins, weights=counts, histtype='stepfilled', label='Inner region, at the boundary', color=colors[1])
    ax[0,0].legend(prop=legprops)
    #outer region energy distributions
    counts, bins = ( files[1]['en3_layer_sign;1'] ).numpy()
    ax[1,0].hist(bins[:-1], bins, weights=counts, histtype='stepfilled', label='Outer region, fully contained', color=colors[2])
    counts, bins = ( files[1]['en3_layers_bckg3;1'] ).numpy()
    ax[1,0].hist(bins[:-1], bins, weights=counts, histtype='stepfilled', label='Outer region, at the boundary', color=colors[3])
    ax[1,0].legend(prop=legprops)
    #text
    plt.text(.01, 1.015, 'Mask '+mask, transform=ax[0,0].transAxes, fontsize=20)
    plt.text(.44, 1.015, 'Integration cylinder radius: 5.3cm', transform=ax[0,0].transAxes, fontsize=20)
    ax[0,0].set_yticks([0.025, 0.05, 0.075, 0.1, 0.125])
    fig.text(0.45, 0.08, 'Layers', fontsize=20, rotation='horizontal')
    fig.text(0.015, 0.64, 'Fraction of deposited energy', fontsize=20, rotation='vertical')
    plt.subplots_adjust(hspace=0.001)
    plt.savefig('/eos/user/b/bfontana/www/PartialWafers/longitudinal_profiles_selection.png')

##################################################################################
######################### MAIN ###################################################
##################################################################################
poster_weight_plot()
