import os, sys
import numpy as np
from UserCode.HGCalMaskResolutionAna import Argparser
from array import array as Carray
from collections import OrderedDict

from ROOT import TFile, TGraph, TMath
from ROOT import Double
from ROOT import gStyle, gROOT, kTemperatureMap
from UserCode.HGCalMaskVisualProd.RootPlotting import RootPlotting
from UserCode.HGCalMaskVisualProd.RootObjects import RootHistograms, RootGraphs
from UserCode.HGCalMaskVisualProd.RootUtils import PyDoubleBufferToList as toList

parser = Argparser.Argparser()
FLAGS = parser.get_flags()
parser.print_args()

def getCorrectionWeights(fname, hnames, nlayers=28, nregions=3):
    fIn = TFile.Open(fname)
    weights = []
    for i in hnames:
        h = fIn.Get(i)
        weights.append([])
        for j in range(1,nlayers+1):
            b = h.FindBin(j)
            weights[hnames.index(i)].append(h.GetBinContent(b))
    return weights
            
def getHistogramMaxima(hvalues, nlayers=28):
    maxima = []
    for ih in hvalues:
        m = max(ih)
        maxima.append((ih.index(m)+1,m))
    return maxima

def correctIncompleteShowers(fname, hnames, weights, limits):
    """
    Arguments:
    -> fname: Name of the root file where the incomplete energy per layer 
    distributions are stored.
    -> hnames: Name of the histograms to be corrected.
    -> weights: weights to correct the incomplete showers
    -> limits: limits on the x axis where the correction should stop being applied
    """
    assert len(hnames) == len(weights)
    fIn = TFile.Open(fname)
    h = []
    for i in hnames:
        h.append(fIn.Get(i))
    correctedGraphs = [TGraph(h[i].GetNbinsX()) for i in range(len(h))]
    g = RootHistograms(h).toGraph()
    g = RootGraphs(g)
    ilim = -1
    for it,ig in enumerate(g.getObjects()):
        if ilim%3 == 0: 
            ilim += 1
        for j in range(ig.GetN()):
            x, y = Double(0.), Double(0.) #pass by reference
            ig.GetPoint(j,x,y)
            if x <= limits[ilim] and weights[it][j] != 0:
                correctedGraphs[it].SetPoint(j, x, y/weights[it][j])
            else:
                correctedGraphs[it].SetPoint(j, x, y)
    return correctedGraphs

def calculateLowStatisticsFactor(fname, hnames, limits):
    fIn = TFile.Open(fname)
    f = []
    ilim = -1
    for ih in hnames:
        if hnames.index(ih)%3 == 0:
            ilim += 1
        h = fIn.Get(ih)
        if h.Integral() == 0:
            raise ValueError('The histograms area is zero.')
        limitb = h.FindBin(limits[ilim])
        f.append( h.Integral(limitb, h.GetNbinsX()) / h.Integral() )
    return f

def main():
    wn = ['weight1_sr1', 'weight2_sr1', 'weight3_sr1',
          'weight1_sr2', 'weight2_sr2', 'weight3_sr2',
          'weight1_sr3', 'weight2_sr3', 'weight3_sr3']
    hn_sig =  ['en1_per_layer_signal', 'en2_per_layer_signal', 'en3_per_layer_signal']
    hn_bckg = ['en1_per_layer_bckg1', 'en1_per_layer_bckg2', 'en1_per_layer_bckg3',
               'en2_per_layer_bckg1', 'en2_per_layer_bckg2', 'en2_per_layer_bckg3',
               'en3_per_layer_bckg1', 'en3_per_layer_bckg2', 'en3_per_layer_bckg3']
    weights = getCorrectionWeights(FLAGS.noPUFile, hnames=wn)
    fIn = TFile.Open(FLAGS.noPUFile)
    
    hvalues = []
    for i in hn_sig:
        tmp = []
        h = fIn.Get(i)
        for ib in range(h.GetNbinsX()):
            tmp.append(h.GetBinContent(ib))
        hvalues.append(tmp)
    m = getHistogramMaxima(hvalues=hvalues)
    m = [m[x][0] for x in range(len(m))] #y values only
    print("Maxima: ", m)

    correctedGraphs = correctIncompleteShowers(fname=FLAGS.noPUFile, hnames=hn_bckg,
                                               weights=weights, limits=m)
    pcoords = [[[0.01,0.755,0.33,0.995],
                [0.34,0.755,0.66,0.995],
                [0.67,0.755,0.99,0.995],
                [0.01,0.505,0.33,0.745],
                [0.34,0.505,0.66,0.745],
                [0.67,0.505,0.99,0.745],
                [0.01,0.255,0.33,0.495],
                [0.34,0.255,0.66,0.495],
                [0.67,0.255,0.99,0.495],
                [0.01,0.005,0.33,0.245],
                [0.34,0.005,0.66,0.245],
                [0.67,0.005,0.99,0.245]]]
    cdims = [[2000,1600]]
    picname = 'corrected_graphs_'+FLAGS.samples

    #Plotting    
    fIn = TFile.Open(FLAGS.noPUFile)
    extrahistos = []
    extrahistos.append(fIn.Get('en1_per_layer_signal'))
    extrahistos.append(fIn.Get('en2_per_layer_signal'))
    extrahistos.append(fIn.Get('en3_per_layer_signal'))
    extragraphs = RootHistograms(extrahistos).toGraph()

    correctedGraphs = extragraphs + correctedGraphs
    with RootPlotting(ncanvas=1, npads=12, cdims=cdims, pcoords=pcoords) as plot:
        for i,ig in enumerate(correctedGraphs):
            plot.plotGraph(cpos=0, ppos=i, g=ig, lw=2, msize=1.2, mstyle=20, 
                           draw_options='APL')
        plot.save(cpos=0, name=picname)

    factors = calculateLowStatisticsFactor(FLAGS.noPUFile, hn_sig, m)
    print("Factors: ", factors)

if __name__ == "__main__":
    gStyle.SetOptStat(0)
    gROOT.SetBatch(True)
    gStyle.SetPalette(kTemperatureMap)
    main()
