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

class IncompleteShowersCorrection:
    """
    Provides all methods required to correct incomplete showers.
    """
    def __init__(self, fname):
        self.nlayers_ = 28
        self.nreg_ = 3
        self.wn = ['weight1_sr1', 'weight2_sr1', 'weight3_sr1',
                   'weight1_sr2', 'weight2_sr2', 'weight3_sr2',
                   'weight1_sr3', 'weight2_sr3', 'weight3_sr3']
        self.hn_sig =  ['en1_per_layer_signal', 
                        'en2_per_layer_signal', 
                        'en3_per_layer_signal']
        self.hn_bckg = ['en1_per_layer_bckg1', 'en1_per_layer_bckg2', 'en1_per_layer_bckg3',
                        'en2_per_layer_bckg1', 'en2_per_layer_bckg2', 'en2_per_layer_bckg3',
                        'en3_per_layer_bckg1', 'en3_per_layer_bckg2', 'en3_per_layer_bckg3']
        self.f_ = TFile.Open(fname)
        self.weights_ = self.getCorrectionWeights()

        factors = calculateLowStatisticsFactor(FLAGS.noPUFile+'.root', hn_sig, m)

    def getCorrectionWeights(self):
        weights = []
        for i in self.wn:
            h = self.f_.Get(i)
            weights.append([])
            for j in range(1,self.nlayers_+1):
                b = h.FindBin(j)
                weights[self.wn.index(i)].append(h.GetBinContent(b))
        return weights
            
    def getHistogramsMaxima(self):
        hvals = []
        for i in self.hn_sig:
            tmp = []
            h = self.f_.Get(i)
            for ib in range(h.GetNbinsX()):
                tmp.append(h.GetBinContent(ib))
            hvals.append(tmp)
        maxima = []
        for ih in hvals:
            m = max(ih)
            maxima.append((ih.index(m)+1,m))
        return [maxima[x][0] for x in range(len(maxima))] #y values only

    def correctIncompleteShowers(self, limits):
        """
        Arguments:
        -> fname: Name of the root file where the incomplete energy per layer 
        distributions are stored.
        -> hnames: Name of the histograms to be corrected.
        -> weights: weights to correct the incomplete showers
        -> limits: limits on the x axis where the correction should stop being applied
        """
        assert len(self.hn_bckg) == len(self.weights)
        h = []
        for i in hnames:
            h.append(self.f_.Get(i))
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
                if x <= limits[ilim] and self.weights[it][j] != 0:
                    correctedGraphs[it].SetPoint(j, x, y/self.weights[it][j])
                else:
                    correctedGraphs[it].SetPoint(j, x, y)
        return correctedGraphs

    def calculateLowStatisticsFactor(self, limits):
        assert len(limits) == len(self.weights) / self.nreg_
        assert len(limits) == len(self.hn_sig)
        f = []
        ilim = -1
        for ih in self.hn_sig:
            if self.hn_sig.index(ih)%3 == 0:
                ilim += 1
            h = self.f_.Get(ih)
            if h.Integral() == 0:
                raise ValueError('The histograms area is zero.')
            limitb = h.FindBin(limits[ilim])
            f.append( h.Integral(limitb, h.GetNbinsX()) / h.Integral() )
        return f
    
