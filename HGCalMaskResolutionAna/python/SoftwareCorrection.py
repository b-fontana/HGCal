import os, sys
import numpy as np

from ROOT import TFile, TGraph
from ROOT import Double
from UserCode.HGCalMaskVisualProd.RootUtils import PyDoubleBufferToList as toList
from UserCode.HGCalMaskResolutionAna.PartialWafersStudies import PartialWafersStudies

def DifferentiateShowersByEnergy(current, standard, thresholds, min_val):
        """
        Returns:
        -> 0: the shower is complete.
        -> 1,2,...: each background category.
        """
        assert len(current) == len(standard)
        cumdiff = 0.
        for i in range(len(current)):
            if standard[i] < min_val:
                continue
            cumdiff += abs(standard[i] - current[i])
        if cumdiff < thresholds[0]:
            return 0
        th_shift = np.roll(thresholds, shift=-1)[:-1]
        for it1,it2 in zip(thresholds[:-1],th_shift):
            if cumdiff >= it1 and cumdiff < it2:
                return np.where(th_shift == it2)[0][0] + 1
        if cumdiff >= thresholds[-1]:
            return np.where(th_shift == it2)[0][0] + 2

class IncompleteShowersCorrection(PartialWafersStudies, object):
    """
    Provides all methods required to correct incomplete showers.
    """
    def __init__(self, fname, discrvals):
        """
        Arguments:
        -> fname: file where the weights are stored.
        -> discrvals: discriminatory values used to create background categories.
        """
        super(IncompleteShowersCorrection, self).__init__()
        self.discrvals_ = discrvals
        self.ndiscrintervals_ = len(self.discrvals_)
        self.wn = [['weight1_sr1', 'weight2_sr1', 'weight3_sr1'],
                   ['weight1_sr2', 'weight2_sr2', 'weight3_sr2'],
                   ['weight1_sr3', 'weight2_sr3', 'weight3_sr3']]
        self.hn_sig =  ['en1_per_layer_signal', 
                        'en2_per_layer_signal', 
                        'en3_per_layer_signal']
        self.hn_bckg = [['en1_per_layer_bckg1','en1_per_layer_bckg2','en1_per_layer_bckg3'],
                        ['en2_per_layer_bckg1','en2_per_layer_bckg2','en2_per_layer_bckg3'],
                        ['en3_per_layer_bckg1','en3_per_layer_bckg2','en3_per_layer_bckg3']]
        for ireg in range(self.nsr):
            assert len(self.hn_bckg[ireg]) == self.ndiscrintervals_
        self.f_ = TFile.Open(fname)

        self.weights_ = []
        for ireg in range(self.nsr):
            assert len(self.wn[ireg]) == self.ndiscrintervals_
            self.weights_.append([])
            for iw in range(self.ndiscrintervals_):
                h = self.f_.Get(self.wn[ireg][iw])
                self.weights_[ireg].append([])
                for j in range(1,self.nlayers+1):
                    b = h.FindBin(j)
                    self.weights_[ireg][iw].append(h.GetBinContent(b))

    def CorrectionWeights(self):
        return self.weights_

    def buildCorrectionWeightsGraphs(self, region):
        g = []
        for il in range(self.nlayers):
            g.append(TGraph(self.ndiscrintervals_))
            for iw in range(self.ndiscrintervals_):
                x = Double(self.discrvals_[iw])
                y = Double(self.weights_[region-1][iw][il])
                g[il].SetPoint(iw, x, y)
        return g

    def HistogramsMaxima(self):
        hvals = []
        for i in range(self.nsr):
            tmp = []
            h = self.f_.Get(self.hn_sig[i])
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
        h, corr_graphs = ([] for _ in range(2))
        for ireg in range(self.nsr):
            h.append([])
            corr_graphs.append([])
            for iw in range(self.ndiscrintervals_):
                h[ireg].append(self.f_.Get(self.hn_bckg[ireg][iw]))
                corr_graphs[ireg].append(TGraph(h[ireg][iw].GetNbinsX()))

        for ireg in range(self.nsr):
            graphs = RootHistograms(h[ireg]).toGraph()
            assert len(graphs) == self.ndiscrintervals_
            for iw in range(self.ndiscrintervals_):
                for j in range(graphs[iw].GetN()):
                    x, y = Double(0.), Double(0.) #pass by reference
                    graphs[iw].GetPoint(j,x,y)
                    if x < limits[ireg] and self.weights_[ireg][iw][j] != 0:
                        corr_graphs[ireg][iw].SetPoint(j, x, y/self.weights_[ireg][iw][j])
                    else:
                        corr_graphs[ireg][iw].SetPoint(j, x, y)
            return corr_graphs

    def calculateLowStatisticsFactor(self, limits, mode):
        """
        Arguments:
        -> limits: boundary of the integral
        -> mode: wether to perform the integral on the 'left' or 'right' side
        """
        f = []
        assert len(limits) == len(self.hn_sig) 
        for ireg in range(self.nsr):
            assert len(limits) == len(self.wn[ireg])
            h = self.f_.Get(self.hn_sig[ireg])
            if h.Integral() == 0:
                raise ValueError('The histograms area is zero.')
            if mode == 'left':
                limita = h.FindBin(1) 
                limitb = h.FindBin(limits[ireg])
            elif mode == 'right':
                limita = h.FindBin(limits[ireg])
                limitb = h.FindBin(h.GetNbinsX())
            f.append( h.Integral(limita, limitb) / h.Integral() )
            h.Delete()
        return f
