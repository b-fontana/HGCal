import abc
from array import array as Carray
from ROOT import TFile, TCanvas, TPad, TMath, TStyle, TExec, TLatex, TGraph
from ROOT import gSystem, gDirectory, gStyle, kBlue

from RootUtils import PyDoubleBufferToList as toList

class _RootObjects:
    """
    Base class.
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, objects):
        self.o_ = objects
        
    def getObject(self, idx):
        return self.o_[idx]

    def save(self, fname):
        f = TFile.Open(fname+'.root', mode)
        for io in range(len(self.o_)):
            self.o_[io].Write()
        f.Write()

class RootHistograms(_RootObjects):
    """
    Root util functions for 1D histograms.
    """
    def __init__(self, *args, **kwargs):
        super(RootHistograms, self).__init__(*args, **kwargs)

    def toGraph(self):
        graphs = []
        for o in self.o_:
            xval, yval = ([] for _ in range(2))
            for ib in range(1,o.GetNbinsX()+1):
                xval.append(ib)
                yval.append(o.GetBinContent( o.FindBin(ib) ))
            assert len(xval) == len(yval)
            graphs.append(TGraph(ib, Carray('d', xval[:]), Carray('d', yval[:])))
        return graphs

class RootGraphs(_RootObjects):
    def __init__(self, *args, **kwargs):
        """
        Root util functions for graphs.
        """
        super(RootGraphs, self).__init__(*args, **kwargs)
    
    def getMax(self):
        """
        Returns the x and y coordinates of the maxima of the graph.
        It only captures one maximum.
        """
        maxs = []
        for i in range(len(self.o_)):
            xvals, yvals = toList(self.o_[i].GetX()), toList(self.o_[i].GetY())
            xmax, ymax = xvals[0], yvals[0]
            xvals, yvals = xvals[1:], yvals[1:]
            for x,y in zip(xvals,yvals):
                if y > ymax: 
                    xmax, ymax = x, y
            maxs.append((xmax,ymax))
        return maxs
