import abc
from array import array as Carray
from ROOT import TFile, TCanvas, TPad, TMath, TStyle, TExec, TLatex, TGraph, THStack
from ROOT import gSystem, gDirectory, gStyle, kBlue

from RootUtils import PyDoubleBufferToList as toList

class _RootObjects:
    """
    Base class.
    """
    __metaclass__ = abc.ABCMeta
    def __init__(self, objects):
        if not isinstance(objects, list):
            raise TypeError('The Root objects must be stored as a list.')
        self.o = objects
        self.N = len(self.o)
        
    def save(self, fname, mode='RECREATE'):
        f = TFile.Open(fname+'.root', mode)
        for io in range(len(self.o)):
            self.o[io].Write()
        f.Write()

class RootHistograms(_RootObjects):
    """
    Root util functions for 1D histograms.
    """
    def __init__(self, *args, **kwargs):
        super(RootHistograms, self).__init__(*args, **kwargs)
        self.hstack = {}

    def toGraph(self):
        graphs = []
        for o in self.o:
            xval, yval = ([] for _ in range(2))
            for ib in range(1,o.GetNbinsX()+1):
                xval.append(ib)
                yval.append(o.GetBinContent( o.FindBin(ib) ))
            assert len(xval) == len(yval)
            graphs.append(TGraph(i, Carray('d', xval[:]), Carray('d', yval[:])))
        return graphs

    def stack(self, name, idx=[None], title=''):
        if title == '':
            self.hstack[name] = THStack(name, name)
        else:
            self.hstack[name] = THStack(name, title)
        if idx[0] == None:
            for o in self.o:
                self.hstack[name].Add(o)
        else:
            for i in idx:
                self.hstack[name].Add(self.o[i])
        return self.hstack[name]

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
        for i in range(len(self.o)):
            xvals, yvals = toList(self.o[i].GetX()), toList(self.o[i].GetY())
            xmax, ymax = xvals[0], yvals[0]
            xvals, yvals = xvals[1:], yvals[1:]
            for x,y in zip(xvals,yvals):
                if y > ymax: 
                    xmax, ymax = x, y
            maxs.append((xmax,ymax))
        return maxs
