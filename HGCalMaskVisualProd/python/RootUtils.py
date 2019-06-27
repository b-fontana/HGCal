from array import array as Carray
from ROOT import TFile, TCanvas, TPad, TMath, TStyle, TExec, TLatex, TGraph
from ROOT import gSystem, gDirectory, gStyle, kBlue

def PyDoubleBufferToList(b):
    n = len(b)
    l = [b[i] for i in range(n)]
    return l
