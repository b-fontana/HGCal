import os, sys
import numpy as np
from UserCode.HGCalMaskResolutionAna import Argparser
from array import array as Carray
from collections import OrderedDict

from ROOT import TCanvas, TLatex, TFile, TMath, TH1F
from ROOT import TLegend, TH2F, TLorentzVector, TProfile, TH1D
from ROOT import gStyle, gROOT, kTemperatureMap
from UserCode.HGCalMaskVisualProd.RootUtils import RootPlotting, RootHistograms

parser = Argparser.Argparser()
FLAGS = parser.get_flags()
parser.print_args()

def getCorrectionWeights(fname, nlayers=28):
    fIn = TFile.Open(fname)
    h = []
    hn = ['weight1_sr0', 'weight2_sr0', 'weight3_sr0',
          'weight1_sr1', 'weight2_sr1', 'weight3_sr1',
          'weight1_sr2', 'weight2_sr2', 'weight3_sr2']
    for i in hn:
        h.append(fIn.Get(i))
    g = RootHistograms(h).toGraph()
    
    pcoords = [[[0.01,0.67,0.33,0.99],   #canvas0, pad0
                [0.34,0.67,0.66,0.99],   #canvas0, pad1
                [0.67,0.67,0.99,0.99],
                [0.01,0.34,0.33,0.66],   #canvas0, pad0
                [0.34,0.34,0.66,0.66],   #canvas0, pad1
                [0.67,0.34,0.99,0.66],   #canvas0, pad2
                [0.01,0.01,0.33,0.33],   #canvas0, pad3
                [0.34,0.01,0.66,0.33],   #canvas0, pad4
                [0.67,0.01,0.99,0.33]]]  #canvas0, pad5
    cdims = [[2000,1200]]
    picname = 'graphs_'+FLAGS.samples

    with RootPlotting(ncanvas=1, npads=len(h), cdims=cdims, pcoords=pcoords) as plot:
        for i,ig in enumerate(g):
            plot.plotGraph(cpos=0, ppos=i, g=ig, title=hn[i], 
                           lw=2, msize=1.2, mstyle=20, draw_options='APL')
        plot.save(cpos=0, name=picname)


def main():
    getCorrectionWeights(FLAGS.noPUFile+'.root')

if __name__ == "__main__":
    gStyle.SetOptStat(0)
    gROOT.SetBatch(True)
    gStyle.SetPalette(kTemperatureMap)
    main()
