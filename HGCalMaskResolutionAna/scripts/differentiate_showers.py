import os, sys
import numpy as np
import pickle
from UserCode.HGCalMaskResolutionAna import Argparser
from array import array as Carray
from collections import OrderedDict

from ROOT import TCanvas, TLatex, TFile, TMath, TH1F
from ROOT import TLegend, TH2F, TLorentzVector, TProfile, TH1D
from ROOT import gStyle, gROOT, kTemperatureMap
from UserCode.HGCalMaskVisualProd.SystemUtils import averageContiguousVals as Av
from UserCode.HGCalMaskResolutionAna.SoftwareCorrection import IncompleteShowersCorrection
from UserCode.HGCalMaskResolutionAna.PartialWafersStudies import PartialWafersStudies
from UserCode.HGCalMaskVisualProd.RootPlotting import RootPlotting
from UserCode.HGCalMaskVisualProd.RootObjects import RootHistograms

def DifferentiateShowers(current, standard, threshold, min_val):
    """Returns True if the shower is complete."""
    assert len(current) == len(standard)
    cumdiff = 0.
    for i in range(len(current)):
        if standard[i] < min_val:
            continue
        cumdiff += abs(standard[i] - current[i])
    if cumdiff > threshold:
        return False
    return True

def main():
    #gStyle.SetOptStat(0)
    gROOT.SetBatch(True)
    gStyle.SetPalette(kTemperatureMap)

    fIn=TFile.Open(FLAGS.noPUFile)
    data=fIn.Get('data')

    with open('calib_'+FLAGS.samples+'_nopu.pck','r') as cachefile:
        calib = pickle.load(cachefile)

    fracEn = np.zeros((NREG,NLAYERS), dtype=float)
    countfracEn = np.zeros((NREG,NLAYERS), dtype=int)
    nentries = data.GetEntriesFast()
    for i in range(0,nentries):
        data.GetEntry(i)
        genen    = getattr(data,'genen')        
        geneta   = abs(getattr(data,'geneta'))
        genphi   = getattr(data,'genphi')

        for ireg in range(1,NREG+1):
            recen = getattr(data,'en_sr{}_ROI'.format(ireg))
            avgnoise = getattr(data,'noise_sr3_ROI')*A[ireg-1]/A[2]

            #Calibration factors. f2 is used for PU.
            f1, f2 = 1., 0.
            if 'L0' in calib:
                f1 /= calib['L0'][ireg].Eval(geneta)+1.0
                if 'L1' in calib:
                    f1 /= calib['L1'][ireg].Eval(f1*recen)+1.0    
                    if 'L2' in calib and ireg in calib['L2']:
                        f2 = calib['L2'][ireg].Eval(avgnoise)
            recen = f1*recen - f2 

            for il in range(1,NLAYERS+1):
                v = f1*getattr(data,'en_sr{}_layer{}'.format(ireg,il)) - f2
                if ( (FLAGS.samples == "inner" and geneta < 2.75 or
                      FLAGS.samples == "outer" and geneta > 1.6) and
                     recen != 0 ):
                    fracEn[ireg-1,il-1] += v / recen
                    countfracEn[ireg-1,il-1] += 1

    fracEn /= countfracEn
    print(np.sum(fracEn, axis=1))
    """
    #Write energy distribution of complete showers for complete/incomplete discrimination
    bins = Carray('d', np.arange(0.5,NLAYERS+1,1))
    for ireg in range(NREG):
        h = TH1F('complete_showers_standard_sr'+str(ireg+1), 
                 'complete_showers_standard_sr'+str(ireg+1),
                 len(bins)-1, bins)
        for il in range(NLAYERS):
            b = h.FindBin(il+1)
            h.Fill(b,fracEn[ireg,il])
        h.Write()
    """

    histos=OrderedDict()
    hn = ['fraction_complete_sr{}', 'complete_sr{}']
    for ireg in range(1,NREG+1):
            bins = Carray('d', np.arange(2.69, 3.04, 0.001))
            strings = ";|#eta|;Complete showers' fraction"
            for ih in range(len(hn)):
                histos[hn[ih].format(ireg)] = TH1F(hn[ih].format(ireg), strings,
                                                   len(bins)-1, bins)
    for h in histos:
        histos[h].Sumw2()
        histos[h].SetMarkerStyle(20)
        histos[h].SetDirectory(0)

    for i in range(0,nentries):
        data.GetEntry(i)
        genen    = getattr(data,'genen')        
        geneta   = abs(getattr(data,'geneta'))
        genphi   = getattr(data,'genphi')

        for ireg in range(1,NREG+1):
            recen = getattr(data,'en_sr{}_ROI'.format(ireg))

            #Calibration factors. f2 is used for PU.
            f1, f2 = 1., 0.
            if 'L0' in calib:
                f1 /= calib['L0'][ireg].Eval(geneta)+1.0
                if 'L1' in calib:
                    f1 /= calib['L1'][ireg].Eval(f1*recen)+1.0    
                    if 'L2' in calib and ireg in calib['L2']:
                        f2 = calib['L2'][ireg].Eval(avgnoise)
            recen = f1*recen - f2 
            #differentiate complete from incomplete showers
            ROI_en = np.zeros((NLAYERS), dtype=float)
            for il in range(1,NLAYERS+1):
                v = f1*getattr(data,'en_sr{}_layer{}'.format(ireg,il)) - f2
                try:
                    ROI_en[il-1] = v/recen
                except ZeroDivisionError:
                    ROI_en[il-1] = 0.
            bool_sig = DifferentiateShowers(ROI_en, fracEn[ireg-1,:], 
                                            threshold=.5, min_val=0.05)
            if bool_sig:
                histos[hn[0].format(ireg)].Fill(geneta)
            histos[hn[1].format(ireg)].Fill(geneta)
 
    #end of tree loop    
    fIn.Close()

    pcoords = [[[0.01,0.01,0.33,0.99],
                [0.34,0.01,0.66,0.99],
                [0.67,0.01,0.99,0.99]]]
    cdims = [[2000,700]]
    picname = 'diffshowers'    
    correct_order = []
    for i in range(len(hn)):
        for ireg in range(1,NREG+1):
            correct_order.append(hn[i].format(ireg))
    assert len(correct_order) == len(histos.keys())
    histos = [histos[correct_order[i]] for i in range(len(correct_order))]

    for ireg in range(NREG):
        histos[ireg].Divide(histos[ireg+NREG])
    histos = histos[:3]
    
    with RootPlotting(ncanvas=1, npads=3, cdims=cdims, pcoords=pcoords) as plot:
        for ireg in range(NREG):
            plot.plotHistogram(cpos=0, ppos=ireg, h=histos[ireg], draw_options='HIST')
        plot.save(cpos=0, name=picname)


if __name__ == "__main__":
    parser = Argparser.Argparser()
    FLAGS = parser.get_flags()
    parser.print_args()
    base = PartialWafersStudies()
    NREG, NLAYERS, A = base.nsr, base.nlayers, base.sr_area
    main()
