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

def DifferentiateShowers(current, standard):




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

    for i in range(0,data.GetEntriesFast()):
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
                if ( (FLAGS.samples == "inner" and geneta < FLAGS.maxgeneta or
                     FLAGS.samples == "outer" and geneta > FLAGS.mingeneta) and
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
    etabins = 12
    etacuts = FLAGS.etacuts
    hn = ['fraction_complete_sr{}']
    for ireg in range(1,NREG+1):
            bins = Carray('d', np.arange(2.69, 3.04, 0.01))
            strings = ';#Eta;Fraction'
            for ih in range(len(hn)):
                histos[hn[ih].format(ireg)] = TH1F(hn[ih].format(ireg), strings,
                                                   len(bins)-1, bins)
    for h in histos:
        histos[h].Sumw2()
        histos[h].SetMarkerStyle(20)
        histos[h].SetDirectory(0)

    for i in range(0,data.GetEntriesFast()):
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
            ROI_en = np.zeros((NLAYERS))
            for il in range(1,NLAYERS+1):
                v = f1*getattr(data,'en_sr{}_layer{}'.format(ireg,il)) - f2
                ROI_en.append(v/recen)
            DifferentiateShowers(ROI_en, fracEn[ireg-1.il-1])

            quit()
            bool_sig = ( (FLAGS.samples == 'inner' and geneta < FLAGS.maxgeneta) or 
                         (FLAGS.samples == 'outer' and geneta >= FLAGS.mingeneta) )
            bools_bckg = []
            check_bools = int(bool_sig)
            for ic in range(len(etacuts)-1):
                check_bools += int(bools_bckg[ic])
            assert check_bools == 1
                
    #end of tree loop    
    fIn.Close()

    """
    pcoords = [[[0.01,0.51,0.33,0.99],
                [0.34,0.51,0.66,0.99],  
                [0.67,0.51,0.99,0.99],
                [0.01,0.01,0.33,0.49],   
                [0.34,0.01,0.66,0.49],   
                [0.67,0.01,0.99,0.49]]]
    cdims = [[1000,600]]
    picname = '2comp_'+FLAGS.samples
    plotHistograms(histos, cdims, pcoords, 
                   os.path.join(FLAGS.outpath,picname))
    """

if __name__ == "__main__":
    parser = Argparser.Argparser()
    FLAGS = parser.get_flags()
    parser.print_args()
    base = PartialWafersStudies()
    NREG, NLAYERS, A = base.nsr, base.nlayers, base.sr_area
    main()
