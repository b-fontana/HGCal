import sys
import ROOT
import array
import numpy as np
from UserCode.HGCalMaskResolutionAna import Argparser
from array import array as Carray

from ROOT import TFile, TNtuple, TH1F
from UserCode.HGCalMaskResolutionAna.ROISummary import ROISummary
from UserCode.HGCalMaskResolutionAna.PartialWafersStudies import PartialWafersStudies

parser = Argparser.Argparser()
FLAGS = parser.get_flags()
parser.print_args()

def main():
    """
    Creates a simple summary of the ROI for fast calibration.
    """
    fOut = TFile(FLAGS.outpath, 'RECREATE')
    varnames  = ['genen', 'geneta', 'genphi']
    for i in range(1,NREG+1):
        varnames += ['en_sr{}_ROI'.format(i), 
                     'noise_sr{}_ROI'.format(i)]
        for il in range(1,NLAYERS+1):        
            varnames += ['en_sr{}_layer{}'.format(i,il), 
                         'noise_sr{}_layer{}'.format(i,il)]
    output_tuple = TNtuple("data","data",":".join(varnames))

    fracEn = np.zeros((NREG,NLAYERS), dtype=float)
    countfracEn = np.zeros((NREG,NLAYERS), dtype=int)

    fIn = TFile.Open(FLAGS.noPUFile)
    t = fIn.Get('an_mask/data')
    for i in range(0,t.GetEntriesFast()):
        t.GetEntry(i)

        #define the ROIs
        roiList={}
        for ir in range(0,t.ROIs.size()):
            if t.ROIs[ir].id() < 0: 
                continue
            roiList[ir] = ROISummary(t.ROIs[ir].p4())

        for h in t.Hits:
            roiIdx  = h.associatedROI()
            rid     = t.ROIs[roiIdx].id()
            roiKey  = roiIdx if rid>0 else abs(rid)-1
            en      = h.en()
            layer   = h.layerId()
            isNoise = True if rid<0 else False
            regIdx  = h.signalRegion()
            roiList[roiKey].AddHit(en=en, layer=layer, isNoise=isNoise, regIdx=regIdx)

        for r in roiList:
            varvals = []
            genP4 = roiList[r].genP4
            varvals += [genP4.E(),genP4.Eta(),genP4.Phi()]
            
            for ireg in range(1,NREG+1):
                recP4 = roiList[r].RecoP4(ireg)
                noiseROI = roiList[r].NoiseInROI(ireg)
                varvals += [recP4.E(),noiseROI]
                for il in range(1,NLAYERS+1):
                    recEn = roiList[r].RecoEnergyDeposited(ireg, il)
                    noiseLayer = roiList[r].NoiseInLayer(ireg, il)
                    varvals += [recEn, noiseLayer]

                    if ( FLAGS.samples == "inner" and genP4.Eta() < FLAGS.maxgeneta or
                         FLAGS.samples == "outer" and genP4.Eta() > FLAGS.mingeneta ):
                        fracEn[ireg-1,il-1] += roiList[r].FractionEnergyDeposited(ireg, il)
                        countfracEn[ireg-1,il-1] += 1

            output_tuple.Fill(array.array("f", varvals))
    fOut.cd()
    
    #Write energy distribution of complete showers for complete/incomplete discrimination
    fracEn /= countfracEn
    bins = Carray('d', np.arange(0.5,NLAYERS+1,1))
    for ireg in range(NREG):
        h = TH1F('complete_showers_standard_sr'+str(ireg+1), 
                 'complete_showers_standard_sr'+str(ireg+1),
                 len(bins)-1, bins)
        for il in range(NLAYERS):
            b = h.FindBin(il+1)
            h.Fill(b,fracEn[ireg,il])
        h.Write()

    #output_tuple.Write()
    fOut.Write()
    fOut.Close()

if __name__ == "__main__":
    studies = PartialWafersStudies()
    NREG, NLAYERS = studies.nsr, studies.nlayers
    main()
