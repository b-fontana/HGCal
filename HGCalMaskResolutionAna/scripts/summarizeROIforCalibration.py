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
subdetectors = ["CEE", "HEF", "HEB"]

def main():
    """
    Creates a simple summary of the ROI for fast calibration.
    """
    fIn = TFile.Open(FLAGS.noPUFile, 'READ')
    fOut = TFile(FLAGS.outpath, 'RECREATE')
    output_tuples = []

    for idet,subdet in enumerate(subdetectors):
        if(idet==0):
            NLAYERS = 28
        else:
            NLAYERS = 24
        varnames  = ['genen', 'geneta', 'genphi']
        for i in range(1,NREG+1):
            varnames += ['en_sr{}_ROI'.format(i), 
                         'noise_sr{}_ROI'.format(i),
                         'Cglob1_sr{}_ROI'.format(i),
                         'Cglob2_sr{}_ROI'.format(i),
                         'Cglob3_sr{}_ROI'.format(i)]
            for il in range(1,NLAYERS+1):
                varnames += ['en_sr{}_layer{}'.format(i,il), 
                             'noise_sr{}_layer{}'.format(i,il)]
        output_tuples.append(TNtuple(subdet,subdet,":".join(varnames)))
        
        tree = fIn.Get('an_mask/'+subdet)
        #for i in range(t.GetEntriesFast()):
        for t in tree:
            #define the ROIs
            roiList={}
            for ir in range(0,t.ROIs.size()):
                if t.ROIs[ir].id() < 0 and t.ROIs[ir].id() != -211: 
                    continue
                roiList[ir] = ROISummary(t.ROIs[ir].p4(), Nlayers=NLAYERS)

            for h in t.Hits:
                roiIdx  = h.associatedROI()
                rid     = t.ROIs[roiIdx].id()
                roiKey  = roiIdx if rid>0 or rid==-211 else abs(rid)-1
                mipflag = True
                en      = h.en(mipflag)
                layer   = int(h.layerId()) #originally it is long
                isNoise = True if rid<0 and rid!=-211 else False
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

                    Cglob1, Cglob2, Cglob3 = roiList[r].Cglob(ireg) 
                    varvals += [Cglob1,Cglob2,Cglob3]

                    for il in range(1,NLAYERS+1):
                        recEn = roiList[r].RecoEnergyDeposited(ireg, il)
                        noiseLayer = roiList[r].NoiseInLayer(ireg, il)
                        varvals += [recEn, noiseLayer]

                output_tuples[idet].Fill(array.array("f", varvals))
    fOut.cd()
    fOut.Write()
    fOut.Close()

if __name__ == "__main__":
    studies = PartialWafersStudies()
    NREG = studies.nsr
    main()
