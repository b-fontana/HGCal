import sys
import ROOT
import array
from UserCode.HGCalMaskResolutionAna import Argparser

from ROOT import TFile, TNtuple
from UserCode.HGCalMaskResolutionAna.ROISummary import ROISummary

parser = Argparser.Argparser()
FLAGS = parser.get_flags()
parser.print_args()

def main():
    """
    Creates a simple summary of the ROI for fast calibration.
    """
    NREG, NLAYERS = 3, 28
    fOut = TFile(FLAGS.outpath+'.root', 'RECREATE')
    fOut.cd()
    varnames  = ['genen','geneta','genphi']
    for i in range(1,NREG+1):
        varnames += ['en_sr{}_ROI'.format(i), 
                     'noise_sr{}_ROI'.format(i)]
        for il in range(1,NLAYERS+1):        
            varnames += ['en_sr{}_layer{}'.format(i,il), 
                         'noise_sr{}_layer{}'.format(i,il)]
    output_tuple = TNtuple("data","data",":".join(varnames))

    fIn = TFile.Open(FLAGS.noPUFile+'.root')
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
            roiList[roiKey].addHit(en=en, layer=layer, isNoise=isNoise, regIdx=regIdx)

        """
        CHECK
        import numpy as np
        all1 = [roiList[0].getRecoEnergyDeposited(1,x) for x in range(1,29)]
        if not (np.isclose(roiList[0].getRecoP4(1).E(),sum(all1))):
            print(roiList[0].getRecoP4(1).E(), sum(all1))
        all2 = [roiList[0].getRecoEnergyDeposited(2,x) for x in range(1,29)]
        if not (np.isclose(roiList[0].getRecoP4(2).E(),sum(all2))):
            print(roiList[0].getRecoP4(2).E(), sum(all2))
        all3 = [roiList[0].getRecoEnergyDeposited(3,x) for x in range(1,29)]
        if not (np.isclose(roiList[0].getRecoP4(3).E(), sum(all3))):
            print(roiList[0].getRecoP4(3).E(), sum(all3))
        continue
        """
        for r in roiList:
            varvals = []
            genP4 = roiList[r].getGenP4()
            varvals += [genP4.E(),genP4.Eta(),genP4.Phi()]
            
            for ireg in range(1,NREG+1):
                recP4 = roiList[r].getRecoP4(ireg)
                noiseROI = roiList[r].getNoiseInROI(ireg)
                varvals += [recP4.E(),noiseROI]
                for il in range(1,NLAYERS+1):
                    recEn = roiList[r].getRecoEnergyDeposited(ireg, il)
                    noiseLayer = roiList[r].getNoiseInLayer(ireg, il)
                    varvals += [recEn, noiseLayer]
            
            output_tuple.Fill(array.array("f", varvals))
        
    fOut.cd()
    output_tuple.Write()
    fOut.Close()

if __name__ == "__main__":
    main()
