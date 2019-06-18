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
    fOut = TFile(FLAGS.outpath, "RECREATE")
    fOut.cd()
    varnames  = ['genen','geneta','genphi']
    for i in range(1,NREG+1):
        for il in range(1,NLAYERS+1):
            varnames += ['en_sr{}_layer{}'.format(i,il), 
                         'noise_sr{}_layer{}'.format(i,il)]
    output_tuple = TNtuple("data","data",":".join(varnames))

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
            roiList[roiKey].addHit(en=en, layer=layer, isNoise=isNoise, regIdx=regIdx)

        #mc truth
        for r in roiList:
            varvals = []
            genP4 = roiList[r].getGenP4()
            varvals += [genP4.E(),genP4.Eta(),genP4.Phi()]
            
            for ireg in range(1,NREG+1):
                for il in range(1,NLAYERS+1):
                    recP4    = roiList[r].getReconstructedP4(ireg, layer)
                    noise    = roiList[r].getNoiseInRegion(ireg, layer)
                    varvals += [recP4.E(),noise]
            output_tuple.Fill(array.array("f", varvals))

    fOut.cd()
    output_tuple.Write()
    fOut.Close()

if __name__ == "__main__":
    main()
