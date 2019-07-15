import os, sys
from UserCode.HGCalMaskResolutionAna import Argparser

from ROOT import TFile, gROOT, TLegend
from UserCode.HGCalMaskResolutionAna.PartialWafersStudies import PartialWafersStudies
from UserCode.HGCalMaskVisualProd.RootPlotting import RootPlotting
from UserCode.HGCalMaskVisualProd.RootPlotting import standardPadCoords as Coords


def main():
    fIn = [TFile.Open('allplots_'+FLAGS.samples+'_'+str(i)+'.root', 'READ') for i in range(3,7)]
    hn = ['bias_vs_eta_after{}', 'indep_vs_eta_after{}']
    hrms = [fIn[i].Get(hn[0].format(ireg)) for ireg in range(1,NREG+1) for i in range(len(fIn))]
    hindep = [fIn[i].Get(hn[1].format(ireg)) for ireg in range(1,NREG+1) for i in range(len(fIn))]

    with RootPlotting(ncanvas=1, npads=6, cdims=[[1000,600]], pcoords=Coords(len(hn))) as plot:
        lc = [4, 2, 3, 7]
        legends1 = [TLegend(0.12, 0.12, 0.44, 0.25) for  _ in range(3)]
        legends2 = [TLegend(0.12, 0.76, 0.44, 0.89) for  _ in range(3)]
        for ireg in range(NREG):
            for ih in range(len(fIn)):
                opt = 'colz' if ih==0 else 'colz same'
                plot.plotHistogram(cpos=0, ppos=ireg, lc=lc[ih], mc=lc[ih], title='Bias vs eta', 
                                   h=hrms[ih*NREG+ireg], yranges=(-2.,0.3), draw_options=opt)
                legends1[ireg].AddEntry(hrms[ih*NREG+ireg], 'Mask '+str(ih+3), 'L')
                legends1[ireg].Draw()
                plot.plotHistogram(cpos=0, ppos=ireg+3, lc=lc[ih], mc=lc[ih], title='RMS/(1+Bias) vs eta', 
                                   h=hindep[ih*NREG+ireg], yranges=(-0.5, 7.), draw_options=opt)
                legends2[ireg].AddEntry(hindep[ih*NREG+ireg], 'Mask '+str(ih+3), 'L')
                legends2[ireg].Draw()
        plot.save(cpos=0, name='RMS_Bias.png')

if __name__ == "__main__":
    gROOT.SetBatch(True)
    parser = Argparser.Argparser()
    FLAGS = parser.get_flags()
    parser.print_args()
    base = PartialWafersStudies()
    NREG, NLAYERS, A = base.nsr, base.nlayers, base.sr_area
    main()
