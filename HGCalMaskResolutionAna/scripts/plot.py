import os, sys
from UserCode.HGCalMaskResolutionAna import Argparser

from ROOT import TFile, gROOT, TLegend
from UserCode.HGCalMaskResolutionAna.PartialWafersStudies import PartialWafersStudies
from UserCode.HGCalMaskVisualProd.RootPlotting import RootPlotting
from UserCode.HGCalMaskVisualProd.RootPlotting import standardPadCoords as Coords

def RMS_Bias(f):
    hn = ['biasVSeta{}', 'indepVSeta{}']
    hrms = [f[i].Get(hn[0].format(ireg)) 
            for i in range(len(f)) for ireg in range(1,NREG+1)]
    hindep = [f[i].Get(hn[1].format(ireg)) 
              for i in range(len(f)) for ireg in range(1,NREG+1)]
    with RootPlotting(ncanvas=1, npads=len(hn)*NREG, 
                      cdims=[[2000,1200]], pcoords=Coords(len(hn))) as plot:
        lc = [4, 2, 3, 7]
        legends1 = [TLegend(0.61, 0.12, 0.84, 0.32) for  _ in range(3)]
        legends2 = [TLegend(0.61, 0.69, 0.84, 0.89) for  _ in range(3)]
        for ireg in range(NREG):
            for ih in range(len(f)):
                opt = 'AP' if ih==0 else 'P'
                plot.plotGraph(cpos=0, ppos=ireg,
                               lc=lc[ih], mc=lc[ih], title='Bias vs eta', 
                               g=hrms[ih*NREG+ireg], yranges=(-1.,1.), 
                               xaxis_title='|#eta|', yaxis_title='Bias',
                               draw_options=opt)
                legends1[ireg].AddEntry(hrms[ih*NREG+ireg], 'Mask '+str(ih+3), 'L')
                legends1[ireg].Draw()
                plot.plotGraph(cpos=0, ppos=ireg+3, 
                               lc=lc[ih], mc=lc[ih], title='RMS/(1+Bias) vs eta', 
                               g=hindep[ih*NREG+ireg], yranges=(0., 1.),
                               xaxis_title='|#eta|', yaxis_title='RMS/(1+Bias)', 
                               draw_options=opt)
                legends2[ireg].AddEntry(hindep[ih*NREG+ireg], 'Mask '+str(ih+3), 'L')
                legends2[ireg].Draw()
        plot.save(cpos=0, name='RMS_Bias_'+FLAGS.samples+'_'+FLAGS.method+'.png')

def Resolution(f):
    hn = ['res_complete_before{}', 'res_complete_after{}',
          'res_incomplete_before{}', 'res_incomplete_after{}',
          'res_vs_eta_before{}', 'res_vs_eta_after{}']
    h1 =   [f[i].Get(hn[0].format(ireg)) 
            for i in range(len(f)) for ireg in range(1,NREG+1)]
    h2 =   [f[i].Get(hn[1].format(ireg)) 
            for i in range(len(f)) for ireg in range(1,NREG+1)]
    h3 =   [f[i].Get(hn[2].format(ireg)) 
            for i in range(len(f)) for ireg in range(1,NREG+1)]
    h4 =   [f[i].Get(hn[3].format(ireg)) 
            for i in range(len(f)) for ireg in range(1,NREG+1)]
    h5 =   [f[i].Get(hn[4].format(ireg)) 
            for i in range(len(f)) for ireg in range(1,NREG+1)]
    h6 =   [f[i].Get(hn[5].format(ireg)) 
            for i in range(len(f)) for ireg in range(1,NREG+1)]
    
    legends1 = [TLegend(0.12, 0.76, 0.44, 0.89) for  _ in range(3)]
    for imask in range(len(f)):
        with RootPlotting(ncanvas=1, npads=len(hn)*NREG,
                          cdims=[[1600,2000]], pcoords=Coords(len(hn))) as plot:

            for ireg in range(NREG):
                idx = ireg+imask*NREG
                plot.plotHistogram(cpos=0, ppos=ireg, 
                                   h=h1[idx],
                                   lw=3,mc=4,msize=.5,lc=4, 
                                   draw_options='colz')
                plot.plotHistogram(cpos=0, ppos=ireg+NREG, 
                                   h=h2[idx], 
                                   lw=3,mc=4,msize=.5,lc=4, 
                                   draw_options='colz')
                plot.plotHistogram(cpos=0, ppos=ireg+2*NREG, 
                                   h=h3[idx], 
                                   lw=3,mc=4,msize=.5,lc=4, 
                                   draw_options='colz')
                plot.plotHistogram(cpos=0, ppos=ireg+3*NREG, 
                                   h=h4[idx], 
                                   lw=3,mc=4,msize=.5,lc=4, 
                                   draw_options='colz')
                plot.plotHistogram(cpos=0, ppos=ireg+4*NREG, 
                                   h=h5[idx], 
                                   lw=3,mc=4,msize=.5,lc=4, 
                                   draw_options='colz')
                plot.plotHistogram(cpos=0, ppos=ireg+5*NREG, 
                                   h=h6[idx], 
                                   lw=3,mc=4,msize=.5,lc=4, 
                                   draw_options='colz')
            plot.save(cpos=0, name='Resolution'+str(imask+3)+'.png')

def main():
    end = '_'+FLAGS.method+'.root'
    fIn = [TFile.Open('allplots_'+FLAGS.samples+'_'+str(i)+end, 'READ') for i in range(3,7)]
    RMS_Bias(fIn)
    #Resolution(fIn)

if __name__ == "__main__":
    gROOT.SetBatch(True)
    parser = Argparser.Argparser()
    FLAGS = parser.get_flags()
    parser.print_args()
    base = PartialWafersStudies()
    NREG, NLAYERS, A = base.nsr, base.nlayers, base.sr_area
    main()
