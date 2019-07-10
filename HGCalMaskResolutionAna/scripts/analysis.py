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
from UserCode.HGCalMaskVisualProd.SystemUtils import EtaStr as ES
from UserCode.HGCalMaskResolutionAna.SoftwareCorrection import IncompleteShowersCorrection
from UserCode.HGCalMaskResolutionAna.Calibration import Calibration
from UserCode.HGCalMaskResolutionAna.PartialWafersStudies import PartialWafersStudies
from UserCode.HGCalMaskVisualProd.RootPlotting import RootPlotting
from UserCode.HGCalMaskVisualProd.RootObjects import RootHistograms

def plotHistograms(histos, cdims, pcoords, cname):
    """
    Plots histograms from a list.
    Arguments:
    -> histos: list of histograms
    -> cdims & pccords: as described in RooUtils
    -> cname: name of the canvas to be created
    """
    if not isinstance(histos, list):
        raise TypeError('The histograms have to be passed in a list.,')

    npads = 6 if FLAGS.mode == 2 else 12
    hdiv = []
    legends1 = [TLegend(0.56, 0.66, 0.86, 0.89) for  _ in range(3)]
    legends2 = [TLegend(0.68, 0.75, 0.86, 0.89) for _ in range(3)]
    with RootPlotting(ncanvas=1, npads=npads, cdims=cdims, pcoords=pcoords) as plot:
        if FLAGS.mode == 1:
            titles = ['Resolution', 'Resolution', 'Nevents', 'Resolution / Nevents']
            it = -1
            for ih in range(len(histos)):
                if ih%3==0: it += 1
                h = histos[ih]
                plot.plotHistogram(cpos=0, ppos=ih, h=h, 
                                   title=titles[it], draw_options='colz')
                tex = plot.setLatex(ts=0.04)
                if FLAGS.samples == 'inner':
                    pass
                    #tex.DrawLatex(0.75,0.93,' |#eta| < '+str(FLAGS.maxgeneta))
                elif FLAGS.samples == 'outer':
                    pass
                    #tex.DrawLatex(0.75,0.93,' |#eta| > '+str(FLAGS.mingeneta))

        elif FLAGS.mode == 2:
            for ih in range(NREG):
                if FLAGS.apply_weights:
                    for ixx in range(2):
                        idx = ih+NREG*ixx
                        plot.plotHistogram(cpos=0, ppos=idx, h=histos[idx], 
                                           lw=3,mc=4,msize=.5,lc=4, 
                                           draw_options='E')
                        tex = plot.setLatex()
                        plot.fitHistogram(h=histos[idx], fname='crystalball', 
                                          frange=(-1.,1.), tex=tex)
                        """
                        if False:
                            if FLAGS.samples == 'inner':
                                tex.DrawLatex(0.58,0.92,'Inner radius;  SR{}'
                                              .format((ih%3)+1))
                                tex.DrawLatex(0.7,0.84,' |#eta| < '+str(FLAGS.maxgeneta))
                            elif FLAGS.samples == 'outer':
                                tex.DrawLatex(0.58,0.92,'Outer radius;  SR{}'
                                              .format((ih%3)+1))
                                tex.DrawLatex(0.7,0.84,' |#eta| > '+str(FLAGS.mingeneta))
                        """
                        if FLAGS.samples == 'inner':
                            tex.DrawLatex(0.58,0.92,'Inner radius;  SR{}'
                                          .format((ih%3)+1))
                            tex.DrawLatex(0.62,0.84,str(FLAGS.maxgeneta)+' < |#eta| < '+str(FLAGS.etacuts[-2]))
                        elif FLAGS.samples == 'outer':
                            tex.DrawLatex(0.58,0.92,'Outer radius;  SR{}'
                                          .format((ih%3)+1))
                            tex.DrawLatex(0.62,0.84,str(FLAGS.etacuts[1])+' > |#eta| > '+str(FLAGS.mingeneta))
                        tex.DrawLatex(0.11,0.92,'#bf{CMS} #it{simulation preliminary}')
                        tex.SetTextAlign(31)

                elif not FLAGS.apply_weights:
                    linec = [4, 2, 3, 7]
                    plot.plotHistogram(cpos=0, ppos=ih, h=histos[ih], 
                                       lw=3,mc=linec[0],msize=.5,lc=linec[0], 
                                       draw_options='E')
                    plot.plotHistogram(cpos=0, ppos=ih, h=histos[ih+3], 
                                       lw=3, mc=linec[1], msize=.5, lc=linec[1],
                                       draw_options='same E')
                    plot.plotHistogram(cpos=0, ppos=ih, h=histos[ih+6], 
                                       lw=3, mc=linec[2], msize=.5, lc=linec[2],
                                       draw_options='same E')
                    plot.plotHistogram(cpos=0, ppos=ih, h=histos[ih+9], 
                                       lw=3, mc=linec[3], msize=.5, lc=linec[3],
                                       draw_options='same E')
                    tex = plot.setLatex()
                    if FLAGS.samples == 'inner':
                        tex.DrawLatex(0.58,0.92,'Inner radius;  SR{}'.format((ih%3)+1))
                        legends1[ih].AddEntry(histos[ih], 
                                              ' |#eta| < '+str(FLAGS.etacuts[0]), 'L')
                        for ic in range(len(FLAGS.etacuts)-1):
                            legends1[ih].AddEntry(histos[ih+3*(ic+1)], 
                        str(FLAGS.etacuts[ic])+' #leq |#eta| < '+str(FLAGS.etacuts[ic+1]), 
                                                  'L')
                    elif FLAGS.samples == 'outer':
                        tex.DrawLatex(0.58,0.92,'Outer radius;  SR{}'.format((ih%3)+1))
                        legends1[ih].AddEntry(histos[ih], 
                                              ' |#eta| #geq '+str(FLAGS.etacuts[-1]), 
                                              'L')
                        for ic in range(len(FLAGS.etacuts)-1):
                            legends1[ih].AddEntry(histos[ih+3*(ic+1)], 
                str(FLAGS.etacuts[-(ic+2)])+' #leq |#eta| < '+str(FLAGS.etacuts[-(ic+1)]), 
                                                  'L')
                    tex.DrawLatex(0.11,0.92,'#bf{CMS} #it{simulation preliminary}')
                    tex.SetTextAlign(31)
                    legends1[ih].Draw()

                    hdiv.append(histos[ih+3].Clone('weight1_sr{}'.format(ih+1)))
                    hdiv.append(histos[ih+6].Clone('weight2_sr{}'.format(ih+1)))
                    hdiv.append(histos[ih+9].Clone('weight3_sr{}'.format(ih+1)))
                    for idiv in range(len(FLAGS.etacuts)-1):
                        hdiv[-3+idiv].Divide(histos[ih])
                        extrastr = '' if idiv==0 else 'same'
                        hdiv[-3+idiv].GetYaxis().SetRangeUser(0., 2.)
                        plot.plotHistogram(cpos=0, ppos=ih+3, h=hdiv[-3+idiv],
                                           yaxis_title='Weight',
                                           lw=3, mc=linec[idiv+1], msize=.5, 
                                           lc=linec[idiv+1],
                                           draw_options='HIST'+extrastr, copy=True)
                    tex = plot.setLatex()
                    if FLAGS.samples == 'inner':
                        tex.DrawLatex(0.58,0.92,'Inner radius;  SR{}'.format((ih%3)+1))
                    elif FLAGS.samples == 'outer':
                        tex.DrawLatex(0.58,0.92,'Outer radius;  SR{}'.format((ih%3)+1))
                    for iv in range(len(FLAGS.etacuts)-1):
                        legends2[ih].AddEntry(hdiv[iv], 'weight'+str(iv+1), 'L')
                        legends2[ih].Draw()
                    tex.DrawLatex(0.11,0.92,'#bf{CMS} #it{simulation preliminary}')
                    tex.SetTextAlign(31)

        plot.save(cpos=0, name=cname)

    if not FLAGS.apply_weights:
        save_str = ( 'calibshowers_mask'+str(FLAGS.mask)+'_'+
                     FLAGS.samples+'_mode'+str(FLAGS.mode) )
        RootHistograms(histos).save(save_str)
        RootHistograms(hdiv).save(save_str, mode='UPDATE')

def main():
    #gStyle.SetOptStat(0)
    gROOT.SetBatch(True)
    gStyle.SetPalette(kTemperatureMap)

    fIn=TFile.Open(FLAGS.noPUFile)
    data=fIn.Get('data')

    calibration = Calibration(FLAGS.mingenen, etaregions,
                              FLAGS.plotLabel, FLAGS.samples, FLAGS.mask, FLAGS.outpath)
    calibration.L0L1Calibration(FLAGS.noPUFile)

    with open('calib_'+FLAGS.samples+"_"+str(FLAGS.mask)+'_nopu.pck','w') as cachefile:
        pickle.dump(calibration.calib, cachefile, pickle.HIGHEST_PROTOCOL)
    with open('calib_'+FLAGS.samples+"_"+str(FLAGS.mask)+'_nopu.pck','r') as cachefile:
        calib = pickle.load(cachefile)

    if FLAGS.apply_weights:
        calibshowers_str = ( 'calibshowers_mask'+str(FLAGS.mask)+'_'+
                             FLAGS.samples+'_mode2')
        showercorr = IncompleteShowersCorrection(calibshowers_str+'.root',
                                                 Av(FLAGS.etacuts))
        weights = showercorr.getCorrectionWeights()
        boundaries = [5, 5, 5]
        corr_mode = 'right' if FLAGS.samples == 'outer' else 'left'
        lowstats_factors = showercorr.calculateLowStatisticsFactor(boundaries, corr_mode)
        weights_graphs = [showercorr.buildCorrectionWeightsGraphs(region=i+1)
                          for i in range(NREG)]

    histos=OrderedDict()
    etabins = 12
    etacuts = FLAGS.etacuts

    if FLAGS.mode == 1:
        if FLAGS.samples == 'inner':
            phibins, etainf, etasup = 12, 2.69, 3.04
        elif FLAGS.samples == 'outer':
            phibins, etainf, etasup = 12, 1.44, 1.66
        hn = ['den{}', 'den{}_2D_res', 'den{}_2D_events']
        for ireg in range(1,NREG+1):
            histos[hn[0].format(ireg)] = TH1F(hn[0].format(ireg),';#Delta E/E;PDF',
                                                  100, -1.05, .8)
            histos[hn[1].format(ireg)] = TH2F(hn[1].format(ireg), ';|#eta|;#phi',
                                                  50, etainf, etasup,
                                                  phibins, -TMath.Pi(), TMath.Pi())
            histos[hn[2].format(ireg)] = TH2F(hn[2].format(ireg), ';|#eta|;#phi',
                                                  50, etainf, etasup,
                                                  phibins, -TMath.Pi(), TMath.Pi())
    elif FLAGS.mode == 2:
        hn = ['res_complete_before{}',          'res_complete_after{}', 
              'res_incomplete_before_left{}',   'res_incomplete_after_left{}',
              'res_incomplete_before_right{}',  'res_incomplete_after_right{}',
              'en{}_per_layer_signal',          'en{}_per_layer_bckg1',
              'en{}_per_layer_bckg2',           'en{}_per_layer_bckg3',
              'noise{}_per_layer_signal',       'noise{}_per_layer_bckg1',
              'noise{}_per_layer_bckg2',        'noise{}_per_layer_bckg3']
        for ireg in range(1,NREG+1):
            bins = Carray('d', np.arange(-1.05, .8, 0.01))
            strings = ';#Delta E/E_{gen};PDF'
            for ih in range(6):
                histos[hn[ih].format(ireg)] = TH1F(hn[ih].format(ireg), strings,
                                               len(bins)-1, bins)
            bins = Carray('d', np.arange(0.5,29,1.))
            strings = ';Layer;E_{reco} / E_{gen}'
            for ih in range(6,14):
                histos[hn[ih].format(ireg)] = TProfile(hn[ih].format(ireg), strings,
                                                       len(bins)-1, bins)

    else:
        raise ValueError('The number of components has to be 1 or 2.')

    for h in histos:
        histos[h].Sumw2()
        histos[h].SetMarkerStyle(20)
        histos[h].SetDirectory(0)

    counts = np.zeros((3,28))

    for i in range(0,data.GetEntriesFast()):
        data.GetEntry(i)
        genen    = getattr(data,'genen')        
        geneta   = abs(getattr(data,'geneta'))
        genphi   = getattr(data,'genphi')

        sc = (FLAGS.maxgeneta, FLAGS.mingeneta)
        singlecut = ( geneta < sc[0] if FLAGS.samples == 'inner' 
                      else geneta > sc[1] )
        insidecut = ( geneta < FLAGS.etacuts[-2] if FLAGS.samples == 'inner' 
                      else geneta > FLAGS.etacuts[1] )

        for ireg in range(1,NREG+1):
            recen = getattr(data,'en_sr{}_ROI'.format(ireg))
            avgnoise = getattr(data,'noise_sr3_ROI')*A[ireg-1]/A[2]

            #Calibration factors. f2 is used for PU.
            f1, f2 = 1., 0.
            etaregions_shifted = np.roll(etaregions, shift=-1)[:-1]
            for ieta1,ieta2 in zip(etaregions[:-1], etaregions_shifted):
                if (geneta < ieta1 or geneta > ieta2): 
                    continue
                idstr = 'sr{}_from{}to{}'.format(ireg, ES(ieta1), ES(ieta2))
                if 'L0' in calib:
                    f1 /= calib['L0'][idstr].Eval(geneta)+1.0
                    if 'L1' in calib:
                        f1 /= calib['L1'][idstr].Eval(f1*recen)+1.0    
                        if 'L2' in calib and ireg in calib['L2']:
                            f2 = calib['L2'][idstr].Eval(avgnoise)

            recen = f1*recen - f2 
            recen_shift = recen / .49
            deltaE = recen/genen-1.
            deltaE_shift = recen_shift/genen-1.

            ###Store the energy resolution###
            if FLAGS.mode == 1:
                histos[hn[0].format(ireg)].Fill(deltaE)
                histos[hn[1].format(ireg)].Fill(geneta, genphi, deltaE)
                histos[hn[2].format(ireg)].Fill(geneta, genphi)

            elif FLAGS.mode == 2:
                #differentiate complete from incomplete showers
                bool_sig = ( (FLAGS.samples == 'inner' and geneta < FLAGS.maxgeneta) or 
                             (FLAGS.samples == 'outer' and geneta >= FLAGS.mingeneta) )
                bools_bckg = []
                for ic in range(len(etacuts)-1):
                    bools_bckg.append((FLAGS.samples=='inner' 
                                       and geneta >= etacuts[ic] 
                                       and geneta < etacuts[ic+1])
                                      or (FLAGS.samples == 'outer' 
                                          and geneta < etacuts[-(ic+1)] 
                                          and geneta >= etacuts[-(ic+2)]) )
                check_bools = int(bool_sig)
                for ic in range(len(etacuts)-1):
                    check_bools += int(bools_bckg[ic])
                assert check_bools == 1

                ###Calculate and calibrate the energy per layer###
                recen_corr = 0 
                for il in range(1,NLAYERS+1):
                    #if: 'signal-like': complete showers
                    #elif: 'background-like': incomplete showers
                    if FLAGS.apply_weights:
                        if FLAGS.samples == 'inner':
                            weight_limit = il > boundaries[ireg-1] 
                        else:
                            weight_limit = il < boundaries[ireg-1] 
                    
                    b = histos[hn[6].format(ireg)].FindBin(il)
                    if bool_sig:
                        v = f1*getattr(data,'en_sr{}_layer{}'.format(ireg,il)) - f2
                        histos[hn[6].format(ireg)].Fill(b,v/recen)
                        if FLAGS.apply_weights:
                            recen_corr += v
                        v = (f1*getattr(data,'noise_sr3_layer{}'.format(il))
                             *A[ireg-1]/A[2] - f2)
                        histos[hn[10].format(ireg)].Fill(b,v/recen)
                    else:
                        lshift = [.9, .85, .7] #layer shift
                        for w in range(len(etacuts)-1):
                            if bools_bckg[w]:      
                                v = f1*getattr(data,'en_sr{}_layer{}'.format(ireg,il)) - f2
                                if w==1: counts[ireg-1,il-1] += 1
                                histos[hn[7+w].format(ireg)].Fill(b*lshift[w],v/recen)
                                if ( FLAGS.apply_weights and 
                                     weights[ireg-1][w][il-1]!=0 and
                                     weight_limit):
                                    recen_corr += v/weights[ireg-1][w][int(round((il-1)*lshift[w],0))]
                                    #weight_graphs[ireg][il].SetBit(weight_graphs[ireg][il].klsSortedX)
                                    #weight_graphs[ireg][il].Eval(geneta, spline=0, 'S')
                                v = (f1*getattr(data,'noise_sr3_layer{}'.format(il))
                                     *A[ireg-1]/A[2] - f2)
                                histos[hn[11+w].format(ireg)].Fill(b,v/recen)

                if FLAGS.apply_weights:
                    if not singlecut and insidecut:
                        recen_corr *= (1 / (1-lowstats_factors[ireg-1]) )
                        recen_corr *= 1/(.75)
                        deltaE_corr = recen_corr/genen-1.
                        if deltaE<-0.3:
                            histos[hn[2].format(ireg)].Fill(deltaE_shift)
                            histos[hn[3].format(ireg)].Fill(deltaE_corr)
                        else:
                            histos[hn[4].format(ireg)].Fill(deltaE)
                            histos[hn[5].format(ireg)].Fill(deltaE_corr)
                    elif singlecut: #complete shower
                        deltaE_corr = recen_corr/genen-1.
                        histos[hn[0].format(ireg)].Fill(deltaE)
                        histos[hn[1].format(ireg)].Fill(deltaE_corr)
                
    #end of tree loop    
    fIn.Close()

    if FLAGS.mode == 1:
        pcoords = [[[0.01,0.76,0.33,0.99],   #canvas0, pad0
                    [0.34,0.76,0.66,0.99],   #canvas0, pad1
                    [0.67,0.76,0.99,0.99],   #canvas0, pad2
                    [0.01,0.51,0.33,0.74],   #canvas0, pad3
                    [0.34,0.51,0.66,0.74],   #canvas0, pad4
                    [0.67,0.51,0.99,0.74],   #canvas0, pad5
                    [0.01,0.26,0.33,0.49],   #canvas0, pad6
                    [0.34,0.26,0.66,0.49],   #canvas0, pad7
                    [0.67,0.26,0.99,0.49],   #canvas0, pad8
                    [0.01,0.01,0.33,0.24],   #canvas0, pad9
                    [0.34,0.01,0.66,0.24],   #canvas0, pad10
                    [0.67,0.01,0.99,0.24]]]  #canvas0, pad11
        cdims = [[2000,2000]]
        picname = '1comp_'+FLAGS.samples
        if FLAGS.apply_weights:
            picname += '_corrected' 
    else:
        #pcoords = [[[0.01,0.755,0.33,0.99],
        #            [0.34,0.755,0.66,0.99],  
        #            [0.67,0.755,0.99,0.99],
        #            [0.01,0.505,0.33,0.745],   
        #            [0.34,0.505,0.66,0.745],   
        #            [0.67,0.505,0.99,0.745],   
        #            [0.01,0.255,0.33,0.495],   
        #            [0.34,0.255,0.66,0.495],   
        #            [0.67,0.255,0.99,0.495],
        #            [0.01,0.005,0.33,0.245],   
        #            [0.34,0.005,0.66,0.245],   
        #            [0.67,0.005,0.99,0.245]]]
        pcoords = [[[0.01,0.51,0.33,0.99],
                    [0.34,0.51,0.66,0.99],  
                    [0.67,0.51,0.99,0.99],
                    [0.01,0.01,0.33,0.49],   
                    [0.34,0.01,0.66,0.49],   
                    [0.67,0.01,0.99,0.49]]]
        cdims = [[1000,600]]
        picname = '2comp_'+FLAGS.samples
        if FLAGS.apply_weights:
            picname += '_corrected' 
    
    correct_order = []
    for i in range(len(hn)):
        for ireg in range(1,NREG+1):
            correct_order.append(hn[i].format(ireg))
    assert len(correct_order) == len(histos.keys())
    histos = [histos[correct_order[i]] for i in range(len(correct_order))]
    if FLAGS.mode == 1:
        histos.append(histos[3].Clone())
        histos.append(histos[4].Clone())
        histos.append(histos[5].Clone())
        histos[-3].Divide(histos[6])
        histos[-2].Divide(histos[7])
        histos[-1].Divide(histos[8])
        plotHistograms(histos, cdims, pcoords, 
                       os.path.join(FLAGS.outpath,picname))
    elif FLAGS.mode == 2:
        if FLAGS.apply_weights:
            histos_left = histos[6:12]
            histos_right = histos[12:18]
            plotHistograms(histos_left, cdims, pcoords, 
                           os.path.join(FLAGS.outpath,picname+'_left'))
            plotHistograms(histos_right, cdims, pcoords, 
                           os.path.join(FLAGS.outpath,picname+'_right'))
        else:
            histos = histos[18:42]
            plotHistograms(histos, cdims, pcoords, 
                           os.path.join(FLAGS.outpath,picname))

if __name__ == "__main__":
    parser = Argparser.Argparser()
    FLAGS = parser.get_flags()
    parser.print_args()
    base = PartialWafersStudies()
    NREG, NLAYERS, A = base.nsr, base.nlayers, base.sr_area
    etaregions = np.round(np.arange(2.7,3.031,0.001).tolist(), 3).tolist()
    main()
