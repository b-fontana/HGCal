import os, sys
import numpy as np
import pickle
from UserCode.HGCalMaskResolutionAna import Argparser
from array import array as Carray
from collections import OrderedDict

from ROOT import TCanvas, TLatex, TFile, TMath, TH1F
from ROOT import TLegend, TH2F, TLorentzVector, TProfile, TH1D
from ROOT import gStyle, gROOT, kTemperatureMap
from UserCode.HGCalMaskResolutionAna.SoftwareCorrection import IncompleteShowersCorrection
from UserCode.HGCalMaskVisualProd.RootPlotting import RootPlotting
from UserCode.HGCalMaskVisualProd.RootObjects import RootHistograms

def valAverage(l):
    l = np.array(l)
    l2 = np.roll(l, shift=-1)[:-1]
    l = l[:-1]
    l += (l2 - l) / 2
    return l

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
                    tex.DrawLatex(0.75,0.93,' |#eta| < '+str(FLAGS.maxgeneta))
                elif FLAGS.samples == 'outer':
                    tex.DrawLatex(0.75,0.93,' |#eta| > '+str(FLAGS.mingeneta))

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
                            tex.DrawLatex(0.62,0.84,str(FLAGS.etacuts[0])+' |#eta| > '+str(FLAGS.mingeneta))
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
    gStyle.SetOptStat(0)
    gROOT.SetBatch(True)
    gStyle.SetPalette(kTemperatureMap)

    fIn=TFile.Open(FLAGS.noPUFile)
    data=fIn.Get('data')

    with open('calib_nopu.pck','r') as cachefile:
        calib = pickle.load(cachefile)

    if FLAGS.apply_weights:
        calibshowers_str = ( 'calibshowers_mask'+str(FLAGS.mask)+'_'+
                             FLAGS.samples+'_mode2')
        showercorr = IncompleteShowersCorrection(calibshowers_str+'.root',
                                                 valAverage(FLAGS.etacuts))
        weights = showercorr.getCorrectionWeights()
        boundaries = [8., 8., 8.]
        corr_mode = 'right' if FLAGS.samples == 'outer' else 'left'
        lowstats_factors = showercorr.calculateLowStatisticsFactor(boundaries, corr_mode)
        weights_graphs = [showercorr.buildCorrectionWeightsGraphs(region=i+1)
                          for i in range(NREG)]
    
    histos=OrderedDict()
    etabins = 12
    etacuts = FLAGS.etacuts
    if FLAGS.samples == 'inner':
        phibins, etainf, etasup = 12, 2.69, 3.04
    elif FLAGS.samples == 'outer':
        phibins, etainf, etasup = 30, 1.34, 1.68
    else:
        raise ValueError("Specify a valid value for the '--samples' option.")

    if FLAGS.mode == 1:
        hn = ['den{}', 'den{}_2D_res', 'den{}_2D_events']
        for ireg in range(1,NREG+1):
            histos[hn[0].format(ireg)] = TH1F(hn[0].format(ireg),';#Delta E/E;PDF',
                                                  100, -1.05, .8)
            histos[hn[1].format(ireg)] = TH2F(hn[1].format(ireg), ';|#eta|;#phi',
                                                  50, etainf, etasup,
                                                  12, -TMath.Pi(), TMath.Pi())
            histos[hn[2].format(ireg)] = TH2F(hn[2].format(ireg), ';|#eta|;#phi',
                                                  50, etainf, etasup,
                                                  12, -TMath.Pi(), TMath.Pi())
    elif FLAGS.mode == 2:
        hn = ['res_complete_before{}',          'res_complete_after{}', 
              'res_incomplete_before_left{}',   'res_incomplete_after_left{}',
              'res_incomplete_before_right{}',  'res_incomplete_after_right{}',
              'en{}_per_layer_signal',          'en{}_per_layer_bckg1',
              'en{}_per_layer_bckg2',           'en{}_per_layer_bckg3',
              'noise{}_per_layer_signal',       'noise{}_per_layer_bckg1',
              'noise{}_per_layer_bckg2',        'noise{}_per_layer_bckg3']
        for ireg in range(1,NREG+1):
            histos[hn[0].format(ireg)] = TH1F(hn[0].format(ireg), ';#Delta E/E_{gen};PDF',
                                              100, -1.05, .8)
            histos[hn[1].format(ireg)] = TH1F(hn[1].format(ireg), ';#Delta E/E_{gen};PDF',
                                              100, -1.05, .8)
            histos[hn[2].format(ireg)] = TH1F(hn[2].format(ireg), ';#Delta E/E_{gen};PDF',
                                              200, -2.05, 2.)
            histos[hn[3].format(ireg)] = TH1F(hn[3].format(ireg), ';#Delta E/E_{gen};PDF',
                                              200, -2.05, 2.)
            histos[hn[4].format(ireg)] = TH1F(hn[4].format(ireg), ';#Delta E/E_{gen};PDF',
                                              200, -2.05, 2.)
            histos[hn[5].format(ireg)] = TH1F(hn[5].format(ireg), ';#Delta E/E_{gen};PDF',
                                              200, -2.05, 2.)
            
            bins = Carray('d', np.arange(0.5,29,1.))
            assert len(bins) == NLAYERS+1
            strings = ';Layer;E_{reco} / E_{gen}'
            for ih in range(6,len(hn)):
                histos[hn[ih].format(ireg)] = TProfile(hn[ih].format(ireg), strings,
                                                       NLAYERS, bins)
    else:
        raise ValueError('The number of components has to be 1 or 2.')

    for h in histos:
        histos[h].Sumw2()
        histos[h].SetMarkerStyle(20)
        histos[h].SetDirectory(0)

    #hc1 = TH2F('check1','en_eta_larger', 100, 15, 130, 100, FLAGS.etacuts[-3], FLAGS.etacuts[-2])
    #hc2 = TH2F('check2','eta_phi_larger', 30, FLAGS.etacuts[-3], FLAGS.etacuts[-2], 100, -TMath.Pi(), TMath.Pi())
    #hc3 = TH2F('check3','en_phi_larger', 100, 30, 130, 100, -TMath.Pi(), TMath.Pi())
    #hc4 = TH2F('check4','en_eta_smaller', 100, 30, 130, 30, 2.93, FLAGS.etacuts[-3])
    #hc5 = TH2F('check5','eta_phi_smaller', 30, 2.93, FLAGS.etacuts[-3], 100, TMath.Pi(), TMath.Pi())
    #hc6 = TH2F('check6','en_phi_smaller', 100, 30, 130, 100, -TMath.Pi(), TMath.Pi())

    for i in range(0,data.GetEntriesFast()):
        data.GetEntry(i)
        genen    = getattr(data,'genen')        
        geneta   = abs(getattr(data,'geneta'))
        genphi   = getattr(data,'genphi')

        if FLAGS.mode == 2:
            bool_sig = ( (FLAGS.samples == 'inner' and geneta < etacuts[0]) or 
                         (FLAGS.samples == 'outer' and geneta >= etacuts[-1]) )
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
            deltaE = recen/genen-1.
            ###Store the energy resolution###
            sc = (FLAGS.maxgeneta, FLAGS.mingeneta)
            singlecut = ( geneta < sc[0] if FLAGS.samples == 'inner' 
                          else geneta > sc[1] )
            insidecut = ( geneta < FLAGS.etacuts[-2] if FLAGS.samples == 'inner' 
                          else geneta > FLAGS.etacuts[0] )

            if FLAGS.mode == 1:
                #if singlecut:
                histos[hn[0].format(ireg)].Fill(deltaE)
                histos[hn[1].format(ireg)].Fill(geneta, genphi, deltaE)
                histos[hn[2].format(ireg)].Fill(geneta, genphi)

            elif FLAGS.mode == 2:
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
                        histos[hn[6].format(ireg)].Fill(b,v/genen)
                        if FLAGS.apply_weights:
                            recen_corr += v
                        v = (f1*getattr(data,'noise_sr3_layer{}'.format(il))
                             *A[ireg-1]/A[2] - f2)
                        histos[hn[10].format(ireg)].Fill(b,v/genen)
                    else:
                        for w in range(len(etacuts)-1):
                            if bools_bckg[w]:      
                                v = f1*getattr(data,'en_sr{}_layer{}'.format(ireg,il)) - f2
                                histos[hn[7+w].format(ireg)].Fill(b,v/genen)
                                if ( FLAGS.apply_weights and 
                                     weights[ireg-1][w][il-1]!=0 and
                                     weight_limit):
                                    recen_corr += v/weights[ireg-1][w][il-1]
                                    #weight_graphs[ireg][il].SetBit(weight_graphs[ireg][il].klsSortedX)
                                    #weight_graphs[ireg][il].Eval(geneta, spline=0, 'S')
                                v = (f1*getattr(data,'noise_sr3_layer{}'.format(il))
                                     *A[ireg-1]/A[2] - f2)
                                histos[hn[11+w].format(ireg)].Fill(b,v/genen)

                if FLAGS.apply_weights:
                    if not singlecut and insidecut:
                        recen_corr *= (1 / (1-lowstats_factors[ireg-1]) )
                        deltaE_corr = recen_corr/genen-1.
                        if deltaE<-0.2:
                            histos[hn[2].format(ireg)].Fill(deltaE)
                            histos[hn[3].format(ireg)].Fill(deltaE_corr)
                        else:
                            histos[hn[4].format(ireg)].Fill(deltaE)
                            histos[hn[5].format(ireg)].Fill(deltaE_corr)
                    elif singlecut: #complete shower
                        deltaE_corr = recen_corr/genen-1.
                        histos[hn[0].format(ireg)].Fill(deltaE)
                        histos[hn[1].format(ireg)].Fill(deltaE_corr)
                
                """
                if ( abs(deltaE) < 0.3 and geneta < FLAGS.etacuts[-2] 
                    and FLAGS.etacuts>[-3] ):
                    #if round(deltaE_corr-deltaE,4) > 0.2:
                    hc1.Fill(genen, geneta)
                    hc2.Fill(geneta, genphi)
                    hc3.Fill(genen, genphi)

                elif round(deltaE_corr-deltaE,4) < -0.2:
                    hc4.Fill(genen, geneta)
                    hc5.Fill(geneta, genphi)
                    hc6.Fill(genen, genphi)
                """
    """
    ptest = [[[0.01,0.51,0.33,0.99],   #canvas0, pad0
              [0.34,0.51,0.66,0.99],   #canvas0, pad1
              [0.67,0.51,0.99,0.99],
              [0.01,0.01,0.33,0.49],   #canvas0, pad0
              [0.34,0.01,0.66,0.49],   #canvas0, pad1
              [0.67,0.01,0.99,0.49]]]   #canvas0, pad2

    with RootPlotting(ncanvas=1, npads=6, cdims=[[2000,1200]], pcoords=ptest) as plot:
        plot.plotHistogram(cpos=0, ppos=0, h=hc1, title=hc1.GetTitle(), draw_options='colz')
        plot.plotHistogram(cpos=0, ppos=1, h=hc2, title=hc2.GetTitle(), draw_options='colz')
        plot.plotHistogram(cpos=0, ppos=2, h=hc3, title=hc3.GetTitle(), draw_options='colz')
        plot.plotHistogram(cpos=0, ppos=3, h=hc4, title=hc4.GetTitle(), draw_options='colz')
        plot.plotHistogram(cpos=0, ppos=4, h=hc5, title=hc5.GetTitle(), draw_options='colz')
        plot.plotHistogram(cpos=0, ppos=5, h=hc6, title=hc6.GetTitle(), draw_options='colz')
        plot.save(cpos=0, name='pic')
    """

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
    elif FLAGS.mode == 2:
        if FLAGS.apply_weights:
            histos1 = histos[6:12]
            histos2 = histos[12:18]
        else:
            histos = histos[18:]
            print("LEN: ", len(histos))
    plotHistograms(histos, cdims, pcoords, 
                   os.path.join(FLAGS.outpath,picname))
    #plotHistograms(histos2, cdims, pcoords, 
    #               os.path.join(FLAGS.outpath,picname+'_right'))

if __name__ == "__main__":
    parser = Argparser.Argparser()
    FLAGS = parser.get_flags()
    parser.print_args()
    NREG=3
    NLAYERS=28
    A=[TMath.Pi()*1.3**2, TMath.Pi()*2.6**2, TMath.Pi()*5.3**2]
    main()
