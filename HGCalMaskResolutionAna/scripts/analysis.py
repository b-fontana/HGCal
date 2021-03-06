import os, sys
import numpy as np
import pickle
from UserCode.HGCalMaskResolutionAna import Argparser
from array import array as Carray
from collections import OrderedDict

from ROOT import TCanvas, TLatex, TFile, TMath, TH1F
from ROOT import TLegend, TH2F, TLorentzVector, TProfile, TH1D, TGraphErrors
from ROOT import Double, gStyle, gROOT, kTemperatureMap
from UserCode.HGCalMaskVisualProd.SystemUtils import averageContiguousVals as Av
from UserCode.HGCalMaskVisualProd.SystemUtils import EtaStr as ES
from UserCode.HGCalMaskResolutionAna.SoftwareCorrection import IncompleteShowersCorrection
from UserCode.HGCalMaskResolutionAna.SoftwareCorrection import DifferentiateShowersByEnergy
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

    npads = len(pcoords[0])
    hdiv, stack = ([] for _ in range(2))
    with RootPlotting(ncanvas=1, npads=npads, cdims=cdims, pcoords=pcoords) as plot:
        titles = ['Resolution', 'Resolution vs Eta', 
                  'RMS vs Eta', 'Bias vs Eta', 'RMS/(1+Bias) vs Eta',
                  'Resolution', 'Nevents', 'Resolution / Nevents']
        if FLAGS.mode == 1:
            legends1 = [TLegend(0.12, 0.76, 0.44, 0.89) for  _ in range(3)]
            it = -1
            for ih in range(len(histos[:-3])):
                if ih%3==0: it += 1
                h = histos[ih]
                if ih<6:
                    plot.plotHistogram(cpos=0, ppos=ih, h=h,
                                       lw=3,mc=1,msize=.5,lc=1, 
                                       title=titles[it], draw_options='colz')
                    if ih<3:
                        tex = plot.setLatex(ts=0.04)
                        plot.fitHistogram(h=h, fname='crystalball', 
                                          frange=(-1.,1.), tex=tex)
                else:
                    plot.plotGraph(cpos=0, ppos=ih, g=h, 
                                   lw=3,mc=4,msize=.5,lc=4, 
                                   yranges=(-0.5,3.5),
                                   title=titles[it],
                                   draw_options='AP')

        elif FLAGS.mode == 2:
            if FLAGS.apply_weights:
                legends1 = [TLegend(0.12, 0.76, 0.44, 0.89) for  _ in range(3)]
            else:
                legends1 = [TLegend(0.56, 0.66, 0.86, 0.89) for  _ in range(3)]
                legends2 = [TLegend(0.68, 0.75, 0.86, 0.89) for _ in range(3)]

            for ih in range(NREG):
                if FLAGS.apply_weights:
                    for ixx in range(int(npads/NREG)):
                        idx = ih+NREG*ixx
                        if ixx == 0:
                            plot.plotHistogram(cpos=0, ppos=idx, h=histos[idx], 
                                               lw=3,mc=4,msize=.5,lc=4, 
                                               draw_options='colz')
                        elif ixx == 1:
                            plot.plotHistogram(cpos=0, ppos=idx, h=histos[idx], 
                                               lw=3,mc=4,msize=.5,lc=4, 
                                               title=titles[ixx],
                                               draw_options='colz')
                        else:
                            plot.plotGraph(cpos=0, ppos=idx, g=histos[idx], 
                                           lw=3,mc=4,msize=.5,lc=4, 
                                           yranges=(-0.5,3.5),
                                           title=titles[ixx],
                                           name = str(ixx*ih),
                                           draw_options='AP')
                        tex = plot.setLatex()
                        if ixx == 0:
                            plot.fitHistogram(h=histos[idx], fname='crystalball', 
                                              frange=(-1.,1.), tex=tex)   
                            if FLAGS.samples == 'inner':
                                tex.DrawLatex(0.58,0.92,'Inner radius;  SR{}'
                                              .format((ih%3)+1))
                            elif FLAGS.samples == 'outer':
                                tex.DrawLatex(0.58,0.92,'Outer radius;  SR{}'
                                              .format((ih%3)+1))

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
                    bckgcuts_str = [str(i) for i in bckgcuts]
                    if FLAGS.samples == 'inner':
                        tex.DrawLatex(0.58,0.92,'Inner radius;  SR{}'.format((ih%3)+1))
                    elif FLAGS.samples == 'outer':
                        tex.DrawLatex(0.58,0.92,'Outer radius;  SR{}'.format((ih%3)+1))

                    legends1[ih].AddEntry(histos[ih], 
                                          'Cumdiff < '+bckgcuts_str[0], 'L')
                    for it in range(len(bckgcuts_str)-1):
                        legends1[ih].AddEntry(histos[ih+3*(it+1)], 
                                    bckgcuts_str[it]+'< Cumdiff < '+bckgcuts_str[it+1], 'L')
                    legends1[ih].AddEntry(histos[ih+3*(it+2)], 
                                    'Cumdiff > '+bckgcuts_str[it+1], 'L')

                    tex.DrawLatex(0.11,0.92,'#bf{CMS} #it{simulation preliminary}')
                    tex.SetTextAlign(31)
                    legends1[ih].Draw()

                    hdiv.append(histos[ih+3].Clone('weight1_sr{}'.format(ih+1)))
                    hdiv.append(histos[ih+6].Clone('weight2_sr{}'.format(ih+1)))
                    hdiv.append(histos[ih+9].Clone('weight3_sr{}'.format(ih+1)))
                    for idiv in range(len(bckgcuts)):
                        hdiv[-3+idiv].Divide(histos[ih])
                        extrastr = '' if idiv==0 else 'same'
                        hdiv[-3+idiv].GetYaxis().SetRangeUser(0., 6.)
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
                    for iv in range(len(bckgcuts)):
                        legends2[ih].AddEntry(hdiv[iv], 'weight'+str(iv+1), 'L')
                        legends2[ih].Draw()
                    tex.DrawLatex(0.11,0.92,'#bf{CMS} #it{simulation preliminary}')
                    tex.SetTextAlign(31)

        plot.save(cpos=0, name=cname)

    if not FLAGS.apply_weights:
        save_str = base.paths.weights
        RootHistograms(histos).save(save_str)
        RootHistograms(hdiv).save(save_str, mode='UPDATE')

def main():
    #gStyle.SetOptStat(0)
    gROOT.SetBatch(True)
    gStyle.SetPalette(kTemperatureMap)

    fIn=TFile.Open(FLAGS.noPUFile)
    data=fIn.Get('data')
    calib_str = base.paths.calibrations_nopu
    """
    calibration = Calibration(FLAGS)
    calibration.nopu_calibration()
    calibration.save(calib_str)
    """
    with open(calib_str, 'r') as cachefile:
        calib = pickle.load(cachefile)

    if FLAGS.apply_weights:
        calibshowers_str = base.paths.weights
        bckgcuts_extended = np.append(bckgcuts, 0.9)
        showercorr = IncompleteShowersCorrection(calibshowers_str,
                                                 discrvals=Av(bckgcuts_extended))
        weights = showercorr.CorrectionWeights()
        if FLAGS.samples == 'inner':
            boundaries = [5, 5, 5]
            corr_mode = 'left' 
        elif FLAGS.samples == 'outer':
            boundaries = [23, 23, 23]
            corr_mode = 'right' 
        lowstats_factors = showercorr.calculateLowStatisticsFactor(boundaries, corr_mode)
        weights_graphs = [showercorr.buildCorrectionWeightsGraphs(region=i+1)
                          for i in range(NREG)]

    histos=OrderedDict()
    enbins, eninf, ensup = 200, -2.01, 1.99
    if FLAGS.samples == 'inner':
        phibins, etabins, etainf, etasup = 12, 10, 2.69, 3.04
    elif FLAGS.samples == 'outer':
        phibins, etabins, etainf, etasup = 12, 10, 1.44, 1.66

    if FLAGS.mode == 1:
        hn = ['den{}', 'den_eta{}', 'den{}_2D_res', 'den{}_2D_events']
        for ireg in range(1,NREG+1):
            bins = Carray('d', np.arange(-1.05, .8, 0.01))
            strings = ';#Delta E/E_{gen};PDF'
            histos[hn[0].format(ireg)] = TH1F(hn[0].format(ireg), strings,
                                              len(bins)-1, bins)
            histos[hn[1].format(ireg)] = TH2F(hn[1].format(ireg), ';|#eta|;#Delta E/E',
                                              etabins, etainf, etasup,
                                              enbins, eninf, ensup)
            histos[hn[2].format(ireg)] = TH2F(hn[2].format(ireg), ';|#eta|;#phi',
                                              50, etainf, etasup,
                                              phibins, -TMath.Pi(), TMath.Pi())
            histos[hn[3].format(ireg)] = TH2F(hn[3].format(ireg), ';|#eta|;#phi',
                                              50, etainf, etasup,
                                              phibins, -TMath.Pi(), TMath.Pi())
    elif FLAGS.mode == 2:
        fracEn = np.zeros((NREG,NLAYERS), dtype=float)
        countfracEn = np.zeros((NREG,NLAYERS), dtype=int)
        
        for i in range(0, data.GetEntriesFast()):
            data.GetEntry(i)
            genen    = getattr(data,'genen')        
            geneta   = abs(getattr(data,'geneta'))
            genphi   = getattr(data,'genphi')

            for ireg in range(1,NREG+1):
                recen = getattr(data,'en_sr{}_ROI'.format(ireg))
                avgnoise = getattr(data,'noise_sr3_ROI')*A[ireg-1]/A[2]

                #Calibration factors. f2 is used for PU.
                f1, f2 = 1., 0.
                etaregions_shifted = np.roll(etaregions, shift=-1)[:-1]
                for ieta1,ieta2 in zip(etaregions[:-1], etaregions_shifted):
                    #in case it lies outside the limits of the calibration
                    #the event is calibrated with the full calibration region
                    if geneta < etaregions[0] or geneta > etaregions[-1]:
                        idstr = 'sr{}_from{}to{}'.format(ireg,
                                                         ES(etaregions[0]), 
                                                         ES(etaregions[-1]))
                    elif (geneta < ieta1 or geneta >= ieta2): 
                        continue
                    else:
                        idstr = 'sr{}_from{}to{}'.format(ireg, ES(ieta1), ES(ieta2))
                    if 'L0' in calib:
                        f1 /= calib['L0'][idstr].Eval(geneta)+1.0
                        if 'L1' in calib:
                            f1 /= calib['L1'][idstr].Eval(f1*recen)+1.0    
                            if 'L2' in calib and ireg in calib['L2']:
                                f2 = calib['L2'][idstr].Eval(avgnoise)
                    recen = f1*recen - f2 
                for il in range(1,NLAYERS+1):
                    v = f1*getattr(data,'en_sr{}_layer{}'.format(ireg,il)) - f2
                    if ( (FLAGS.samples == "inner" and geneta < etaregions[0]+0.05 or
                          FLAGS.samples == "outer" and geneta > etaregions[-1]-0.05) and
                         recen != 0 ):
                        fracEn[ireg-1,il-1] += v / recen
                        countfracEn[ireg-1,il-1] += 1
        fracEn /= countfracEn

        hn = ['res_complete_before{}',          'res_complete_after{}', 
              'res_incomplete_before{}',        'res_incomplete_after{}',
              'res_total_before{}',             'res_total_after{}',
              'res_vs_eta_before{}',            'res_vs_eta_after{}',
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

            histos[hn[6].format(ireg)] = TH2F(hn[6].format(ireg), ';|#eta|;#Delta E/E',
                                              etabins, etainf, etasup,
                                              enbins, eninf, ensup)
            histos[hn[7].format(ireg)] = TH2F(hn[7].format(ireg), ';|#eta|;#Delta E/E',
                                              etabins, etainf, etasup,
                                              enbins, eninf, ensup)

            bins = Carray('d', np.arange(0.5,29,1.))
            strings = ';Layer;E_{reco} / E_{gen}'
            for ih in range(8,16):
                histos[hn[ih].format(ireg)] = TProfile(hn[ih].format(ireg), strings,
                                                       len(bins)-1, bins)
    for h in histos:
        histos[h].Sumw2()
        histos[h].SetMarkerStyle(20)
        histos[h].SetDirectory(0)

    for i in range(0, data.GetEntriesFast()):
        data.GetEntry(i)
        genen    = getattr(data,'genen')        
        geneta   = abs(getattr(data,'geneta'))
        genphi   = getattr(data,'genphi')

        for ireg in range(1,NREG+1):
            recen = getattr(data,'en_sr{}_ROI'.format(ireg))
            avgnoise = getattr(data,'noise_sr3_ROI')*A[ireg-1]/A[2]

            #Calibration factors. f2 is used for PU.
            f1, f2 = 1., 0.
            etaregions_shifted = np.roll(etaregions, shift=-1)[:-1]
            for ieta1,ieta2 in zip(etaregions[:-1], etaregions_shifted):
                if geneta < etaregions[0] or geneta > etaregions[-1]:
                    idstr = 'sr{}_from{}to{}'.format(ireg, 
                                                     ES(etaregions[0]), ES(etaregions[-1]))
                elif (geneta < ieta1 or geneta > ieta2): 
                    continue
                else:
                    idstr = 'sr{}_from{}to{}'.format(ireg, ES(ieta1), ES(ieta2))
                if 'L0' in calib:
                    f1 /= calib['L0'][idstr].Eval(geneta)+1.0
                    if 'L1' in calib:
                        f1 /= calib['L1'][idstr].Eval(f1*recen)+1.0    
                        if 'L2' in calib and ireg in calib['L2']:
                            f2 = calib['L2'][idstr].Eval(avgnoise)
            assert f1 != 1.
            recen = f1*recen - f2 
            deltaE = recen/genen-1.

            ###Store the energy resolution###
            if FLAGS.mode == 1:
                if deltaE > -1:
                    histos[hn[0].format(ireg)].Fill(deltaE)
                    histos[hn[1].format(ireg)].Fill(geneta, deltaE)
                histos[hn[2].format(ireg)].Fill(geneta, genphi, deltaE)
                histos[hn[3].format(ireg)].Fill(geneta, genphi)

            elif FLAGS.mode == 2:
                #differentiate complete from incomplete showers
                ROI_en = np.zeros((NLAYERS), dtype=float)
                for il in range(1,NLAYERS+1):
                    v = f1*getattr(data,'en_sr{}_layer{}'.format(ireg,il)) - f2
                    try:
                        ROI_en[il-1] = v/recen
                    except ZeroDivisionError:
                        ROI_en[il-1] = 0.

                lshift = [1., 1., 1.] if FLAGS.samples == 'outer' else [.65, .59, .48] #layer shift
                assert len(bckgcuts) == len(lshift)
                showerid = DifferentiateShowersByEnergy(ROI_en, fracEn[ireg-1,:], 
                                                        thresholds=bckgcuts, min_val=0.05)

                ###Calculate andc calibrate the energy per layer###
                recen_corr = 0 
                for il in range(1,NLAYERS+1):
                    if FLAGS.apply_weights:
                        if FLAGS.samples == 'inner':
                            weight_limit = il > boundaries[ireg-1] 
                        else:
                            weight_limit = il < boundaries[ireg-1] 
                    
                    b = histos[hn[8].format(ireg)].FindBin(il)
                    if showerid==0: #complete shower
                        v = f1*getattr(data,'en_sr{}_layer{}'.format(ireg,il)) - f2
                        try:
                            histos[hn[8].format(ireg)].Fill(b,v/recen)
                        except ZeroDivisionError:
                            histos[hn[8].format(ireg)].Fill(b,0.)
                        if FLAGS.apply_weights:
                            recen_corr += v
                        v = (f1*getattr(data,'noise_sr3_layer{}'.format(il))
                             *A[ireg-1]/A[2] - f2)
                        try:
                            histos[hn[12].format(ireg)].Fill(b,v/recen)
                        except ZeroDivisionError:
                            histos[hn[12].format(ireg)].Fill(b,0.)
                    else:
                        w = showerid-1
                        v = f1*getattr(data,'en_sr{}_layer{}'.format(ireg,il)) - f2
                        try:
                            histos[hn[9+w].format(ireg)].Fill(b*lshift[w],v/recen)
                        except ZeroDivisionError:
                            histos[hn[9+w].format(ireg)].Fill(b*lshift[w],0.)
                        if ( FLAGS.apply_weights and 
                             weights[ireg-1][w][il-1]!=0 and
                             weight_limit):
                            recen_corr += v/weights[ireg-1][w][int(round((il-1)*lshift[w],0))]
                            #weight_graphs[ireg][il].SetBit(weight_graphs[ireg][il].klsSortedX)
                            #weight_graphs[ireg][il].Eval(geneta, spline=0, 'S')
                        v = (f1*getattr(data,'noise_sr3_layer{}'.format(il))
                             *A[ireg-1]/A[2] - f2)
                        try:
                            histos[hn[13+w].format(ireg)].Fill(b,v/recen)
                        except ZeroDivisionError:
                            histos[hn[13+w].format(ireg)].Fill(b,0.)

                if FLAGS.apply_weights and FLAGS.method == 'ed':
                    if showerid==0: #complete shower
                        deltaE_corr = recen_corr/genen-1.
                        histos[hn[0].format(ireg)].Fill(deltaE)
                        histos[hn[1].format(ireg)].Fill(deltaE_corr)
                    else:
                        recen_corr *= (1 / (1-lowstats_factors[ireg-1]) )
                        distshift = 0.09 if FLAGS.samples == 'inner' else 0.08
                        recen_corr *= 1/distshift
                        deltaE_corr = recen_corr/genen-1.
                        if deltaE>-.95 and deltaE<-0.1:
                            histos[hn[2].format(ireg)].Fill(deltaE)
                            histos[hn[3].format(ireg)].Fill(deltaE_corr)
                    histos[hn[6].format(ireg)].Fill(geneta, deltaE)
                    histos[hn[7].format(ireg)].Fill(geneta, deltaE_corr)
    #end of tree loop    
    fIn.Close()

    if FLAGS.mode == 1:
        pcoords = [[[0.01,0.755,0.33,0.995],   
                    [0.34,0.755,0.66,0.995],   
                    [0.67,0.755,0.99,0.995],
                    [0.01,0.505,0.33,0.745],   
                    [0.34,0.505,0.66,0.745],
                    [0.67,0.505,0.99,0.745],
                    [0.01,0.255,0.33,0.495],
                    [0.34,0.255,0.66,0.495],
                    [0.67,0.255,0.99,0.495],
                    [0.01,0.005,0.33,0.245],
                    [0.34,0.005,0.66,0.245],
                    [0.67,0.005,0.99,0.245]]]
        cdims = [[1600,2000]]
        picname = '1comp_'+FLAGS.samples+'_'+FLAGS.method
        if FLAGS.apply_weights:
            picname += '_corrected' 
    else:
        pcoords = [[[0.01,0.51,0.33,0.99],
                    [0.34,0.51,0.66,0.99],  
                    [0.67,0.51,0.99,0.99],
                    [0.01,0.01,0.33,0.49],   
                    [0.34,0.01,0.66,0.49],   
                    [0.67,0.01,0.99,0.49]]]
        cdims = [[1000,600]]
        picname = '2comp_'+FLAGS.samples+'_'+FLAGS.method
        if FLAGS.apply_weights:
            picname += '_corrected' 
    
    correct_order = []
    for i in range(len(hn)):
        for ireg in range(1,NREG+1):
            correct_order.append(hn[i].format(ireg))
    assert len(correct_order) == len(histos.keys())
    histos = [histos[correct_order[i]] for i in range(len(correct_order))]

    if FLAGS.mode == 1:
        histos.append(histos[6].Clone())
        histos.append(histos[7].Clone())
        histos.append(histos[8].Clone())
        histos[-3].Divide(histos[9])
        histos[-2].Divide(histos[10])
        histos[-1].Divide(histos[11])

        htmp = []
        for ireg in range(NREG):
            h = histos[3+ireg]
            xbins, exbins, rms, erms, bias, ebias = ([] for _ in range(6))
            for xbin in xrange(1,h.GetNbinsX()+1):
                tmp = h.ProjectionY('tmp', xbin, xbin)
                xbins.append(h.GetXaxis().GetBinCenter(xbin))
                horizerror = ( h.GetXaxis().GetBinCenter(1) - 
                               h.GetXaxis().GetBinLowEdge(1) )
                exbins.append( horizerror )
                rms.append(tmp.GetRMS())
                erms.append(tmp.GetRMSError())
                """
                xq = Carray('d', [0.16,0.5,0.84]) 
                yq = Carray('d', [0.0,0.0,0.0 ])
                tmp.GetQuantiles(3,yq,xq)
                bias.append(yq[1])
                """
                bias.append(tmp.GetMean())
                #ebias.append((yq[0]+yq[2])/2)
                ebias.append(tmp.GetMeanError())
                tmp.Delete()
            xbins, exbins = np.array(xbins), np.array(exbins)
            rms, erms = np.array(rms), np.array(erms)
            bias, ebias = np.array(bias), np.array(ebias)
            indep = rms/(1.+bias)
            eindep = indep * np.sqrt( erms**2/rms**2 + ebias**2/(1+bias**2)  )
            htmp.append( TGraphErrors(etabins, xbins, rms, exbins, erms) )
            htmp.append( TGraphErrors(etabins, xbins, bias, exbins, ebias) )
            htmp.append( TGraphErrors(etabins, xbins, indep, exbins, eindep) )

        ht_tmp = [htmp[0],htmp[3],htmp[6],htmp[1],htmp[4],htmp[7],htmp[2],htmp[5],htmp[8]]
        histos = histos[:6] + ht_tmp
        if FLAGS.method == 'fineeta':
            htitles = ['rmsVSeta1',  'rmsVSeta2',  'rmsVSeta3',
                       'biasVSeta1', 'biasVSeta2', 'biasVSeta3',
                       'indepVSeta1','indepVSeta2','indepVSeta3']
            indices = [1, 2, 3] * 3
            indices.sort()
            fOut = TFile( base.paths.plots, 'RECREATE' )
            fOut.cd()
            for ih,h in enumerate(ht_tmp):
                h.SetName(htitles[ih])
                h.Write(htitles[ih])
            fOut.Write()
            fOut.Close()

        plotHistograms(histos, cdims, pcoords, 
                       os.path.join(FLAGS.outpath,picname+'.png'))

    elif FLAGS.mode == 2:
        if FLAGS.apply_weights:
            fOut = TFile( base.paths.plots, 'RECREATE' )
            fOut.cd()
            for ireg in range(NREG):
                str1 = hn[4].format(ireg+1)
                str2 = hn[5].format(ireg+1)
                histos[12+ireg] = histos[ireg].Clone(str1)
                histos[15+ireg] = histos[3+ireg].Clone(str2)      
                histos[12+ireg].Add(histos[6+ireg])
                histos[15+ireg].Add(histos[9+ireg])
            for h in histos:
                h.Write()
            histos_complete = histos[:6]
            histos_incomplete = histos[6:12]
            histos_total = histos[12:18]
            histos_res2D_before = histos[21:24]
            histos_res2D_after = histos[21:24]
            plotHistograms(histos_complete, cdims, pcoords, 
                           os.path.join(FLAGS.outpath,picname+'_complete.png'))
            plotHistograms(histos_incomplete, cdims, pcoords, 
                           os.path.join(FLAGS.outpath,picname+'_incomplete.png'))

            ht = histos_total[3:] + histos_res2D_after #res 1D + res 2D
            ht_tmp = []
            for ireg in range(NREG):
                h = ht[3+ireg]
                xbins, exbins, rms, erms, bias, ebias = ([] for _ in range(6))
                for xbin in xrange(1,h.GetNbinsX()+1):
                    tmp = h.ProjectionY('tmp', xbin, xbin)
                    xbins.append(h.GetXaxis().GetBinCenter(xbin))
                    horizerror = ( h.GetXaxis().GetBinCenter(1) - 
                                   h.GetXaxis().GetBinLowEdge(1) )
                    exbins.append( horizerror )
                    rms.append(tmp.GetRMS())
                    erms.append(tmp.GetRMSError())
                    """
                    xq = Carray('d', [0.16,0.5,0.84]) 
                    yq = Carray('d', [0.0,0.0,0.0 ])
                    tmp.GetQuantiles(3,yq,xq)
                    bias.append(yq[1])
                    """
                    bias.append(tmp.GetMean())
                    #ebias.append((yq[0]+yq[2])/2)
                    ebias.append(tmp.GetMeanError())
                    tmp.Delete()
                xbins, exbins = np.array(xbins), np.array(exbins)
                rms, erms = np.array(rms), np.array(erms)
                bias, ebias = np.array(bias), np.array(ebias)
                indep = rms/(1.+bias)
                eindep = indep * np.sqrt( erms**2/rms**2 + ebias**2/(1+bias**2)  )
                ht_tmp.append( TGraphErrors(etabins, xbins, rms, exbins, erms) )
                ht_tmp.append( TGraphErrors(etabins, xbins, bias, exbins, ebias) )
                ht_tmp.append( TGraphErrors(etabins, xbins, indep, exbins, eindep) )

            fOut.cd()
            ht_tmp = [ht_tmp[-9],ht_tmp[-6],ht_tmp[-3],
                      ht_tmp[-8],ht_tmp[-5],ht_tmp[-2],
                      ht_tmp[-7],ht_tmp[-4],ht_tmp[-1]]
            ht_titles = ['rmsVSeta1',  'rmsVSeta2',  'rmsVSeta3',
                         'biasVSeta1', 'biasVSeta2', 'biasVSeta3',
                         'indepVSeta1','indepVSeta2','indepVSeta3']
            indices = [1, 2, 3] * 3
            indices.sort()
            for ih,h in enumerate(ht_tmp):
                h.SetName(ht_titles[ih])
                h.Write(ht_titles[ih])
            pcoords = [[[0.01,0.805,0.33,0.995],   
                        [0.34,0.805,0.66,0.995],   
                        [0.67,0.805,0.99,0.995],
                        [0.01,0.605,0.33,0.795],   
                        [0.34,0.605,0.66,0.795],
                        [0.67,0.605,0.99,0.795],
                        [0.01,0.405,0.33,0.595],
                        [0.34,0.405,0.66,0.595],
                        [0.67,0.406,0.99,0.595],
                        [0.01,0.205,0.33,0.395],
                        [0.34,0.205,0.66,0.395],
                        [0.67,0.205,0.99,0.395],
                        [0.01,0.005,0.33,0.195],
                        [0.34,0.005,0.66,0.195],
                        [0.67,0.005,0.99,0.195]]]
            cdims = [[1600,2000]]
            ht += ht_tmp
            plotHistograms(ht, cdims, pcoords, 
                           os.path.join(FLAGS.outpath,picname+'_total.png'))
            fOut.Write()
            fOut.Close()

        else:
            histos = histos[24:]
            plotHistograms(histos, cdims, pcoords, 
                           os.path.join(FLAGS.outpath,picname+'.png'))

if __name__ == "__main__":
    parser = Argparser.Argparser()
    FLAGS = parser.get_flags()
    parser.print_args()
    if FLAGS.apply_weights and FLAGS.mode != 2:
        raise ValueError('The weights can only be used when mode==2.')
    base = PartialWafersStudies(FLAGS)
    NREG, NLAYERS, A = base.nsr, base.nlayers, base.sr_area
    etaregions = base.etaregions
    if FLAGS.method == 'ed':
        bckgcuts = np.array(FLAGS.bckgcuts)
    main()
