import os, sys
import numpy as np
import pickle
from UserCode.HGCalMaskResolutionAna import Argparser
from array import array as Carray
from collections import OrderedDict

from ROOT import TCanvas, TLatex, TFile, TMath, TH1F
from ROOT import TLegend, TH2F, TLorentzVector, TProfile, TH1D
from ROOT import gStyle, gROOT, kTemperatureMap
from UserCode.HGCalMaskResolutionAna.RootTools import buildMedianProfile
from UserCode.HGCalMaskResolutionAna.SoftwareCorrection import IncompleteShowersCorrection
from UserCode.HGCalMaskVisualProd.RootPlotting import RootPlotting
from UserCode.HGCalMaskVisualProd.RootObjects import RootHistograms

parser = Argparser.Argparser()
FLAGS = parser.get_flags()
parser.print_args()

A=[TMath.Pi()*1.3**2, TMath.Pi()*2.6**2, TMath.Pi()*5.3**2]
NREG=3
NLAYERS=28

def getEnergiesForCalibration(data,ireg):
    """
    Returns an array of reconstructed energies, eta and average noise
    """
    x=[]
    for i in range(0,data.GetEntriesFast()):
        data.GetEntry(i)
        genen=getattr(data,'genen')
        geneta=abs(getattr(data,'geneta'))
        if genen<FLAGS.mingenen: continue
        if geneta<FLAGS.mingeneta: continue
        if geneta>FLAGS.maxgeneta: continue
        for il in range(1,NREG+1):
            recen = getattr(data,'en_sr{}_ROI'.format(ireg,il))
            avgnoise = getattr(data,'noise_sr3_layer{}'.format(il))*A[ireg-1]/A[2]
            x.append( [genen,geneta,recen,avgnoise] )
    return np.array(x)

def calibrateSpectrum(h,title,proc,func='pol1'):
    """
    Calibrates a DeltaE/E versus x spectrum based on the profile.
    """
    c=TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.1)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.1)
    #prof=h.ProfileX()
    prof=buildMedianProfile(h)
    h.Draw('colz')
    prof.Draw('e2p')
    prof.Fit(func)
    calibGr=prof.GetListOfFunctions().At(0).Clone(h.GetName()+title+'_calib')        
    tex=TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.SetNDC()
    tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')
    tex.DrawLatex(0.15,0.88,title)
    tex.SetTextAlign(31)
    tex.DrawLatex(0.97,0.96,proc)        
    #c.SaveAs(os.path.join(FLAGS.outpath,h.GetName()+title+'.png'))
    prof.Delete()        
    c.Delete()
    return calibGr

def doL0L1Calibration(url, calib, nq=6):
    """
    Performs the L0 (uniform eta response) and L1 (absolute scale) calibrations.
    """
    #derive L0 and L1 calibrations for 3 different signal regions
    fIn=TFile.Open(url+'.root')
    data=fIn.Get('data')
    for ireg in [1,2,3]:
        x=getEnergiesForCalibration(data,ireg)
        xq=np.percentile(x, [i*100./nq for i in range(0,nq+1)], axis=0)

        #relative calibration versus eta        
        resolVsEta = TH2F('resolvseta',';|#eta|;#DeltaE/E', nq, 
                          Carray('d',xq[:,1]), 100,-1,1)
        for i in range(0,len(x)):
            genEn, genEta, recEn,_ = x[i]
            deltaE=recEn/genEn-1.
            resolVsEta.Fill(genEta,deltaE)
        calib['L0'][ireg]=calibrateSpectrum(resolVsEta,'SR%d'%ireg,
                                            FLAGS.plotLabel+' (PU=0)','pol2')
    
        #relative calibration versus energy
        resolVsEn = TH2F('resolvsen',';Reconstructed energy [GeV];#DeltaE/E', nq, 
                         Carray('d',xq[:,0]), 100,-1,1)
        for i in range(0,len(x)):
            genEn,genEta,recEn,_=x[i]
            recEn=recEn/(calib['L0'][ireg].Eval(genEta)+1.0)
            deltaE=recEn/genEn-1.
            resolVsEn.Fill(recEn,deltaE)
            resolVsEta.Fill(genEta,deltaE)
        calib['L1'][ireg]=calibrateSpectrum(resolVsEn,'SR%d'%ireg,
                                            FLAGS.plotLabel+' (PU=0)','pol1')
    fIn.Close()
    return calib

def doPUCalibration(url, calib, nq=10):
    """
    Parametrizes the absolute shift in energy as function of the 
    average noise in the SR cone.
    """
    fIn=TFile.Open(url+'.root')
    data=fIn.Get('data')

    for ireg in [1,2,3]:
        x=getEnergiesForCalibration(data,ireg)
        xq=np.percentile(x, [i*100./nq for i in range(0,nq+1)], axis=0)

        #relative calibration versus eta
        resolVsNoiseFrac = TH2F('resolvsnoisefrac',';<Noise> [GeV];#DeltaE [GeV]', nq, 
                                Carray('d',xq[:,3]),50,-100,100)

        for i in range(0,len(x)):
            genEn,genEta,recEn,noiseEst=x[i]
            recEn=recEn/(calib['L0'][ireg].Eval(genEta)+1.0)
            recEn=recEn/(calib['L1'][ireg].Eval(recEn)+1.0)
            deltaE=recEn-genEn 
            resolVsNoiseFrac.Fill(noiseEst,deltaE)
    
        calib['L2'][ireg]=calibrateSpectrum(resolVsNoiseFrac,'SR%d'%ireg,
                                            FLAGS.plotLabel+' (PU=140)', 'pol2')
    return calib

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
                tex = plot.setLatex(ts=0.06)
                if FLAGS.samples == 'inner':
                    tex.DrawLatex(0.68,0.82,' |#eta| < '+str(FLAGS.etacuts[-2]))
                elif FLAGS.samples == 'outer':
                    tex.DrawLatex(0.68,0.92,' |#eta| > '+str(FLAGS.etacuts[-2]))

        elif FLAGS.mode == 2:
            for ih in range(NREG):
                if FLAGS.apply_weights:
                    plot.plotHistogram(cpos=0, ppos=ih, h=histos[ih], 
                                       lw=3,mc=4,msize=.5,lc=4, 
                                       draw_options='E')
                    tex = plot.setLatex()
                    plot.fitHistogram(h=histos[ih], fname='crystalball', frange=(-1.,1.),
                                      tex=tex)
                    if FLAGS.samples == 'inner':
                        tex.DrawLatex(0.58,0.92,'Inner radius;  SR{}'.format((ih%3)+1))
                    elif FLAGS.samples == 'outer':
                        tex.DrawLatex(0.58,0.92,'Outer radius;  SR{}'.format((ih%3)+1))
                    tex.DrawLatex(0.11,0.92,'#bf{CMS} #it{simulation preliminary}')
                    tex.SetTextAlign(31)
                    plot.plotHistogram(cpos=0, ppos=ih+3, h=histos[ih+3], 
                                       lw=3,mc=4,msize=.5,lc=4,
                                       draw_options='E')
                    tex = plot.setLatex()
                    plot.fitHistogram(h=histos[ih+3], fname='crystalball', frange=(-1.,1.),
                                      tex=tex)
                    if FLAGS.samples == 'inner':
                        tex.DrawLatex(0.58,0.92,'Inner radius;  SR{}'.format((ih%3)+1))
                    elif FLAGS.samples == 'outer':
                        tex.DrawLatex(0.58,0.92,'Outer radius;  SR{}'.format((ih%3)+1))
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

    if FLAGS.apply_weights:
        calibshowers_str = ( 'calibshowers_mask'+str(FLAGS.mask)+'_'+
                             FLAGS.samples+'_mode2')
        showercorr = IncompleteShowersCorrection(calibshowers_str+'.root')
        weights = showercorr.getCorrectionWeights()
        lowstats_factors = showercorr.calculateLowStatisticsFactor([8., 8., 8])

    fIn=TFile.Open(FLAGS.noPUFile+'.root')
    data=fIn.Get('data')

    calib=OrderedDict()
    calib['L0']={}
    calib['L1']={}
    doL0L1Calibration(url=FLAGS.noPUFile, calib=calib)
    if FLAGS.PUFile != '':
        calib['L2']={}
        doPUCalibration(url=FLAGS.PUFile, calib=calib)
        applyCalibrationTo(url=FLAGS.PUFile, calib=calib,
                           title=FLAGS.plotLabel+' (PU={})'.format(FLAGS.puTag))
        with open('calib_pu{}.pck'.format(FLAGS.puTag),'w') as cachefile:
            pickle.dump(calib,cachefile, pickle.HIGHEST_PROTOCOL)

    histos=OrderedDict()
    limsup, liminf = .6, -1.1
    etabins = 12
    etacuts = FLAGS.etacuts
    if FLAGS.samples == 'inner':
        phibins, etainf, etasup = 12, 2.69, 3.51
    elif FLAGS.samples == 'outer':
        phibins, etainf, etasup = 30, 1.34, 1.68
    else:
        raise ValueError("Specify a valid value for the '--samples' option.")

    if FLAGS.mode == 1:
        hn = ['den{}', 'den{}_2D_res', 'den{}_2D_events']
        for ireg in range(1,NREG+1):
            histos[hn[0].format(ireg)] = TH1F(hn[0].format(ireg),';#Delta E/E;PDF',
                                                  100, liminf, limsup)
            histos[hn[1].format(ireg)] = TH2F(hn[1].format(ireg), ';|#eta|;#phi',
                                                  50, etainf, etasup,
                                                  12, -TMath.Pi(), TMath.Pi())
            histos[hn[2].format(ireg)] = TH2F(hn[2].format(ireg), ';|#eta|;#phi',
                                                  50, etainf, etasup,
                                                  12, -TMath.Pi(), TMath.Pi())
    elif FLAGS.mode == 2:
        hn = ['res{}',                    'en{}', 
              'en{}_per_layer_signal',    'en{}_per_layer_bckg1',
              'en{}_per_layer_bckg2',     'en{}_per_layer_bckg3',
              'noise{}_per_layer_signal', 'noise{}_per_layer_bckg1',
              'noise{}_per_layer_bckg2',  'noise{}_per_layer_bckg3']
        for ireg in range(1,NREG+1):
            histos[hn[0].format(ireg)] = TH1F(hn[0].format(ireg), ';#Delta E/E_{gen};PDF',
                                              50, liminf, limsup)
            histos[hn[1].format(ireg)] = TH1F(hn[1].format(ireg), ';#Delta E/E_{gen};PDF',
                                              50, liminf, limsup)
            #histos[hn[1].format(ireg)] = TH1F(hn[1].format(ireg), ';E_{reco};PDF',
            #                                  80, 0., 145.)
            bins = Carray('d', np.arange(0.5,29,1.))
            assert len(bins) == NLAYERS+1
            strings = ';Layer;E_{reco} / E_{gen}'
            for ih in range(2,len(hn)):
                histos[hn[ih].format(ireg)] = TProfile(hn[ih].format(ireg), strings,
                                                       NLAYERS, bins)
    else:
        raise ValueError('The number of components has to be 1 or 2.')

    for h in histos:
        histos[h].Sumw2()
        histos[h].SetMarkerStyle(20)
        histos[h].SetDirectory(0)

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
            if FLAGS.mode == 1:
                if geneta < FLAGS.etacuts[-2]:
                    histos[hn[0].format(ireg)].Fill(deltaE)
                    histos[hn[1].format(ireg)].Fill(geneta, genphi, deltaE)
                    histos[hn[2].format(ireg)].Fill(geneta, genphi)

            elif FLAGS.mode == 2:
                ###Calculate and calibrate the energy per layer###
                recen_corr = 0 
                for il in range(1,NLAYERS+1):
                    #if: 'signal-like': complete showers
                    #elif: 'background-like': incomplete showers
                    b = histos[hn[2].format(ireg)].FindBin(il)
                    if bool_sig:
                        v = f1*getattr(data,'en_sr{}_layer{}'.format(ireg,il)) - f2
                        histos[hn[2].format(ireg)].Fill(b,v/genen)
                        if FLAGS.apply_weights:
                            recen_corr += v
                        v = (f1*getattr(data,'noise_sr3_layer{}'.format(il))
                             *A[ireg-1]/A[2] - f2)
                        histos[hn[6].format(ireg)].Fill(b,v/genen)
                    else:
                        for w in range(len(etacuts)-1):
                            if bools_bckg[w]:      
                                v = f1*getattr(data,'en_sr{}_layer{}'.format(ireg,il)) - f2
                                histos[hn[3+w].format(ireg)].Fill(b,v/genen)
                                if ( FLAGS.apply_weights and 
                                     weights[w + (ireg-1)*NREG][il-1]!=0):
                                    recen_corr += v/weights[w + (ireg-1)*NREG][il-1]
                                v = (f1*getattr(data,'noise_sr3_layer{}'.format(il))
                                     *A[ireg-1]/A[2] - f2)
                                histos[hn[7+w].format(ireg)].Fill(b,v/genen)
                if geneta < FLAGS.etacuts[-2]:
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
        pcoords = [[[0.01,0.505,0.33,0.99],   #canvas0, pad0
                    [0.34,0.505,0.66,0.99],   #canvas0, pad1
                    [0.67,0.505,0.99,0.99],   #canvas0, pad2
                    [0.01,0.01,0.33,0.495],   #canvas0, pad3
                    [0.34,0.01,0.66,0.495],   #canvas0, pad4
                    [0.67,0.01,0.99,0.495]]]  #canvas0, pad5
        cdims = [[2000,1200]]
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
            histos = histos[:6]
        else:
            histos = histos[6:]
    plotHistograms(histos, cdims, pcoords, 
                   os.path.join(FLAGS.outpath,picname))
    fOut=TFile.Open('calib{}.root'.format(''.join(calib.keys())),'RECREATE')
    for h in histos: 
        h.Write()
    fOut.Close()

if __name__ == "__main__":
    main()
