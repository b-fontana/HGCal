import os, sys
import numpy as np
from UserCode.HGCalMaskResolutionAna import Argparser
from array import array as Carray
from collections import OrderedDict

from ROOT import TCanvas, TLatex, TFile, TMath, TH1F, TH2F, TLorentzVector, TF1, TH1D
from ROOT import gStyle, gROOT, kTemperatureMap
from UserCode.HGCalMaskResolutionAna.RootTools import buildMedianProfile
from UserCode.HGCalMaskVisualProd.RootUtils import RootPlotting

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
    c.SaveAs(os.path.join(FLAGS.outpath,h.GetName()+title+'.png'))
    prof.Delete()        
    c.Delete()
    return calibGr

def doL0L1Calibration(url, calib, nq=6):
    """
    Performs the L0 (uniform eta response) and L1 (absolute scale) calibrations.
    """
    #derive L0 and L1 calibrations for 3 different signal regions
    fIn=TFile.Open(url)
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
    fIn=TFile.Open(url)
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

def CalibratedResolution(url, calib, title, ncomponents=1, ncands=1):
    """
    Applies the calibration to the photons and shows the energy resolution.
    """
    pfix=''.join(calib.keys())
    fIn=TFile.Open(url)
    data=fIn.Get('data')

    histos=OrderedDict()
    limsup, liminf = 4., -2.
    etabins = 12
    etacut = 1.7#2.95
    en_depo_sig = [[0. for _ in range(NLAYERS)] for _ in range(NREG)] 
    en_depo_bckg = [[0. for _ in range(NLAYERS)] for _ in range(NREG)] 
    avgnoise_depo_sig = [[0. for _ in range(NLAYERS)] for _ in range(NREG)] 
    avgnoise_depo_bckg = [[0. for _ in range(NLAYERS)] for _ in range(NREG)] 
    if FLAGS.samples == 'inner':
        phibins, etainf, etasup = 12, 2.75, 3.15
    elif FLAGS.samples == 'outer':
        phibins, etainf, etasup = 30, 1.42, 1.65
    else:
        raise ValueError("Specify a valid value for the '--samples' option.")

    if FLAGS.ncomponents == 1:
        hnames = ['den{}', 'den{}_2D_res', 'den{}_2D_events']
        for ireg in range(1,NREG+1):
            histos[hnames[0].format(ireg)] = TH1F(hnames[0].format(ireg),';#Delta E/E;PDF',
                                                  50, liminf, limsup)
            histos[hnames[1].format(ireg)] = TH2F(hnames[1].format(ireg), ';|#eta|;#phi',
                                                  50, etainf, etasup,
                                                  12, -TMath.Pi(), TMath.Pi())
            histos[hnames[2].format(ireg)] = TH2F(hnames[2].format(ireg), ';|#eta|;#phi',
                                                  50, etainf, etasup,
                                                  12, -TMath.Pi(), TMath.Pi())
    elif FLAGS.ncomponents == 2:
        hnames = ['res{}_signal',             'res{}_bckg', 
                  'en{}_signal',              'en{}_bckg',
                  'en{}_per_layer_signal',    'en{}_per_layer_bckg',
                  'noise{}_per_layer_signal', 'noise{}_per_layer_bckg']
        for ireg in range(1,NREG+1):
            histos[hnames[0].format(ireg)] = TH1F(hnames[0].format(ireg), ';#Delta E/E;PDF',
                                                  50, liminf, limsup)
            histos[hnames[1].format(ireg)] = TH1F(hnames[1].format(ireg), ';#Delta E/E;PDF',
                                                  50, liminf, limsup)
            histos[hnames[2].format(ireg)] = TH1F(hnames[2].format(ireg), ';E;PDF',
                                                  80, 0., 180.)
            histos[hnames[3].format(ireg)] = TH1F(hnames[3].format(ireg), ';E;PDF',
                                                  80, 0., 180.)
            bins = Carray('d', np.arange(0.5,29,1.))
            assert len(bins) == NLAYERS+1
            histos[hnames[4].format(ireg)] = TH1F(hnames[4].format(ireg), ';Layer;Energy',
                                                  NLAYERS, bins)
            histos[hnames[5].format(ireg)] = TH1F(hnames[5].format(ireg), ';Layer;Energy',
                                                  NLAYERS, bins)
            histos[hnames[6].format(ireg)] = TH1F(hnames[6].format(ireg), ';Layer;Energy',
                                                  NLAYERS, bins)
            histos[hnames[7].format(ireg)] = TH1F(hnames[7].format(ireg), ';Layer;Energy',
                                                  NLAYERS, bins)
    else:
        raise ValueError('The number of components has to be 1 or 2.')

    for h in histos:
        histos[h].Sumw2()
        histos[h].SetLineColor(1)
        histos[h].SetMarkerColor(1)
        histos[h].SetMarkerStyle(20)
        histos[h].SetDirectory(0)

    for i in range(0,data.GetEntriesFast()):
        data.GetEntry(i)
        genen    = getattr(data,'genen')        
        geneta   = getattr(data,'geneta')
        genphi   = getattr(data,'genphi')
        for ireg in range(1,NREG+1):
            recen = getattr(data,'en_sr{}_ROI'.format(ireg))
            avgnoise = getattr(data,'noise_sr3_ROI')*A[ireg-1]/A[2]

            if 'L0' in calib:
                recen /= calib['L0'][ireg].Eval(abs(geneta))+1.0
                if 'L1' in calib:
                    recen /= calib['L1'][ireg].Eval(recen)+1.0
                    if 'L2' in calib and ireg in calib['L2']:
                        recen -= calib['L2'][ireg].Eval(avgnoise)

            deltaE = recen/genen-1.
            ###Store the energy resolution###
            if FLAGS.ncomponents == 1:
                histos[hnames[0].format(ireg)].Fill(deltaE)
                histos[hnames[1].format(ireg)].Fill(geneta, genphi, deltaE)
                histos[hnames[2].format(ireg)].Fill(geneta, genphi)
            if FLAGS.ncomponents == 2:
                if abs(geneta) > etacut:
                    histos[hnames[0].format(ireg)].Fill(deltaE)
                    histos[hnames[2].format(ireg)].Fill(recen)
                else:
                    histos[hnames[1].format(ireg)].Fill(deltaE)
                    histos[hnames[3].format(ireg)].Fill(recen)
                    
                ###Calculate the energy per layer###
                for il in range(1,NLAYERS+1):
                    if abs(geneta) > etacut:
                        en_depo_sig[ireg-1][il-1] += getattr(data,'en_sr{}_layer{}'
                                                      .format(ireg,il))
                        avgnoise_depo_sig[ireg-1][il-1] += getattr(data,'noise_sr3_layer{}'
                                                      .format(il))*A[ireg-1]/A[2]
                    elif abs(geneta) <= etacut:
                        en_depo_bckg[ireg-1][il-1] += getattr(data,'en_sr{}_layer{}'
                                                      .format(ireg,il))
                        avgnoise_depo_bckg[ireg-1][il-1] += getattr(data,'noise_sr3_layer{}'
                                                      .format(il))*A[ireg-1]/A[2]

                

    #end of tree loop

    ###Store the energy per layer###
    for ireg in range(1,NREG+1):
        for il in range(1,NLAYERS+1):
            if FLAGS.ncomponents == 1:
                pass
            if FLAGS.ncomponents == 2:
                b = histos[hnames[4].format(ireg)].FindBin(il)
                histos[hnames[4].format(ireg)].Fill(b, en_depo_sig[ireg-1][il-1])
                histos[hnames[6].format(ireg)].Fill(b, avgnoise_depo_sig[ireg-1][il-1])
                histos[hnames[5].format(ireg)].Fill(b, en_depo_bckg[ireg-1][il-1])
                histos[hnames[7].format(ireg)].Fill(b, avgnoise_depo_bckg[ireg-1][il-1])
        
        
    fIn.Close()

    if FLAGS.ncomponents == 1:
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
    else:
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
        cdims = [[1100,1400]]
    
    correct_order = []
    for i in range(len(hnames)):
        for ireg in range(1,NREG+1):
            correct_order.append(hnames[i].format(ireg))
    assert len(correct_order) == len(histos.keys())
    histos = [histos[correct_order[i]] for i in range(len(correct_order))]
    if FLAGS.ncomponents == 1:
        histos.append(histos[3].Clone())
        histos.append(histos[4].Clone())
        histos.append(histos[5].Clone())
        histos[-3].Divide(histos[6])
        histos[-2].Divide(histos[7])
        histos[-1].Divide(histos[8])
    elif FLAGS.ncomponents == 2:
        histos = histos[12:]
    plotHistograms(histos, cdims, pcoords, 'pic')
    fOut=TFile.Open('calib{}.root'.format(pfix),'RECREATE')
    for h in histos: 
        h.Write()
    fOut.Close()

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

    with RootPlotting(ncanvas=1, npads=len(histos), cdims=cdims, pcoords=pcoords) as plot:
        for ih in range(len(histos)):
            h = histos[ih]
            #if h.Integral()==0.: continue
            #h.Scale(1./h.Integral())
            plot.plotHistogram(cpos=0, ppos=ih, h=h, draw_options='colz')
            h.GetYaxis().SetRangeUser(0,h.GetMaximum()*1.2)
            tex = TLatex()
            tex.SetTextFont(42)
            tex.SetTextSize(0.04)
            tex.SetNDC()
            if ih < 3 or (ih > 5 and ih < 8):
                """
                f = TF1('f1', 'gaus', 
                        h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
                f.SetParameters(1., 0., .15)
                f.SetParLimits(0, 0.001, 100.) #Constant 
                f.SetParLimits(1, -.5, .5) #Mean
                f.SetParLimits(2, 0.1, 10.) #Sigma
                plot.fitHistogram(f, h)
                gaus = h.GetListOfFunctions().At(0)
                """
                tex.DrawLatex(0.12,0.92,'#bf{CMS} #it{simulation preliminary}')
                tex.DrawLatex(0.15,0.83,'SR{}'.format((ih%3)+1))
                #tex.DrawLatex(0.15,0.80,'#mu=%3.3f#pm%3.3f'
                #              %(gaus.GetParameter(1),gaus.GetParError(1)))
                #tex.DrawLatex(0.15,0.76,'#sigma=%3.3f#pm%3.3f'
                #              %(gaus.GetParameter(2),gaus.GetParError(2)))
                tex.DrawLatex(0.7,0.83,'|#eta| #leq 2.95')
            else:
                tex.DrawLatex(0.12,0.92,'#bf{CMS} #it{simulation preliminary}')
                tex.DrawLatex(0.7,0.83,'|#eta| > 2.95')
            tex.SetTextAlign(31)

            #cutoff = 50000
            #if h1.GetMaximum()>cutoff:
            #    h1.SetMaximum(cutoff)
        plot.save(cpos=0, name=cname)

def main():
    gStyle.SetOptStat(0)
    gROOT.SetBatch(True)
    gStyle.SetPalette(kTemperatureMap)

    calib=OrderedDict()
    calib['L0']={}
    calib['L1']={}
    doL0L1Calibration(url=FLAGS.noPUFile, calib=calib)

    
    CalibratedResolution(url=FLAGS.noPUFile, calib=calib, 
                         title=FLAGS.plotLabel+' (PU=0)', 
                         ncomponents=FLAGS.ncomponents, ncands=1)
    
    if FLAGS.PUFile != '':
        calib['L2']={}
        doPUCalibration(url=FLAGS.PUFile, calib=calib)
        applyCalibrationTo(url=FLAGS.PUFile, calib=calib,
                           title=FLAGS.plotLabel+' (PU={})'.format(FLAGS.puTag))

        #save final calibration
        import pickle
        with open('calib_pu{}.pck'.format(FLAGS.puTag),'w') as cachefile:
            pickle.dump(calib,cachefile, pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    main()
