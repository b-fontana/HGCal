import os, sys
import numpy as np
from UserCode.HGCalMaskResolutionAna import Argparser
from array import array as Carray
from collections import OrderedDict

from ROOT import TCanvas, TLatex, TFile, TMath, TH1F, TH2F, TLorentzVector, TF1
from ROOT import gStyle, gROOT, kTemperatureMap
from UserCode.HGCalMaskResolutionAna.RootTools import buildMedianProfile
from UserCode.HGCalMaskVisualProd.RootUtils import RootPlotting

parser = Argparser.Argparser()
FLAGS = parser.get_flags()
parser.print_args()

A=[TMath.Pi()*1.3**2, TMath.Pi()*2.6**2, TMath.Pi()*5.3**2]
NREG=3

def getEnergiesForCalibration(data,ireg):
    """
    Returns an array of reconstructed energies, eta and average noise
    """
    x=[]
    for i in range(0,data.GetEntriesFast()):
        data.GetEntry(i)
        for ia in range(1,2):
            genen=getattr(data,'genen%d'%ia)
            geneta=abs(getattr(data,'geneta%d'%ia))
            if genen<FLAGS.mingenen: continue
            if geneta<FLAGS.mingeneta: continue
            if geneta>FLAGS.maxgeneta: continue
            recen=getattr(data,'en%d_%d'%(ia,ireg))
            #avgnoise=getattr(data,'noise%d_%d'%(ia,ireg))
            avgnoise=getattr(data,'noise%d_3'%ia)*A[ireg-1]/A[2]
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
        hnames = ['res{}_signal', 'res{}_bckg', 
                  'en{}_signal',  'en{}_bckg']
        for ireg in range(1,NREG+1):
            histos[hnames[0].format(ireg)] = TH1F(hnames[0].format(ireg), ';#Delta E/E;PDF',
                                                  50, liminf, limsup)
            histos[hnames[1].format(ireg)] = TH1F(hnames[1].format(ireg), ';#Delta E/E;PDF',
                                                  50, liminf, limsup)
            histos[hnames[2].format(ireg)] = TH1F(hnames[1].format(ireg), ';E;PDF',
                                                  80, 0., 180.)
            histos[hnames[3].format(ireg)] = TH1F(hnames[1].format(ireg), ';E;PDF',
                                                  80, 0., 180.)

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

        #reconstructed photons in different regions
        for ireg in range(1,NREG+1):
            for ia in range(1,ncands+1):
                genen    = getattr(data,'genen{}'.format(ia))
                geneta   = getattr(data,'geneta{}'.format(ia))
                genphi   = getattr(data,'genphi{}'.format(ia))
                recen    = getattr(data,'en{}_{}'.format(ia,ireg))
                #avgnoise = getattr(data,'noise{}_{}'.format(ia,ireg))
                avgnoise = getattr(data,'noise{}_3'.format(ia))*A[ireg-1]/A[2]

                if 'L0' in calib:
                    recen=recen/(calib['L0'][ireg].Eval(abs(geneta))+1.0)
                    if 'L1' in calib:
                        recen=recen/(calib['L1'][ireg].Eval(recen)+1.0)
                        if 'L2' in calib and ireg in calib['L2']:
                            recen=recen-calib['L2'][ireg].Eval(avgnoise)

                deltaE = recen/genen-1.
                if FLAGS.ncomponents == 1:
                    histos[hnames[0].format(ireg)].Fill(deltaE)
                    histos[hnames[1].format(ireg)].Fill(geneta, genphi, deltaE)
                    histos[hnames[2].format(ireg)].Fill(geneta, genphi)
                if FLAGS.ncomponents == 2:
                    if abs(geneta) <= 2.95:
                        histos[hnames[0].format(ireg)].Fill(deltaE)
                        histos[hnames[2].format(ireg)].Fill(recen)
                    else:
                        histos[hnames[1].format(ireg)].Fill(deltaE)
                        histos[hnames[3].format(ireg)].Fill(recen)
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
        cdims = [[2000,2000]]
        correct_order = ['res1_signal', 'res2_signal', 'res3_signal',
                         'res1_bckg',   'res2_bckg',   'res3_bckg',
                         'en1_signal',  'en2_signal',  'en3_signal',
                         'en1_bckg',    'en2_bckg',    'en3_bckg']
        assert len(correct_order) == len(histos.keys())
        histos = [histos[correct_order[i]] for i in range(len(correct_order))]

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
            if h.Integral()==0.: continue
            h.Scale(1./h.Integral())
            plot.plotHistogram(cpos=0, ppos=ih, h=h, draw_options='colz')
            h.GetYaxis().SetRangeUser(0,h.GetMaximum()*1.2)
            tex = TLatex()
            tex.SetTextFont(42)
            tex.SetTextSize(0.04)
            tex.SetNDC()
            if ih < 3 or (ih > 5 and ih < 8):
                f = TF1('f1', 'gaus', 
                        h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
                f.SetParameters(1., 0., .15)
                f.SetParLimits(0, 0.001, 100.) #Constant 
                f.SetParLimits(1, -.5, .5) #Mean
                f.SetParLimits(2, 0.1, 10.) #Sigma
                plot.fitHistogram(f, h)
                gaus = h.GetListOfFunctions().At(0)
                tex.DrawLatex(0.12,0.92,'#bf{CMS} #it{simulation preliminary}')
                tex.DrawLatex(0.15,0.83,'SR{}'.format(ih+1))
                tex.DrawLatex(0.15,0.80,'#mu=%3.3f#pm%3.3f'
                              %(gaus.GetParameter(1),gaus.GetParError(1)))
                tex.DrawLatex(0.15,0.76,'#sigma=%3.3f#pm%3.3f'
                              %(gaus.GetParameter(2),gaus.GetParError(2)))
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
