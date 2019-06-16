import os, sys
import numpy as np
from UserCode.HGCalMaskResolutionAna import Argparser
from array import array as Carray
from collections import OrderedDict

from ROOT import TCanvas, TLatex, TFile, TMath, TH1F, TH2F, TLorentzVector
from ROOT import gStyle, gROOT, kTemperatureMap
from UserCode.HGCalMaskResolutionAna.RootTools import buildMedianProfile

"""
parser = argparse.ArgumentParser()
FLAGS, _ = Argparser.add_args(parser)
Argparser.print_args(FLAGS)
"""
parser = Argparser.Argparser()
FLAGS = parser.get_flags()
parser.print_args()

A=[TMath.Pi()*1.3**2, TMath.Pi()*2.6**2, TMath.Pi()*5.3**2]

def getEnergiesForCalibration(data,ireg):
    """
    Returns an array of reconstructed energies, eta and average noise
    """
    x=[]
    for i in xrange(0,data.GetEntriesFast()):
        data.GetEntry(i)
        for ia in xrange(1,2):
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
        xq=np.percentile(x, [i*100./nq for i in xrange(0,nq+1)], axis=0)

        #relative calibration versus eta        
        resolVsEta = TH2F('resolvseta',';|#eta|;#DeltaE/E', nq, 
                          Carray('d',xq[:,1]), 100,-1,1)
        for i in xrange(0,len(x)):
            genEn, genEta, recEn,_ = x[i]
            deltaE=recEn/genEn-1.
            resolVsEta.Fill(genEta,deltaE)
        calib['L0'][ireg]=calibrateSpectrum(resolVsEta,'SR%d'%ireg,
                                            FLAGS.plotLabel+' (PU=0)','pol2')
    
        #relative calibration versus energy
        resolVsEn = TH2F('resolvsen',';Reconstructed energy [GeV];#DeltaE/E', nq, 
                         Carray('d',xq[:,0]), 100,-1,1)
        for i in xrange(0,len(x)):
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
        xq=np.percentile(x, [i*100./nq for i in xrange(0,nq+1)], axis=0)

        #relative calibration versus eta
        resolVsNoiseFrac = TH2F('resolvsnoisefrac',';<Noise> [GeV];#DeltaE [GeV]', nq, 
                                Carray('d',xq[:,3]),50,-100,100)

        for i in xrange(0,len(x)):
            genEn,genEta,recEn,noiseEst=x[i]
            recEn=recEn/(calib['L0'][ireg].Eval(genEta)+1.0)
            recEn=recEn/(calib['L1'][ireg].Eval(recEn)+1.0)
            deltaE=recEn-genEn 
            resolVsNoiseFrac.Fill(noiseEst,deltaE)
    
        calib['L2'][ireg]=calibrateSpectrum(resolVsNoiseFrac,'SR%d'%ireg,
                                            FLAGS.plotLabel+' (PU=140)', 'pol2')
    return calib

def applyCalibrationTo(url,calib,title,ncands=2):
    """
    Applies the calibration to the photons and shows the energy and H->gg mass resolutions.
    """
    pfix=''.join(calib.keys())
    fIn=TFile.Open(url)
    data=fIn.Get('data')

    histos={}
    limsup, liminf = 0.4, -0.2
    for ireg in xrange(1,4):
        histos['dm%d'%ireg]  = TH1F('dm%d'%ireg, 
                                    ';#Delta m_{#gamma#gamma}/m_{#gamma#gamma};PDF',
                                    50, liminf, limsup)
        histos['den%d'%ireg] = TH1F('den%d'%ireg,';#Delta E/E;PDF',50,liminf,limsup)
    for h in histos:
        histos[h].Sumw2()
        histos[h].SetLineColor(1)
        histos[h].SetMarkerColor(1)
        histos[h].SetMarkerStyle(20)
        histos[h].SetDirectory(0)

    for i in xrange(0,data.GetEntriesFast()):
        data.GetEntry(i)

        #generator level photons
        genphotons=[]
        for ia in xrange(1,ncands+1):
            genen  = getattr(data,'genen%d'%ia)
            geneta = getattr(data,'geneta%d'%ia)
            genphi = getattr(data,'genphi%d'%ia)
            genphotons.append(TLorentzVector(0,0,0,0))
            genphotons[-1].SetPtEtaPhiM(genen/TMath.CosH(geneta),geneta,genphi,0.)
        genh = None
        if ncands==2 : genh=genphotons[0]+genphotons[1]

        #H->gg fiducial cuts
        if ncands==2:
            if genphotons[0].Pt()<20 or  genphotons[1].Pt()<20 : continue
            if genphotons[0].Pt()<40 and genphotons[1].Pt()<40 : continue
            if abs(genphotons[0].Eta())<1.5 or abs(genphotons[1].Eta())<1.5 : continue
            if abs(genphotons[0].Eta())>2.8 or abs(genphotons[1].Eta())>2.8 : continue

        #reconstructed photons in different regions
        for ireg in xrange(1,4):
            photons=[]
            for ia in xrange(1,ncands+1):
                genen    = getattr(data,'genen%d'%ia)
                geneta   = getattr(data,'geneta%d'%ia)
                genphi   = getattr(data,'genphi%d'%ia)
                recen    = getattr(data,'en%d_%d'%(ia,ireg))
                #avgnoise = getattr(data,'noise%d_%d'%(ia,ireg))
                avgnoise=getattr(data,'noise%d_3'%ia)*A[ireg-1]/A[2]

                if 'L0' in calib:
                    recen=recen/(calib['L0'][ireg].Eval(abs(geneta))+1.0)
                    if 'L1' in calib:
                        recen=recen/(calib['L1'][ireg].Eval(recen)+1.0)
                        if 'L2' in calib and ireg in calib['L2']:
                            recen=recen-calib['L2'][ireg].Eval(avgnoise)

                deltaE = recen/genen-1.
                histos['den%d'%ireg].Fill(deltaE)
                photons.append(TLorentzVector(0,0,0,0))
                photons[-1].SetPtEtaPhiM(recen/TMath.CosH(geneta),geneta,genphi,0.)

            if ncands<2: continue
            h = photons[0]+photons[1]
            deltaM=h.M()/genh.M()-1
            histos['dm%d'%ireg].Fill(deltaM)

    fIn.Close()

    c=TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.1)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.03)
    for ireg in xrange(1,4):
        for k in ['dm','den']:
            h=histos['%s%d'%(k,ireg)]
            if h.Integral()==0.: continue
            h.Scale(1./h.Integral())
            h.Draw()
            h.GetYaxis().SetTitleOffset(0.9)
            h.GetYaxis().SetRangeUser(0,h.GetMaximum()*1.2)
            h.Fit('gaus','M+')
            gaus=h.GetListOfFunctions().At(0)
            tex=TLatex()
            tex.SetTextFont(42)
            tex.SetTextSize(0.04)
            tex.SetNDC()
            tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')
            tex.DrawLatex(0.15,0.88,'SR%d (%s-calibrated)'%(ireg,pfix))
            tex.DrawLatex(0.15,0.84,'#mu=%3.3f#pm%3.3f'
                          %(gaus.GetParameter(1),gaus.GetParError(1)))
            tex.DrawLatex(0.15,0.80,'#sigma=%3.3f#pm%3.3f'
                          %(gaus.GetParameter(2),gaus.GetParError(2)))
            tex.SetTextAlign(31)
            tex.DrawLatex(0.97,0.96,title)
            c.SaveAs(os.path.join(FLAGS.outpath,pfix+h.GetName()+'.png'))

    #save in a local file
    fOut=TFile.Open('calib{}.root'.format(pfix),'RECREATE')
    for h in histos: histos[h].Write()
    fOut.Close()

def main():
    gStyle.SetOptStat(0)
    gStyle.SetOptTitle(0)
    gROOT.SetBatch(True)
    gStyle.SetPalette(kTemperatureMap)

    calib=OrderedDict()
    calib['L0']={}
    calib['L1']={}
    doL0L1Calibration(url=FLAGS.noPUFile, calib=calib)
    applyCalibrationTo(url=FLAGS.noPUFile, calib=calib, 
                       title=FLAGS.plotLabel+' (PU=0)', ncands=1)
    
    if len(sys.argv) > 5:
        puTag=sys.argv[5]
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
