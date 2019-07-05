import os
import numpy as np
from array import array as Carray
from collections import OrderedDict

from ROOT import TCanvas, TFile, TH2F, TLatex
from UserCode.HGCalMaskResolutionAna.PartialWafersStudies import PartialWafersStudies
from UserCode.HGCalMaskResolutionAna.RootTools import buildMedianProfile

class Calibration(PartialWafersStudies, object):
    def __init__(self, mingenen, mingeneta, maxgeneta, 
                 plotLabel, outpath):
        super(Calibration, self).__init__()
        self.calib = OrderedDict({'L0': {}, 'L1': {}, 'L2': {}})
        self.mingenen = mingenen
        self.maxgeneta = maxgeneta
        self.mingeneta = mingeneta
        self.plotLabel = plotLabel
        self.outpath = outpath

    def _calibrateSpectrum(self, h, title, proc, func='pol1', plot):
        """
        Calibrates a DeltaE/E versus x spectrum based on the profile.
        """
        prof = buildMedianProfile(h)
        prof.Fit(func)
        calibGr = prof.GetListOfFunctions().At(0).Clone(h.GetName()+title+'_calib')

        if plot:
            c = TCanvas('c','c',500,500)
            c.SetTopMargin(0.05)
            c.SetBottomMargin(0.1)
            c.SetLeftMargin(0.12)
            c.SetRightMargin(0.1)
            h.Draw('colz')
            prof.Draw('e2p')

            tex = TLatex()
            tex.SetTextFont(42)
            tex.SetTextSize(0.04)
            tex.SetNDC()
            tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')
            tex.DrawLatex(0.15,0.88,title)
            tex.SetTextAlign(31)
            tex.DrawLatex(0.97,0.96,proc)        
            c.SaveAs(os.path.join(self.outpath,h.GetName()+title+'.png'))
            prof.Delete()        
            c.Delete()

        return calibGr

    def _getEnergiesForCalibration(self, data, ireg):
        """
        Returns an array of reconstructed energies, eta and average noise
        """
        x = []
        for i in range(0, data.GetEntriesFast()):
            data.GetEntry(i)
            genen = getattr(data,'genen')
            geneta = abs(getattr(data,'geneta'))
            if genen < self.mingenen: continue
            if geneta < self.mingeneta: continue
            if geneta > self.maxgeneta: continue
            for il in range(1,self.nsr+1):
                recen = getattr(data,'en_sr{}_ROI'.format(ireg,il))
                avgnoise = ( getattr(data,'noise_sr3_layer{}'.format(il))
                             *self.sr_area[ireg-1]/self.sr_area[2] )
                x.append( [genen, geneta, recen, avgnoise] )
        return np.array(x)


    def L0L1Calibration(self, f, nq=6, plot=True):
        """
        Performs the L0 (uniform eta response) and L1 (absolute scale) calibrations.
        """
        fIn = TFile.Open(f)
        data = fIn.Get('data')
        for ireg in [1,2,3]:
            x = self._getEnergiesForCalibration(data, ireg)
            xq = np.percentile(x, [i*100./nq for i in range(0,nq+1)], axis=0)

            #relative calibration versus eta        
            resVsEta = TH2F('resvseta',';|#eta|;#DeltaE/E', nq, 
                              Carray('d',xq[:,1]), 100,-1,1)
            for i in range(0,len(x)):
                genEn, genEta, recEn,_ = x[i]
                deltaE=recEn/genEn-1.
                resVsEta.Fill(genEta,deltaE)
            self.calib['L0'][ireg] = self._calibrateSpectrum(resVsEta,'SR%d'%ireg,
                                                        self.plotLabel+' (PU=0)','pol2',
                                                             plot=plot)
    
            #relative calibration versus energy
            resVsEn = TH2F('resvsen',';Reconstructed energy [GeV];#DeltaE/E', nq, 
                             Carray('d',xq[:,0]), 100,-1,1)
            for i in range(0,len(x)):
                genen, geneta, recen, _ = x[i]
                recen /= (self.calib['L0'][ireg].Eval(geneta)+1.0)
                deltaE = (recen/genen) - 1.
                resVsEn.Fill(recen, deltaE)
                resVsEta.Fill(genEta, deltaE)
            self.calib['L1'][ireg] = self._calibrateSpectrum(resVsEn,'SR%d'%ireg,
                                                        self.plotLabel+' (PU=0)','pol1',
                                                             plot=plot)
        fIn.Close()

    def PUCalibration(f, nq=10):
        """
        Parametrizes the absolute shift in energy as function of the 
        average noise in the SR cone.
        """
        fIn=TFile.Open(f)
        data=fIn.Get('data')

        for ireg in [1,2,3]:
            x = self._getEnergiesForCalibration(data,ireg)
            xq=np.percentile(x, [i*100./nq for i in range(0,nq+1)], axis=0)

            #relative calibration versus eta
            resolVsNoiseFrac = TH2F('resolvsnoisefrac',';<Noise> [GeV];#DeltaE [GeV]', nq, 
                                    Carray('d',xq[:,3]),50,-100,100)

            for i in range(0,len(x)):
                genEn,genEta,recEn,noiseEst=x[i]
                recEn /= self.calib['L0'][ireg].Eval(geneta)+1.0
                recEn /= self.calib['L1'][ireg].Eval(recen)+1.0
                deltaE = recen - genen 
                resolVsNoiseFrac.Fill(noiseEst,deltaE)
    
            self.calib['L2'][ireg] = self._calibrateSpectrum(resolVsNoiseFrac,'SR%d'%ireg,
                                                    FLAGS.plotLabel+' (PU=140)', 'pol2')
        fIn.Close()
