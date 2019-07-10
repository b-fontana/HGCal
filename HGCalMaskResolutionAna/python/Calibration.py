import os
import numpy as np
from array import array as Carray
from collections import OrderedDict

from ROOT import TCanvas, TFile, TH2F, TLatex
from UserCode.HGCalMaskResolutionAna.PartialWafersStudies import PartialWafersStudies
from UserCode.HGCalMaskResolutionAna.RootTools import buildMedianProfile
from UserCode.HGCalMaskVisualProd.SystemUtils import EtaStr as ES

class Calibration(PartialWafersStudies, object):
    def __init__(self, mingenen, eta_regions, 
                 label, samples, mask, outpath):
        super(Calibration, self).__init__()
        self.calib = OrderedDict({'L0': {}, 'L1': {}, 'L2': {}})
        self.mingenen = mingenen
        self.etas_l = eta_regions[:-1]
        self.etas_r = np.roll(eta_regions, shift=-1)[:-1]
        self.label = label
        self.samples = samples
        self.mask = str(mask)
        self.outpath = outpath

    def _CalibrateSpectrum(self, h, title, proc, func='pol1', plot=True):
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
            c.SaveAs(os.path.join(self.outpath,h.GetName()+'_'+title+'.png'))
            prof.Delete()        
            c.Delete()

        return calibGr

    def _EnergiesForCalibration(self, data, ireg):
        """
        Returns an array of reconstructed energies, eta and average noise
        """
        x = [[] for _ in range(len(self.etas_l))]
        for i in range(0, data.GetEntriesFast()):
            data.GetEntry(i)
            genen = getattr(data,'genen')
            if genen < self.mingenen: 
                continue
            geneta = abs(getattr(data,'geneta'))
            region_check = False
            for ieta1,ieta2 in zip(self.etas_l, self.etas_r):
                if (geneta < ieta1 or geneta > ieta2): 
                    continue
                if region_check:
                    raise ValueError('This event was already assigned to an eta region.')
                region_check = True
                recen = getattr(data,'en_sr{}_ROI'.format(ireg))
                avgnoise = ( getattr(data,'noise_sr3_ROI')
                             *self.sr_area[ireg-1]/self.sr_area[2] )
                x[self.etas_l.index(ieta1)].append( [genen, geneta, recen, avgnoise] )
        return [np.array(x[i]) for i in range(len(x))]


    def L0L1Calibration(self, f, nq=6, plot=True):
        """
        Performs the L0 (uniform eta response) and L1 (absolute scale) calibrations.
        """
        fIn = TFile.Open(f)
        data = fIn.Get('data')

        for ireg in [1,2,3]:
            x = self._EnergiesForCalibration(data, ireg)
            for ieta1,ieta2 in zip(self.etas_l, self.etas_r):
                idx = self.etas_l.index(ieta1)
                idstr = 'sr{}_from{}to{}'.format(ireg, ES(ieta1), ES(ieta2))
                xq = np.percentile(x[idx], [i*100./nq for i in range(0,nq+1)], axis=0)
            
                #relative calibration versus eta
                hn = idstr+'_'+self.samples+'_mask'+self.mask
                htmp = TH2F('resVSeta_'+hn, 
                            ';|#eta|;#DeltaE/E', nq, Carray('d',xq[:,1]), 100,-1,1)
                for i in range(0,len(x[idx])):
                    genEn, genEta, recEn,_ = x[idx][i]
                    deltaE = (recEn/genEn) - 1.
                    htmp.Fill(genEta, deltaE)
                self.calib['L0'][idstr] = self._CalibrateSpectrum(htmp,'SR%d'%ireg,
                                                                  self.label+' (PU=0)',
                                                                  func='pol2',
                                                                  plot=plot)
                htmp.Delete()

                #relative calibration versus energy
                htmp = TH2F('resVSen_'+hn,
                            ';Reconstructed energy [GeV];#DeltaE/E', nq, 
                            Carray('d',xq[:,0]), 100,-1,1)
                for i in range(0,len(x[idx])):
                    genen, geneta, recen, _ = x[idx][i]
                    recen /= (self.calib['L0'][idstr].Eval(geneta)+1.0)
                    deltaE = (recen/genen) - 1.
                    htmp.Fill(recen, deltaE)
                self.calib['L1'][idstr] = self._CalibrateSpectrum(htmp,'SR%d'%ireg,
                                                                  self.label+' (PU=0)',
                                                                  func='pol1',
                                                                  plot=plot)
                htmp.Delete()
        fIn.Close()

    

    def PUCalibration(f, eta_region, nq=10):
        """
        Parametrizes the absolute shift in energy as function of the 
        average noise in the SR cone.
        """
        fIn=TFile.Open(f)
        data=fIn.Get('data')

        for ireg in [1,2,3]:
            idstr = ( 'sr{}_from{}to{}'.format(ireg, eta_region[0], eta_region[1])
                      if eta_region!=(0,0) else 'sr{}'.format(ireg) )
            x = self._EnergiesForCalibration(data, ireg, eta_region)
            xq=np.percentile(x, [i*100./nq for i in range(0,nq+1)], axis=0)

            #relative calibration versus eta
            resolVsNoiseFrac = TH2F('resVSnoisefrac_'+self.samples+'_'+self.mask+'_',
                                    ';<Noise> [GeV];#DeltaE [GeV]', nq, 
                                    Carray('d',xq[:,3]),50,-100,100)

            for i in range(0,len(x)):
                genEn,genEta,recEn,noiseEst=x[i]
                recEn /= self.calib['L0'][idstr].Eval(geneta)+1.0
                recEn /= self.calib['L1'][idstr].Eval(recen)+1.0
                deltaE = recen - genen 
                resolVsNoiseFrac.Fill(noiseEst,deltaE)
    
            self.calib['L2'][ireg] = self._CalibrateSpectrum(resolVsNoiseFrac,'SR%d'%ireg,
                                                    self.label+' (PU=140)', func='pol2')
        fIn.Close()
