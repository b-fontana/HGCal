from UserCode.HGCalMaskResolutionAna.PartialWafersStudies import PartialWafersStudies

class Calibration(PartialWafersStudies, object):
    def __init__(self, flags):
        import numpy as np
        from collections import OrderedDict

        super(Calibration, self).__init__(flags)
        self.calib = OrderedDict({'L0': {}, 'L1': {}, 'L2': {}})
        self.mingenen = self.flags.mingenen
        self.etas_l = self.etaregions[:-1]
        self.etas_r = np.roll(self.etaregions, shift=-1)[:-1]
        self.label = self.flags.plotLabel
        self.samples = self.flags.samples
        self.mask = str(self.flags.mask)
        self.outpath = self.flags.outpath

    def _calibrate_spectrum(self, h, title, proc, func='pol1', plot=True):
        """
        Calibrates a DeltaE/E versus x spectrum based on the profile.
        """
        import os
        from UserCode.HGCalMaskResolutionAna.RootTools import buildMedianProfile

        prof = buildMedianProfile(h)
        prof.Fit(func)
        calibGr = prof.GetListOfFunctions().At(0).Clone(h.GetName()+title+'_calib')

        if plot:
            from ROOT import TCanvas, TLatex
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

    def _energies_for_calibration(self, data, ireg):
        """
        Returns an array of reconstructed energies, eta and average noise
        """
        import numpy as np
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
                npidx = np.where(self.etas_l==ieta1)[0][0]
                x[npidx].append( [genen, geneta, recen, avgnoise] )
        return [np.array(x[i]) for i in range(len(x))]


    def nopu_calibration(self, nq=6, plot=True):
        """
        Performs the L0 (uniform eta response) and L1 (absolute scale) calibrations.
        When dealing with pile-up, 'pu_calibration()' has to be called afterwords.
        """
        import numpy as np
        from array import array as Carray
        from ROOT import TFile, TH2F
        from UserCode.HGCalMaskVisualProd.SystemUtils import EtaStr as ES

        fIn = TFile.Open(self.flags.noPUFile)
        data = fIn.Get('data')

        for ireg in [1,2,3]:
            x = self._energies_for_calibration(data, ireg)
            for ieta1,ieta2 in zip(self.etas_l, self.etas_r):
                idx = np.where(self.etas_l==ieta1)[0][0]
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
                self.calib['L0'][idstr] = self._calibrate_spectrum(htmp,'SR%d'%ireg,
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
                self.calib['L1'][idstr] = self._calibrate_spectrum(htmp,'SR%d'%ireg,
                                                                   self.label+' (PU=0)',
                                                                   func='pol1',
                                                                   plot=plot)
                htmp.Delete()
        fIn.Close()

    #correct arguments of np.percentile
    def pu_calibration(self, nq=10, plot=True):
        """
        To be used in pile-up situations only.
        It parametrizes the absolute shift in energy as function of the 
        average noise in the SR cone.
        """
        import numpy as np
        from array import array as Carray
        from ROOT import TFile, TH2F
        from UserCode.HGCalMaskVisualProd.SystemUtils import EtaStr as ES

        fIn=TFile.Open(self.flags.PUFile)
        data=fIn.Get('data')

        for ireg in [1,2,3]:
            x = self._energies_for_calibration(data, ireg)
            for ieta1,ieta2 in zip(self.etas_l, self.etas_r):
                idx = np.where(self.etas_l==ieta1)[0][0]
                idstr = 'sr{}_from{}to{}'.format(ireg, ES(ieta1), ES(ieta2))
                xq=np.percentile(x, [i*100./nq for i in range(0,nq+1)], axis=0)
                print(xq)
                #relative calibration versus eta
                resolVsNoiseFrac = TH2F('resVSnoisefrac_'+self.samples+'_'+self.mask+'_',
                                        ';<Noise> [GeV];#DeltaE [GeV]', nq, 
                                        Carray('d',xq[:,3]),50,-100,100)
                for i in range(0,len(x)):
                    genen, geneta, recen, noiseest = x[i]
                    recen /= (self.calib['L0'][idstr].Eval(geneta)+1.0)
                    recen /= (self.calib['L1'][idstr].Eval(recen)+1.0)
                    deltaE = recen - genen 
                    resolVsNoiseFrac.Fill(noiseEst, deltaE)
    
                self.calib['L2'][ireg] = self._calibrate_spectrum(resolVsNoiseFrac,
                                                                  'SR%d'%ireg,
                                                                  self.label+' (PU=140)', 
                                                                  func='pol2',
                                                                  plot=plot)
        fIn.Close()

    def save(self, name):
        if name[-4:] != '.pck':
            raise ValueError('The calibration has to be saved as a pickle file.')
        import pickle
        with open(name,'w') as cachefile:
            pickle.dump(self.calib, cachefile, pickle.HIGHEST_PROTOCOL)
