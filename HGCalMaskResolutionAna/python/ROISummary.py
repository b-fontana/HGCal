import numpy as np
from ROOT import TLorentzVector, TMath

class ROISummary():
    """
    Summary of the energy collected in a region of interest (ROI). The ROIs correspond
    either to a particle signal (photons/pions) or regions to study the noise close to the signal ROI. Each ROI may thus have many hits in different layers.
    The energy of the hits is stored per layer and per signal region.
    The signal region numering scheme starts at 1 for the smallest one.
    """
    def __init__(self, genP4, Nregions=3, Nlayers=28):
        self.genP4   = genP4
        self.Nregions_ = Nregions
        self.Nlayers_ = Nlayers
        self.recEn = np.zeros((self.Nregions_,self.Nlayers_))
        self.noiseEn = np.zeros((self.Nregions_,self.Nlayers_))
        self.layeroffset_ = 28
        self.signalHits = [[] for _ in range(self.Nregions_)]
        self.mipthresh_ = [4, 5, 6]

    def AddHit(self, en, layer, isNoise, regIdx):
        if layer > self.Nlayers_:
            layer -= self.layeroffset_
        if isNoise:
            self.noiseEn[regIdx-1,layer-1] += en
        else:
            self.recEn[regIdx-1,layer-1] += en
            self.signalHits[regIdx-1].append(en)

    def Cglob(self, regIdx):
        """Factor useful por pion compensation. 
        The energy of the Hits added must be in MIPs.
        It calculates the average using the hits that were already added. """
        cglob = np.empty(3)
        for ith,th in enumerate(self.mipthresh_):
            under_thresh = [it for ireg in range(regIdx) for it in self.signalHits[ireg] if it<th]
            try:
                average = self.AverageHitEn(regIdx)
                under_average = [it for ireg in range(regIdx) for it in self.signalHits[ireg] if it<average]
                cglob[ith] = float(len(under_thresh))/len(under_average)
            except ZeroDivisionError:
                if len(under_thresh)<2: 
                    print("Warning: cglob set to zero.")
                    cglob[ith] = 0
                else:                    
                    print(len(under_thresh))
                    print(len(under_average))
                    print(self.RecoP4(regIdx).E())
                    print(self.AverageHitEn(ireg)) 
                    raise

        return cglob

    def AverageHitEn(self, regIdx):
        """It calculates the average using the hits that were already added."""
        tmp = [it for ireg in range(regIdx) for it in self.signalHits[ireg]]
        tot_en = sum(tmp)
        try:
            assert(np.isclose(self.RecoP4(regIdx).E(),tot_en))
        except AssertionError:
            print(self.RecoP4(regIdx).E())
            print(tot_en)
            raise
        try:
            av = tot_en/len(tmp)
        except ZeroDivisionError:
            print("The average could not be calculated.")
            raise
        return av

    def RecoP4(self, regIdx):
        """
        Sums the energy collected in a given signal region
        and returns the 4-vector.
        """
        p4 = TLorentzVector(0,0,0,0)
        en = np.sum( np.sum(self.recEn, axis=1)[:regIdx], axis=0)
        pt = en / TMath.CosH(self.genP4.Eta())
        p4.SetPtEtaPhiM(pt, self.genP4.Eta(), self.genP4.Phi(), 0)
        return p4

    def RecoEnergyDeposited(self, regIdx, layer):
        #return sum([self.recEn[x,layer-1] for x in range(regIdx)])
        return np.sum(self.recEn, axis=0)[layer-1]

    def FractionEnergyDeposited(self, regIdx, layer):
        try:
            return self.RecoEnergyDeposited(regIdx, layer) / self.RecoP4(regIdx).E()
        except ZeroDivisionError:
            return 0.

    def NoiseInROI(self, regIdx):
        """
        Gets the averaged noise in a region.
        """
        noise = 0
        for ireg in range(regIdx):
            for il in range(self.Nlayers_):
                noise += self.noiseEn[ireg][il] 
        return noise / 5.
        
    def NoiseInLayer(self, regIdx, layer):
        return sum([self.noiseEn[x][layer-1] for x in range(regIdx)])
