from ROOT import TLorentzVector, TMath

class ROISummary():
    """
    Summary of the energy collected in a region of interest.
    """
    def __init__(self, genP4, Nregions=3, Nlayers=28):
        self.genP4_   = genP4
        self.Nregions_ = Nregions
        self.Nlayers_ = Nlayers
        self.recEn_ = [[0. for _ in range(self.Nlayers_)] for _ in range(self.Nregions_)]
        self.noiseEn_ = [[0. for _ in range(self.Nlayers_)] for _ in range(self.Nregions_)]

    def getGenP4(self):
        return self.genP4_

    def addHit(self, en, layer, isNoise, regIdx):
        if isNoise:
            self.noiseEn_[regIdx-1][layer-1] += en
        else:
            self.recEn_[regIdx-1][layer-1] += en

    def getReconstructedP4(self, regIdx, layer):
        """
        Sums the energy collected in a given signal region and layer
        and returns the 4-vector.
        """
        p4 = TLorentzVector(0,0,0,0)
        en = sum([self.recEn_[ireg][layer-1] for ireg in range(regIdx)])
        pt = en / TMath.CosH(self.genP4_.Eta())
        p4.SetPtEtaPhiM(pt, self.genP4_.Eta(), self.genP4_.Phi(), 0)
        return p4

    def getTotalReconstructedP4(self, regIdx):
        """
        Sums the energy collected in a given signal region and returns the 4-vector.
        """
        p4 = TLorentzVector(0,0,0,0)
        en = 0
        for il in range(self.Nlayers):
            en += sum([self.recEn_[ireg][il] for ireg in range(regIdx)])
        pt = totalEn / TMath.CosH(self.genP4_.Eta())
        p4.SetPtEtaPhiM(pt, self.genP4_.Eta(), self.genP4_.Phi(), 0)
        return p4

    def getNoiseInRegion(self, regIdx, layer):
        """
        Gets the noise in a region.
        """
        return sum([self.noiseEn_[ireg][layer-1] for ireg in range(regIdx)])/5.

    def getTotalNoiseInRegion(self, regIdx):
        """
        Gets the noise in a region.
        """
        for il in range(self.Nlayers):
            en += sum([self.noiseEn_[ireg][il] for ireg in range(regIdx)])
        return en/5.
