from ROOT import TLorentzVector, TMath

class ROISummary():
    """
    Summary of the energy collected in a region of interest (ROI). The ROIs correspond
    either a photon or regions to study the noise close to the signal ROI. Each ROI
    may thus have many hits in different layers.
    The energy of the hits is stored per layer.
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

    def getRecoP4(self, regIdx):
        """
        Sums the energy collected in a given signal region
        and returns the 4-vector.
        """
        p4 = TLorentzVector(0,0,0,0)
        en = 0
        for ireg in range(regIdx):
            for il in range(self.Nlayers_):
                en += self.recEn_[ireg][il] 
        pt = en / TMath.CosH(self.genP4_.Eta())
        p4.SetPtEtaPhiM(pt, self.genP4_.Eta(), self.genP4_.Phi(), 0)
        return p4

    def getRecoEnergyDeposited(self, regIdx, layer):
        return sum([self.recEn_[x][layer-1] for x in range(regIdx)])

    def getNoiseInROI(self, regIdx):
        """
        Gets the averaged noise in a region.
        """
        noise = 0
        for ireg in range(regIdx):
            for il in range(self.Nlayers_):
                noise += self.noiseEn_[ireg][il] 
        return noise / 5.
        
    def getNoiseInLayer(self, regIdx, layer):
        return sum([self.noiseEn_[x][layer-1] for x in range(regIdx)])
