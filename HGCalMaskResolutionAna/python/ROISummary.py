from ROOT import TLorentzVector, TMath

class ROISummary():
    """
    Summary of the energy collected in a region of interest (ROI). The ROIs correspond
    either a photon or regions to study the noise close to the signal ROI. Each ROI
    may thus have many hits in different layers.
    The energy of the hits is stored per layer and per signal region.
    """
    def __init__(self, genP4, Nregions=3, Nlayers=28):
        self.genP4   = genP4
        self.Nregions_ = Nregions
        self.Nlayers_ = Nlayers
        self.recEn = [[0. for _ in range(self.Nlayers_)] for _ in range(self.Nregions_)]
        self.noiseEn = [[0. for _ in range(self.Nlayers_)] for _ in range(self.Nregions_)]

    def AddHit(self, en, layer, isNoise, regIdx):
        if isNoise:
            self.noiseEn[regIdx-1][layer-1] += en
        else:
            self.recEn[regIdx-1][layer-1] += en

    def RecoP4(self, regIdx):
        """
        Sums the energy collected in a given signal region
        and returns the 4-vector.
        """
        p4 = TLorentzVector(0,0,0,0)
        en = 0
        for ireg in range(regIdx):
            for il in range(self.Nlayers_):
                en += self.recEn[ireg][il] 
        pt = en / TMath.CosH(self.genP4.Eta())
        p4.SetPtEtaPhiM(pt, self.genP4.Eta(), self.genP4.Phi(), 0)
        return p4

    def RecoEnergyDeposited(self, regIdx, layer):
        return sum([self.recEn[x][layer-1] for x in range(regIdx)])

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
