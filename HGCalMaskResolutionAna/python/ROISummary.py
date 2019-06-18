from ROOT import TLorentzVector, TMath

class ROISummary():
    """
    Summary of the energy collected in a region of interest.
    """
    def __init__(self, genP4):
        self.genP4_   = genP4
        self.recEn_   = [0., 0., 0.]
        self.noiseEn_ = [0., 0., 0.]
        
    def getGenP4():
        return self.genP4_

    def addHit(self, en, isNoise, regIdx):
        if isNoise:
            self.noiseEn_[regIdx-1] += en
        else:
            self.recEn_[regIdx-1] += en

    def getReconstructedP4(self, regIdx):
        """
        Sum the energy collected in a given signal region and return the 4-vector.
        """
        p4 = TLorentzVector(0,0,0,0)
        totalEn = sum(self.recEn_[0:regIdx]) #sum up necessary rings
        pt = totalEn / TMath.CosH(self.genP4_.Eta())
        p4.SetPtEtaPhiM(pt, self.genP4_.Eta(), self.genP4_.Phi(), 0)
        return p4

    def getNoiseInRing(self,regIdx):
        """
        Gets the noise in a region.
        """
        return sum(self.noiseEn_[0:regIdx])
