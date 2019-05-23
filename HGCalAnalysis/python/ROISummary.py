import ROOT

class ROISummary():

    """summary of the energy collected in a region of interest"""

    def __init__(self,genP4):

        """start with MC truth and init energy arrays"""

        self.genP4   = genP4
        self.recEn   = [0., 0., 0.]
        self.noiseEn = [0., 0., 0.]

    def addHit(self,en,isNoise,regIdx):
        
        """add hit"""

        if isNoise:
            self.noiseEn[regIdx-1]+=en
        else:
            self.recEn[regIdx-1]+=en

    def getReconstructedP4(self,regIdx):
        
        """sum up the energy collected in a given signal region and return the 4-vector"""
        
        p4=ROOT.TLorentzVector(0,0,0,0)
        totalEn=sum(self.recEn[0:regIdx]) #sum up necessary rings
        pt=totalEn/ROOT.TMath.CosH(self.genP4.Eta())
        p4.SetPtEtaPhiM(pt,self.genP4.Eta(),self.genP4.Phi(),0)
        return p4

    def getNoiseInRing(self,regIdx):

        """gets noise in a region"""

        return sum(self.noiseEn[0:regIdx])
