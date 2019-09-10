import numpy as np
from ROOT import TLorentzVector, TMath

class ROISummary():
    """
    Summary of the energy collected in a region of interest (ROI). The ROIs correspond
    either to a particle signal (photons/pions) or regions to study the noise close to the signal ROI. Each ROI may thus have many hits in different layers.
    The energy of the hits is stored per layer and per signal region.
    The signal region numbering scheme starts at 1 for the smallest one. The subdetector numbering scheme is: 1->CEE, 2->HEF, 3->HEB. 
    """
    def __init__(self, genP4, Nregions=3, Nlayers=28, Nsubdets=3, PartType='Photon'):
        if PartType not in ('Photon', 'Pion'):
            raise ValueError('Only Photons and Pions are currently supported.')
            exit()

        self.genP4   = genP4
        self.Nregions_ = Nregions
        self.Nlayers_ = Nlayers
        self.Nsubdets_ = Nsubdets
        self.PartType = PartType
        #Nlayers should be 22 or 14 for the hadronic part, but it is easier to work with rectangular arrays
        self.recEn = np.zeros((self.Nregions_, self.Nlayers_, self.Nsubdets_))
        self.noiseEn = np.zeros((self.Nregions_, self.Nlayers_, self.Nsubdets_))
        self.layeroffset_ = 28
        self.mipthresh_ = [4, 5, 6]

    def AddHit(self, en, layer, subdet, isNoise, regIdx):
        if layer > self.Nlayers_:
            layer -= self.layeroffset_
        if isNoise:
            self.noiseEn[regIdx-1,layer-1,subdet-1] += en
        else:
            self.recEn[regIdx-1,layer-1,subdet-1] += en
            #self.signalHits[regIdx-1].append(en)

    def RecoP4(self, regIdx):
        """
        Sums the energy collected in a given signal region
        and returns the 4-vector.
        """
        p4 = TLorentzVector(0,0,0,0)
        en = np.sum(self.recEn[:regIdx], axis=None)
        pt = en / TMath.CosH(self.genP4.Eta())
        p4.SetPtEtaPhiM(pt, self.genP4.Eta(), self.genP4.Phi(), 0)
        return p4

    def RecoEnergyDeposited(self, regIdx, layer):
        return np.sum(self.recEn[:regIdx], axis=(2,0))[layer-1]

    def FractionEnergyDeposited(self, regIdx, layer):
        try:
            return self.RecoEnergyDeposited(regIdx, layer) / self.RecoP4(regIdx).E()
        except ZeroDivisionError:
            return 0.

    def SubdetEnergyDeposited(self, regIdx, det):
        return np.sum(self.recEn[:regIdx,:,det-1], axis=(1,0))

    def SubdetNoiseDeposited(self, regIdx, det):
        return np.sum(self.noiseEn[:regIdx,:,det-1], axis=(1,0)) / (5. if self.PartType=='Photon' else 2.)

    def NoiseInROI(self, regIdx):
        """Gets the averaged noise in a region."""
        noise = 0
        for ireg in range(regIdx):
            for il in range(self.Nlayers_):
                for idet in range(self.Nsubdets_):
                    noise += self.noiseEn[ireg,il,idet] 
        return noise / (5. if self.PartType=='Photon' else 2.)
        
    def NoiseInLayer(self, regIdx, layer):
        return np.sum(self.noiseEn[:regIdx], axis=(2,0))[layer-1]
