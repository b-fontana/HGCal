import os, sys
import numpy as np

from ROOT import TFile, TGraph
from ROOT import Double
from UserCode.HGCalMaskVisualProd.RootUtils import PyDoubleBufferToList as toList

class PartialWafersStudies(object):
    """
    Project base class.
    """
    def __init__(self):
        self.nlayers = 28
        self.nreg = 3

    @property
    def nlayers(self):
        return self.__nlayers
    @nlayers.setter
    def nlayers(self, nlayers):
        self.__nlayers = nlayers
