import os
import numpy as np

class PartialWafersStudies(object):
    """
    Project base class.
    """
    def __init__(self):
        self._nlayers = 28
        self.nsr = 3
        self.sr_dist = (1.3, 2.6, 5.3)
        self.sr_area = (np.pi*self.sr_dist[0]**2,
                        np.pi*self.sr_dist[1]**2,
                        np.pi*self.sr_dist[2]**2)

    @property
    def nlayers(self):
        return self._nlayers
    @nlayers.setter
    def nlayers(self, nlayers):
        self._nlayers = nlayers
