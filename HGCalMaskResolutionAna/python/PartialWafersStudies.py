import numpy as np

class PartialWafersStudies(object):
    """
    Project base class.
    """
    def __init__(self, flags=None):
        self._nlayers = 28
        self.nsr = 3
        self.sr_dist = (1.3, 2.6, 5.3)
        self.sr_area = (np.pi*self.sr_dist[0]**2,
                        np.pi*self.sr_dist[1]**2,
                        np.pi*self.sr_dist[2]**2)
        self._etaregions = list()
        self.flags = flags
        self._paths = None

    @property
    def nlayers(self):
        return self._nlayers
    @nlayers.setter
    def nlayers(self, nlayers):
        self._nlayers = nlayers

    @property
    def etaregions(self):
        if self.flags.method == 'fineeta':
            if self.flags.samples == 'inner':
                self._etaregions = np.linspace(2.7, 3.03, 331)
            elif self.flags.samples == 'outer':
                self._etaregions = np.linspace(1.45, 1.65, 201)
        elif self.flags.method == 'ed':
            if self.flags.samples == 'inner':
                self._etaregions = np.array((2.7, 2.94))
            elif self.flags.samples == 'outer':
                self._etaregions = np.array((1.55, 1.65))
        else:
            raise ValueError(" You cannot retrieve 'etaregions' "
                             " without specifying 'method'." )
        return self._etaregions

    @property
    def paths(self):
        from collections import namedtuple
        p = namedtuple('paths', 'weights plots calibrations_pu calibrations_nopu')
        return p('weights/calibshowers_mask' + str(self.flags.mask) + '_' + 
                 self.flags.samples+'_mode2_' + self.flags.method + '.root', 
                 'allplots_' + self.flags.samples + '_' + str(self.flags.mask)
                 + '_' + self.flags.method+'.root',
                 'calibrations/calib_' + self.flags.samples + '_' + str(self.flags.mask) 
                 + "_" + self.flags.method + '_pu{}.pck'.format(self.flags.puTag), 
                 'calibrations/calib_' + self.flags.samples + str(self.flags.mask)
                 + "_" + self.flags.method + '_nopu.pck')
