import os
import pickle
import numpy as np
from UserCode.HGCalMaskResolutionAna import Argparser
from collections import OrderedDict

from UserCode.HGCalMaskResolutionAna.Calibration import Calibration

def main():
    calibration = Calibration(FLAGS.mingenen, etaregions,
                              FLAGS.plotLabel, FLAGS.samples, FLAGS.mask, FLAGS.outpath)
    calibration.nopu_calibration(FLAGS.noPUFile, plot=False)

    if FLAGS.PUFile != '':
        str_pu = 'calib_'+FLAGS.samples+"_"+str(FLAGS.mask)+"_"+FLAGS.method+'_pu{}.pck'
        calibration.pu_calibration(FLAGS.PUFile, plot=False)
        calibration.save(str_pu.format(FLAGS.puTag))
    else:
        str_nopu = 'calib_'+FLAGS.samples+"_"+str(FLAGS.mask)+"_"+FLAGS.method+'_nopu.pck'
        calibration.save(str_nopu)


if __name__ == "__main__":
    from ROOT import gStyle, gROOT, kTemperatureMap
    gStyle.SetOptStat(0)
    gROOT.SetBatch(True)
    gStyle.SetPalette(kTemperatureMap)

    parser = Argparser.Argparser()
    FLAGS = parser.get_flags()
    parser.print_args()

    if FLAGS.method == 'fineeta':
        if FLAGS.samples == 'inner':
            etaregions = np.linspace(2.7, 3.03, 331)
        elif FLAGS.samples == 'outer':
            etaregions = np.linspace(1.45, 1.65, 201)
    elif FLAGS.method == 'ed':
        if FLAGS.samples == 'inner':
            etaregions = np.array((2.7, 2.94))
        elif FLAGS.samples == 'outer':
            etaregions = np.array((1.55, 1.65))
    main()
