import os
import pickle
from UserCode.HGCalMaskResolutionAna import Argparser
from collections import OrderedDict

from ROOT import gStyle, gROOT, kTemperatureMap
from UserCode.HGCalMaskResolutionAna.Calibration import Calibration

def main():
    calibration = Calibration(FLAGS.mingenen, FLAGS.mingeneta, FLAGS.maxgeneta,
                              FLAGS.plotLabel, FLAGS.outpath)
    calibration.L0L1Calibration(FLAGS.noPUFile)

    if FLAGS.PUFile != '':
        calibration.PUCalibration(FLAGS.PUFile)
        with open('calib_pu{}.pck'.format(FLAGS.puTag),'w') as cachefile:
            pickle.dump(calibration.calib, cachefile, pickle.HIGHEST_PROTOCOL)
    else:
        with open('calib_nopu.pck','w') as cachefile:
            pickle.dump(calibration.calib, cachefile, pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    gStyle.SetOptStat(0)
    gROOT.SetBatch(True)
    gStyle.SetPalette(kTemperatureMap)

    parser = Argparser.Argparser()
    FLAGS = parser.get_flags()
    parser.print_args()

    main()
