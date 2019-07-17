import os
import pickle
from UserCode.HGCalMaskResolutionAna import Argparser
from collections import OrderedDict

from ROOT import gStyle, gROOT, kTemperatureMap
from UserCode.HGCalMaskResolutionAna.Calibration import Calibration

def main():
    calibration = Calibration(FLAGS.mingenen, etaregions,
                              FLAGS.plotLabel, FLAGS.samples, FLAGS.mask, FLAGS.outpath)
    calibration.L0L1Calibration(FLAGS.noPUFile)

    if FLAGS.PUFile != '':
        str_pu = 'calib_'+FLAGS.samples+"_"+str(FLAGS.mask)+"_"+FLAGS.method+'_pu{}.pck'
        calibration.PUCalibration(FLAGS.PUFile)
        with open(calib_str_pu.format(FLAGS.puTag),'w') as cachefile:
            pickle.dump(calibration.calib, cachefile, pickle.HIGHEST_PROTOCOL)
    else:
        str_nopu = 'calib_'+FLAGS.samples+"_"+str(FLAGS.mask)+"_"+FLAGS.method+'_nopu.pck'
        with open(calib_str_nopu,'w') as cachefile:
            pickle.dump(calibration.calib, cachefile, pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    gStyle.SetOptStat(0)
    gROOT.SetBatch(True)
    gStyle.SetPalette(kTemperatureMap)

    parser = Argparser.Argparser()
    FLAGS = parser.get_flags()
    parser.print_args()

    if FLAGS.method == 'fineeta':
        etaregions = np.round(np.arange(2.7,3.031,0.001).tolist(), 3).tolist()
    elif FLAGS.method == 'ed':
        etaregions = [2.7, 2.94]

    main()
