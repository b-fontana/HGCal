from UserCode.HGCalMaskResolutionAna import Argparser
from UserCode.HGCalMaskResolutionAna.Calibration import Calibration

def main():
    calibration = Calibration(FLAGS)
    calibration.nopu_calibration(plot=False)

    if FLAGS.PUFile != '':
        calibration.pu_calibration(plot=False)
        calibration.save( base.paths.calibrations_pu )
    else:
        calibration.save( base.paths.calibrations_nopu )


if __name__ == "__main__":
    from ROOT import gStyle, gROOT, kTemperatureMap
    gStyle.SetOptStat(0)
    gROOT.SetBatch(True)
    gStyle.SetPalette(kTemperatureMap)

    parser = Argparser.Argparser()
    FLAGS = parser.get_flags()
    parser.print_args()
    base = PartialWafersStudies(FLAGS)

    main()
