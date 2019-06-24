#!/usr/bin/env bash
#Arguments: 1: mask; 2: samples
output_file=summary_mask"${1}"_"${2}";
python scripts/summarizeROIforCalibration.py --noPUFile /eos/user/b/bfontana/HGCalMaskResolution/mask"${1}"_"${2}"/hadd_mask"${1}"_"${2}"_nopu --outpath "${output_file}";
python scripts/runROICalibration.py --noPUFile "${output_file}" --samples "${2}" --ncomponents 2 --outpath pics_"${2}"/mask"${1}";
