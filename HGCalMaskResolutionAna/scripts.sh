#!/usr/bin/env bash
#Arguments: 1: mask; 2: samples
output_file=summary_mask"${1}"_"${2}".root
#python scripts/summarizeROIforCalibration.py /eos/user/b/bfontana/HGCalMaskResolution/mask"${1}"_"${2}"/hadd_mask"${1}"_"${2}"_nopu.root "${output_file}"; 
python scripts/runROICalibration.py --noPUFile "${output_file}" --outpath pics_"${2}"/mask"${1}" --samples "${2}"
