#!/usr/bin/env bash
#Etacuts: 1.5 outer; 2.95 inner
#Arguments: 1: mask; 2: samples; 3: mode; 4: apply_weights 
output_file=summary_mask"${1}"_"${2}";
#python scripts/summarizeROIforCalibration.py --noPUFile /eos/user/b/bfontana/HGCalMaskResolution/mask"${1}"_"${2}"/hadd_mask"${1}"_"${2}"_nopu --outpath "${output_file}";
if [ "${2}" == "inner" ]; then
    if [ "${4}" == "False" ]; then
    python scripts/runROICalibration.py --noPUFile "${output_file}" --samples "${2}" --mode "${3}" --mingeneta 1.6 --maxgeneta 2.9 --etacuts 2.9 2.93 2.96 3.6 --outpath pics_"${2}"/mask"${1}" --mask "${1}";
    elif [ "${4}" == "True" ]; then
	python scripts/runROICalibration.py --noPUFile "${output_file}" --samples "${2}" --mode "${3}" --mingeneta 1.6 --maxgeneta 2.9 --etacuts 2.9 2.93 2.96 3.6 --outpath pics_"${2}"/mask"${1}" --mask "${1}" --apply_weights;
    fi
elif [ "${2}" == "outer" ]; then
    if [ "${4}" == "False" ]; then
	python scripts/runROICalibration.py --noPUFile "${output_file}" --samples "${2}" --mode "${3}" --mingeneta 1.6 --maxgeneta 1.9 --etacuts 1.35 1.5 1.52 1.6 --outpath pics_"${2}"/mask"${1}" --mask "${1}";
    elif [ "${4}" == "True" ]; then
	python scripts/runROICalibration.py --noPUFile "${output_file}" --samples "${2}" --mode "${3}" --mingeneta 1.6 --maxgeneta 1.9 --etacuts 1.35 1.5 1.52 1.6 --outpath pics_"${2}"/mask"${1}" --mask "${1}" --apply_weights;
    fi
fi
