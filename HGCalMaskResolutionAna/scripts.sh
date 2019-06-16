#!/usr/bin/env bash
samples=inner
for i in {3..6}; do 
    output_file=summary_mask"${i}"_"${samples}".root
    python scripts/summarizeROIforCalibration.py /eos/user/b/bfontana/HGCalMaskResolution/mask"${i}"_"${samples}"/hadd_mask"${i}"_"${samples}"_nopu.root "${output_file}"; 
    python scripts/runROICalibration.py --noPUFile "${output_file}" --outpath pics_"${samples}"/mask"${i}"
done;
