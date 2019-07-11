#!/usr/bin/env bash
ARGS=`getopt -o "sca" -l "summary,calibrate,analysis,mask:,samples:,mode:,apply_weights" -n "getopts_${0}" -- "$@"`

#Bad arguments
if [ $? -ne 0 ];
then
  exit 1
fi

eval set -- "$ARGS"
fs=0
fc=0
fa=0
fw=0
echo "##### Input options: #####"
while true;
do
  case "$1" in
    -s|--summary) 
      fs=1;
      shift;;

    -c|--calibrate)
      fc=1;
      shift;;

    -a|--analysis)
      fa=1;
      shift;;

    --mask)    
      if [ -n "$2" ];
      then
        MASK="${2}";
	echo "Mask: ${2}";
      fi
      shift 2;;

    --samples)    
      if [ -n "$2" ];
      then
        SAMPLES="${2}";
	echo "Samples: ${2}";
      fi
      shift 2;;

    --mode)    
      if [ -n "$2" ];
      then
        MODE="${2}";
	echo "Mode: ${2}";
      fi
      shift 2;;

    --apply_weights)
      fw=1;
      echo "Weights being applied."; 
      shift;;

    --)
      shift
      break;;
  esac
done
echo "##########################"

#Declare variables
OUTFILE=summary_mask"${MASK}"_"${SAMPLES}".root;
if [ "${SAMPLES}" == "inner" ]; then
    MINGENETA=1.6
    MAXGENETA=2.93
    ETACUTS=("${MAXGENETA}" 2.95 2.96 3.6)
elif [ "${SAMPLES}" == "outer" ]; then
    MINGENETA=1.55
    MAXGENETA=1.9
    ETACUTS=(1.35 1.515 1.525 "${MINGENETA}")
fi

#Summarize RECO NTuples information
if [ "${fs}" -eq 1 ]; then
    echo "Summarize!";
    python scripts/summarizeROIforCalibration.py --noPUFile /eos/user/b/bfontana/HGCalMaskResolution/mask"${MASK}"_"${SAMPLES}"/hadd_mask"${MASK}"_"${SAMPLES}"_nopu.root --mingeneta "${MINGENETA}" --maxgeneta "${MAXGENETA}" --samples "${SAMPLES}" --outpath "${OUTFILE}";
fi

#Create files that store the calibration corretion factors
if [ "${fc}" -eq 1 ]; then
    echo "Calibrate!";
    python scripts/runROICalibration.py --noPUFile "${OUTFILE}" --mingeneta "${MINGENETA}" --maxgeneta "${MAXGENETA}" --outpath pics_"${SAMPLES}"/mask"${MASK}" --samples "${SAMPLES}" --mask "${MASK}";
fi

#Run the analysis
if [ "${fa}" -eq 1 ]; then
    echo "Analyze!";
    python scripts/analysis.py --noPUFile summary_mask"${MASK}"_"${SAMPLES}".root --samples "${SAMPLES}" --mode "${MODE}" --mingeneta "${MINGENETA}" --maxgeneta "${MAXGENETA}" --etacuts "${ETACUTS[@]}" --outpath pics_"${SAMPLES}"/mask"${MASK}" --mask "${MASK}" --apply_weights "${fw}"
fi
