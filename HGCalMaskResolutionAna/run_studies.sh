#!/usr/bin/env bash
ARGS=`getopt -o "sca" -l "summary,calibrate,analysis,mask:,samples:,mode:,method:,apply_weights,all" -n "getopts_${0}" -- "$@"`

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
fall=0
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

    --method)    
      if [ -n "$2" ];
      then
        METHOD="${2}";
	echo "Method: ${2}";
      fi
      shift 2;;

    --apply_weights)
      fw=1;
      echo "Weights being applied."; 
      shift;;

    --all)
      fall=1;
      echo "Weights being applied."; 
      shift;;

    --)
      shift
      break;;
  esac
done
echo "##########################"

{ #the braces protect agains run-time file editing

#Declare variables
BCKGCUTS=(0.5 0.6 0.7)

#Summarize RECO NTuples information
if [ "${fs}" -eq 1 ]; then
    if [ "${fall}" -eq 1 ]; then
	echo "Summarize all!";
	allsamples=("inner" "outer")
	for imask in {3..6}; do
	    for isample in "${allsamples[@]}"; do
		OUTFILE=summary_mask"${imask}"_"${isample}".root
		python scripts/summarizeROIforCalibration.py --noPUFile /eos/user/b/bfontana/HGCalMaskResolution/mask"${imask}"_"${isample}"/hadd_mask"${imask}"_"${isample}"_nopu.root --samples "${isample}" --outpath "${OUTFILE}"  1>/dev/null &
	    done
	done
    else
	OUTFILE=summary_mask"${MASK}"_"${SAMPLES}".root;
	echo "Summarize!";
	python scripts/summarizeROIforCalibration.py --noPUFile /eos/user/b/bfontana/HGCalMaskResolution/mask"${MASK}"_"${SAMPLES}"/hadd_mask"${MASK}"_"${SAMPLES}"_nopu.root --samples "${SAMPLES}" --outpath "${OUTFILE}";
    fi
fi

#Create files that store the calibration corretion factors
if [ "${fc}" -eq 1 ]; then
    if [ "${fall}" -eq 1 ]; then
	echo "Calibrate all!";
	allsamples=("inner" "outer")
	for imask in {3..6}; do
	    for isample in "${allsamples[@]}"; do
		OUTFILE=summary_mask"${imask}"_"${isample}".root
		python scripts/runROICalibration.py --noPUFile "${OUTFILE}"  --outpath pics_"${isample}"/mask"${imask}" --samples "${isample}" --mask "${imask}" --method "${METHOD}"  1>/dev/null &
	    done
	done
    else
	echo "Calibrate!";
	OUTFILE=summary_mask"${MASK}"_"${SAMPLES}".root;
	python scripts/runROICalibration.py --noPUFile "${OUTFILE}"  --outpath pics_"${SAMPLES}"/mask"${MASK}" --samples "${SAMPLES}" --mask "${MASK}" --method "${METHOD}";
    fi
fi

#Run the analysis
if [ "${fa}" -eq 1 ]; then
    if [ "${fall}" -eq 1 ]; then
	echo "Analyze all!";
	allsamples=("inner" "outer")
	for imask in {3..6}; do
	    for isample in "${allsamples[@]}"; do
		OUTFILE=summary_mask"${imask}"_"${isample}".root
		python scripts/analysis.py --noPUFile summary_mask"${imask}"_"${isample}".root --samples "${isample}" --mode "${MODE}" --bckgcuts "${BCKGCUTS[@]}" --outpath pics_"${isample}"/mask"${imask}" --mask "${imask}" --apply_weights "${fw}" --method "${METHOD}"  1>/dev/null &
	    done
	done
    else
	echo "Analyze!";
	OUTFILE=summary_mask"${MASK}"_"${SAMPLES}".root;
	python scripts/analysis.py --noPUFile summary_mask"${MASK}"_"${SAMPLES}".root --samples "${SAMPLES}" --mode "${MODE}" --bckgcuts "${BCKGCUTS[@]}" --outpath pics_"${SAMPLES}"/mask"${MASK}" --mask "${MASK}" --apply_weights "${fw}" --method "${METHOD}"
    fi
fi

#Warn the user after all the jobs have finished
if [ "${fall}" -eq 1 ]; then
    pids=()
    for i in `jobs -p`; do
	pids+=("${i}")
    done
    for pid in ${pids[@]}; do
	wait ${pid}
    done
    echo "All jobs have finished."
fi
}
exit $? #protects against run-time appends
