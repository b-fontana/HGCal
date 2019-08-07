#!/usr/bin/env bash
ARGS=`getopt -o "a" -l "analysis,mask:,samples:,method:,first_part,second_part,all,compile" -n "getopts_${0}" -- "$@"`

#Bad arguments
if [ $? -ne 0 ];
then
  exit 1
fi

eval set -- "$ARGS"
fp1=0
fp2=0
fall=0
fcompile=0

echo "##### Input options: #####"
while true;
do
  case "$1" in
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

    --method)    
      if [ -n "$2" ];
      then
        METHOD="${2}";
	echo "Method: ${2}";
      fi
      shift 2;;

    --first_part)
      fp1=1;
      echo "Run only the first part of the 'ed' method."; 
      shift;;

    --second_part)
      fp2=1;
      echo "Run only the second part of the 'ed' method."; 
      shift;;

    --all)
      fall=1;
      echo "Weights being applied."; 
      shift;;

    --compile)
      fcompile=1;
      echo "The code is going to be compiled."; 
      shift;;

    --)
      shift
      break;;
  esac
done
if [ "${fp1}" -eq 0 ] && [ "${fp2}" -eq 0 ] && [ "${METHOD}"=="ed" ]; then
      echo "Running the full ed method."
fi
if [ "${fcompile}" -eq 0 ]; then
      echo "Warning: The code was not compiled!"; 
fi
echo "##########################"

{ #the braces protect agains run-time file editing

###Energy distribution###
EXEC1_NAME="exe1.o"
EXEC2_NAME="exe2.o"
if [ "${METHOD}" = "ed" ]; then
    if [ "${fcompile}" -eq 1 ]; then
	g++ -std=c++17 -o "${EXEC1_NAME}" `root-config --libs --cflags` bin/analysis_ed.cc src/calibration.cc src/utils.cc src/software_correction.cc src/parser.cc
	g++ -std=c++17 -o "${EXEC2_NAME}" `root-config --libs --cflags` bin/analysis_ed_weights.cc src/calibration.cc src/utils.cc src/software_correction.cc src/parser.cc
    fi
    if [ "${fall}" = 1 ]; then
	if [ "${fp1}" = 1 ]; then
	    ./"${EXEC1_NAME}" --samples "${SAMPLES}" --mask 3
	    ./"${EXEC1_NAME}" --samples "${SAMPLES}" --mask 4	
	    ./"${EXEC1_NAME}" --samples "${SAMPLES}" --mask 5
	    ./"${EXEC1_NAME}" --samples "${SAMPLES}" --mask 6
	elif [ "${fp2}" = 1 ]; then
	    ./"${EXEC2_NAME}" --samples "${SAMPLES}" --mask 3
	    ./"${EXEC2_NAME}" --samples "${SAMPLES}" --mask 4	
	    ./"${EXEC2_NAME}" --samples "${SAMPLES}" --mask 5
	    ./"${EXEC2_NAME}" --samples "${SAMPLES}" --mask 6
	else
	    ./"${EXEC1_NAME}" --samples "${SAMPLES}" --mask 3 
	    ./"${EXEC2_NAME}" --samples "${SAMPLES}" --mask 3
	    echo "Mask 3 finished."
	    ./"${EXEC1_NAME}" --samples "${SAMPLES}" --mask 4 
	    ./"${EXEC2_NAME}" --samples "${SAMPLES}" --mask 4	
	    echo "Mask 4 finished."
	    ./"${EXEC1_NAME}" --samples "${SAMPLES}" --mask 5 
	    ./"${EXEC2_NAME}" --samples "${SAMPLES}" --mask 5
	    echo "Mask 5 finished."
	    ./"${EXEC1_NAME}" --samples "${SAMPLES}" --mask 6 
	    ./"${EXEC2_NAME}" --samples "${SAMPLES}" --mask 6
	    echo "Mask 6 finished."
	fi
    else
	if [ "${fp1}" = 1 ]; then
	    ./"${EXEC1_NAME}" --samples "${SAMPLES}" --mask "${MASK}"
	elif [ "${fp2}" = 1 ]; then
	    ./"${EXEC2_NAME}" --samples "${SAMPLES}" --mask "${MASK}"
	else
	    ./"${EXEC1_NAME}" --samples "${SAMPLES}" --mask "${MASK}"
	    ./"${EXEC2_NAME}" --samples "${SAMPLES}" --mask "${MASK}"
	fi
    fi
###Fine eta - calibrate al events###
EXEC3_NAME="exe3.o"
elif [ "${METHOD}" = "fineeta" ]; then
    echo "Not implemented (yet)."
fi
}
exit $? #protects against run-time appends
