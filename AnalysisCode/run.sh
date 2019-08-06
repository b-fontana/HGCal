#!/usr/bin/env bash
ARGS=`getopt -o "a" -l "analysis,mask:,samples:,method:,apply_weights,all,compile" -n "getopts_${0}" -- "$@"`

#Bad arguments
if [ $? -ne 0 ];
then
  exit 1
fi

eval set -- "$ARGS"
fw=0
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

    --apply_weights)
      fw=1;
      echo "Weights being applied."; 
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
if [ "${fcompile}" -eq 0 ]; then
      echo "Warning: The code was not compiled!"; 
fi
echo "##########################"

{ #the braces protect agains run-time file editing

###Energy distribution###
if [ "${METHOD}" = "ed" ]; then
    if [ "${fcompile}" -eq 1 ]; then
	g++ -std=c++17 -o exe.o `root-config --libs --cflags` analysis_ed.cc src/calibration.cc src/utils.cc src/software_correction.cc src/parser.cc
	g++ -std=c++17 -o exe2.o `root-config --libs --cflags` analysis_ed_weights.cc src/calibration.cc src/utils.cc src/software_correction.cc src/parser.cc
    fi
    if [ "${fall}" = 1 ]; then
	if [ "${fw}" = 1 ]; then
	    ./exe2.o --samples "${SAMPLES}" --mask 3
	    ./exe2.o --samples "${SAMPLES}" --mask 4	
	    ./exe2.o --samples "${SAMPLES}" --mask 5
	    ./exe2.o --samples "${SAMPLES}" --mask 6
	else
	    ./exe.o --samples "${SAMPLES}" --mask 3 
	    ./exe2.o --samples "${SAMPLES}" --mask 3
	    echo "Mask 3 finished."
	    ./exe.o --samples "${SAMPLES}" --mask 4 
	    ./exe2.o --samples "${SAMPLES}" --mask 4	
	    echo "Mask 4 finished."
	    ./exe.o --samples "${SAMPLES}" --mask 5 
	    ./exe2.o --samples "${SAMPLES}" --mask 5
	    echo "Mask 5 finished."
	    ./exe.o --samples "${SAMPLES}" --mask 6 
	    ./exe2.o --samples "${SAMPLES}" --mask 6
	    echo "Mask 6 finished."
	fi
    else
	if [ "${fw}" = 1 ]; then
	    ./exe2.o --samples "${SAMPLES}" --mask "${MASK}"
	else
	    ./exe.o --samples "${SAMPLES}" --mask "${MASK}"
	    ./exe2.o --samples "${SAMPLES}" --mask "${MASK}"
	fi
    fi
###Fine eta - calibrate al events###
elif [ "${METHOD}" = "fineeta" ]; then
    echo "Not implemented (yet)."
fi
}
exit $? #protects against run-time appends
