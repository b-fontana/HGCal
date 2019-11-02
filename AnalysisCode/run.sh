#!/usr/bin/env bash                                                                                                                                                                                                                                                                      
ARGS=`getopt -o "sca" -l "samples:,method:,region:" -n "getopts_${0}" -- "$@"`

#Bad arguments                                                                                                                                                                                                                                                                           
if [ $? -ne 0 ];
then
  exit 1
fi

eval set -- "$ARGS"
echo "##### Input options: #####"
while true;
do
  case "$1" in
    --method)
      if [ -n "$2" ];
      then
        METHOD="${2}";
        echo "Method: ${2}";
      fi
      shift 2;;

    --samples)
      if [ -n "$2" ];
      then
        SAMPLES="${2}";
        echo "Samples: ${2}";
      fi
      shift 2;;

    --region)
      if [ -n "$2" ];
      then
        REGION="${2}";
        echo "Region: ${2}";
      fi
      shift 2;;

    --)
      shift
      break;;
  esac
done

{
#Run the macros related to each mask in parallel
for i in {3..6}; do 
    echo "Running mask " $i;
    photons_exe --mask "$i" --samples "$SAMPLES" --method "$METHOD" && python plotting/uproot_resolution.py "$i" "$SAMPLES" "$METHOD" 0 && python plotting/uproot_resolution.py "$i" "$SAMPLES" "$METHOD" 1 &
done

#Wait for all jobs to complete
pids=()
for i in `jobs -p`; do                                                                                
    pids+=("${i}")                                                          
done
for pid in ${pids[@]}; do                                                                            
    wait ${pid}
done
#Run final plotting script
echo "Producing final plot...";
python plotting/final.py "$METHOD" "$SAMPLES" "$REGION";
echo "Done."
}
exit $? #protects against run-time appends 
