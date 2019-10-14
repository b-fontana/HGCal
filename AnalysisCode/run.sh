#!/usr/bin/env bash                                                                                                                                                                                                                                                                      
ARGS=`getopt -o "sca" -l "samples:,method:,region:" -n "getopts_${0}" -- "$@"`

#Bad arguments                                                                                                                                                                                                                                                                           
if [ $? -ne 0 ];
then
  exit 1
fi

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

for i in {3..6}; do 
    photons_exe --mask "$i" --samples "$SAMPLES" --method "$METHOD";
    python plotting/uproot_resolution.py "$i" "$SAMPLES" "$METHOD" 1;
done
python plotting/final.py "$METHOD" "$SAMPLES" "$REGION"
