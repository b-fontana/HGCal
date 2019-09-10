#!/bin/bash

varExists() { 
    # Checks whether a certain environment variable already exists
    # Arguments:
    # 1. Variable being checked
    local flag=false;
    if [ -z "${1}" ]; then
	flag=true;
    fi
    echo $flag;
}

if [ $(varExists "${INIT_FOLDER}") = true ] && [ $(varExists "${CMSSW_PATH}") = true ] &&
    [ $(varExists "${HOME_DIR}") = true ] && [ $(varExists "${FULL_PATH}") = true ] && 
    [ $(varExists "${CONFIG_FILE}") = true ]; then
    INIT_FOLDER=$(pwd);
    CMSSW_PATH="CMSSW_10_6_0/src/";
    HOME_DIR="/afs/cern.ch/user/b/bfontana/";
    FULL_PATH="${HOME_DIR}""${CMSSW_PATH}";
    CONFIG_FILE="${FULL_PATH}UserCode/HGCalMaskResolutionAna/test/MaskResolutionAna_cfg.py";
else
    echo "Use different variable names.";
    exit 0;
fi

export XRD_NETWORKSTACK=IPv4
export CMSSWVER="CMSSW_10_6_0"
export SCRAM_ARCH="slc7_amd64_gcc700"

cd "${FULL_PATH}";

source /afs/cern.ch/cms/cmsset_default.sh
eval `scramv1 runtime -sh` #cmsenv substitute

#back to the job folder
cd "${INIT_FOLDER}";
pwd

if [ -e "${CONFIG_FILE}" ]; then
    echo "Before command";
    cmsRun "${CONFIG_FILE}" pu=0 fidx="${1}" mask="${2}" samples="${3}"
    echo "After command";
else
    echo "The configuration file was not found. End of program.";
    exit 0;
fi

outfile="${1}_mask${2}_${3}_nopu"
if [ -r "${outfile}.root" ]; then
    mv "${outfile}.root" /eos/user/b/bfontana/HGCalMaskResolution/Pions/mask"${2}"_"${3}"/;
    :
else
    echo "File ${outfile}.root was not produced by the configuration file.";
    exit 0;
fi
