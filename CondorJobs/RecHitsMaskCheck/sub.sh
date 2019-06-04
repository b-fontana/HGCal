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
    CMSSW_PATH="/CMSSW_10_6_0/src/";
    HOME_DIR="/afs/cern.ch/user/b/bfontana";
    FULL_PATH="${HOME_DIR}""${CMSSW_PATH}";
    CONFIG_FILE="test/RecHitsMaskStudies_cfg.py";
else
    echo "Use different variable names.";
    exit 0;
fi

export XRD_NETWORKSTACK=IPv4
export CMSSWVER="CMSSW_10_6_0"
export SCRAM_ARCH="slc6_amd64_gcc700"

echo "${FULL_PATH}";
cd "${FULL_PATH}";

source /afs/cern.ch/cms/cmsset_default.sh
eval `scramv1 runtime -sh` #cmsenv substitute

cd UserCode/HGCalRecHitsMaskStudies/;
scram b -j 8;

if [ -e "${CONFIG_FILE}" ]; then
    echo "Before command";
    cmsRun "${CONFIG_FILE}" pu=0 fidx="${1}" mask="${2}";
    echo "After command";
else
    echo "The configuration file was not found. End of program.";
    exit 0;
fi

cd "${INIT_FOLDER}";
