#!/bin/bash
INIT_FOLDER=$(pwd);
HOME_PATH="/afs/cern.ch/user/b/bfontana/"
CMS_PATH="CMSSW_10_6_0/src/"
FULL_PATH="${HOME_PATH}""${CMS_PATH}"

export XRD_NETWORKSTACK=IPv4
export CMSSWVER="CMSSW_10_6_0"
export SCRAM_ARCH="slc6_amd64_gcc700"

cd "${FULL_PATH}";

source /afs/cern.ch/cms/cmsset_default.sh
eval `scramv1 runtime -sh` #cmsenv substitute

cd UserCode/HGCalAnalysis/;
scram b -j 8;

cmsRun test/config_cfg.py pu=True fidx="${1}"
cd "${INIT_FOLDER}"
