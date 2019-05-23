##Installation

```
cmsrel CMSSW_10_6_0
cd CMSSW_10_6_0/src/
cmsenv
git clone git@github.com:CMS-HGCAL/reco-ntuples.git RecoNtuples
git clone git@github.com:b-fontana/HGCal.git UserCode
scram b -j 8
```

## Region Of Interest analyzer

The ROIanalyzer is currently set to analyze H->gg events but can be easily changed for other topologies.
To create a summary ntuple including the ROI with 3 different SR sizes in HGCAL around the generator level particles
and additional 5 control regions for noise/pileup from a rotation in phi at the same eta run:

```
cmsRun test/roiAnalysisConfig.py maxEvents=-1 inputFiles=/store/cmst3/user/psilva/HGCal/H125gg_EE/CMSSW_9_3_2/DIGI_PU0_0p0/RECO/ outputFile=ROISummary_PU0_0p0.root;
```

A second script is used to collect the energy in each ROI and in the associated noise control regions.
It produces a small ntuple with the basic inputs for calibration.

```
python scripts/summarizeROIforCalibration.py ROISummary_PU0_0p0.root ROISummary_PU0_0p0_forcalib.root 
```

The calibration can be run as follows (first argument is the no pileup file, second argument is the pileup file, last argument is a tag for the calibration).
The L0 (relative - eta), L1 (absolute - E), L2 (pileup - average noise) calibration is stored in local pickle file in a dict.

```
python scripts/runROICalibration.py ROISummary_PU0_0p0_forcalib.root ROISummary_PU140_0p0_forcalib.root 140
```