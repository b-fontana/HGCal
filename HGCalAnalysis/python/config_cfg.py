import os, sys, glob
import FWCore.ParameterSet.Config as cms

#arguments parsing
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.register ('pu',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "string")
options.parseArguments()
print("************")
for k,v in options.__dict__["_singletons"].items():
    print("{}: {}".format(k,v))
print("************")

#package loading
process = cms.Process("postRECO")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.Geometry.GeometryExtended2023D28_cff')
process.load('Configuration.Geometry.GeometryExtended2023D28Reco_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
"""
process.MessageLogger = cms.Service("MessageLogger",
                                    destinations = cms.untracked.vstring(
                                        'detailedInfo'),
                                    detailedInfo = cms.untracked.PSet(
                                        threshold = cms.untracked.string('INFO'),
                                        default = cms.untracked.PSet(
                                            limit = cms.untracked.int32(-1))),
                                    debugModules = cms.untracked.vstring('*'))
"""

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
from RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi import *

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

indir = "/eos/cms/store/group/dpg_hgcal/comm_hgcal/sitong/D41/photon_flatE/FlatRandomEGunProducer_sitong_20190516/RECO/" if not options.pu else "/eos/cms/store/relval/CMSSW_10_6_0/RelValPhotonGunPt8To150/GEN-SIM-RECO/PU25ns_106X_upgrade2023_realistic_v2_2023D41PU200-v1/10000/"
fNames = ["file:" + it for it in glob.glob(os.path.join(indir,"*.root"))[:90]]

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(*fNames),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)

process.prod = cms.EDProducer('HGCalRecHitsMask',
                              Mask = cms.uint32(3))

process.an = cms.EDAnalyzer("HGCROIAnalyzer",
                            dEdXWeights = dEdX.weights,
                            thicknessCorrection = cms.vdouble(1.132,1.092,1.084),
                            byClosest = cms.bool(False))

outdir = "/eos/user/b/bfontana/HGCalAnalysis/"
pu_str = "pu" if options.pu else "nopu"
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outdir+pu_str+"mask.root"))
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string(outdir+pu_str+"output.root"))

process.p = cms.Path(process.prod * process.an)
process.outpath = cms.EndPath(process.out)
