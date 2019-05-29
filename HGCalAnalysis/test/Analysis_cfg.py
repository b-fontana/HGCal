import os, sys, glob
import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

#arguments parsing
from FWCore.ParameterSet.VarParsing import VarParsing
F = VarParsing('analysis')
F.register('pu',
           0,
           F.multiplicity.singleton,
           F.varType.bool,
           "Whether to run with pile-up.")
F.register('fidx',
           0,
           F.multiplicity.singleton,
           F.varType.int,
           "Which file index to consider.")
F.register('outdir',
           '/eos/user/b/bfontana/HGCalAnalysis/',
           F.multiplicity.singleton,
           F.varType.string,
           "Output directory.")
F.parseArguments()
print("********************")
print("Input arguments:")
for k,v in F.__dict__["_singletons"].items():
    print("{}: {}".format(k,v))
    print("********************")

#package loading
process = cms.Process("postRECO", eras.Phase2C8)
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.Geometry.GeometryExtended2023D41_cff')
process.load('Configuration.Geometry.GeometryExtended2023D41Reco_cff')
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

indir = "/eos/cms/store/group/dpg_hgcal/comm_hgcal/sitong/D41/photon_flatE/FlatRandomEGunProducer_sitong_20190516/RECO/" if not F.pu else "/eos/cms/store/relval/CMSSW_10_6_0/RelValPhotonGunPt8To150/GEN-SIM-RECO/PU25ns_106X_upgrade2023_realistic_v2_2023D41PU200-v1/10000/"
fNames = ["file:" + it for it in glob.glob(os.path.join(indir,"*.root"))[:90]][F.fidx]

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

pu_str = "pu" if F.pu else "nopu"
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(F.outdir+str(F.fidx)+
                                                         "_"+pu_str+".root"))
process.out = cms.OutputModule("PoolOutputModule", 
                               fileName = cms.untracked.string(F.outdir+str(F.fidx)+
                                                               "_"+pu_str+"out.root"))

process.p = cms.Path(process.prod * process.an)
process.outpath = cms.EndPath(process.out)
