import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

config_name = 'RecHitsMask'
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('h'+config_name+'.root'))

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('/store/cmst3/group/hgcal/CMG_studies/Production/FlatRandomPtGunProducer_PDGid22_nPart1_Pt100_Eta1p6To2p8_noPU_clange_cmssw1040pre1_20181128/RECO/partGun_PDGid22_x100_Pt100.0To100.0_RECO_12.root'),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"))

process.prod = cms.EDProducer('HGCalMaskProd',
                              LayersAnalysed = cms.vuint32(1,2),
                              Mask = cms.uint32(3))

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('out'+config_name+'.root'))

process.p = cms.Path(process.prod)
process.outpath = cms.EndPath(process.out)
