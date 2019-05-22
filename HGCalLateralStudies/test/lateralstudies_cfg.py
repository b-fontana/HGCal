import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
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
process.MessageLogger = cms.Service("MessageLogger",
                                    destinations = cms.untracked.vstring(
                                        'detailedInfo'),
                                    detailedInfo = cms.untracked.PSet(
                                        threshold = cms.untracked.string('INFO'),
                                        default = cms.untracked.PSet(
                                            limit = cms.untracked.int32(-1))),
                                    debugModules = cms.untracked.vstring('*'))

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/cmst3/group/hgcal/CMG_studies/Production/FlatRandomPtGunProducer_PDGid22_nPart1_Pt100_Eta1p6To2p8_noPU_clange_cmssw1040pre1_20181128/RECO/partGun_PDGid22_x100_Pt100.0To100.0_RECO_12.root'
    ),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)

process.prod = cms.EDProducer('HGCalLateralStudies',
                              CellUVCoordinates = cms.string("CellUVCoordinates"),
                              WaferUVCoordinates = cms.string("WaferUVCoordinates"),
                              LayersAnalysed = cms.vuint32(1,2),
                              lCellFilterCut = cms.double(17.5*1.5),
                              hCellFilterCut = cms.double(0),
                              Mask = cms.uint32(3))
                              
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("histo_mask.root"))

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string("output.root"))

process.p = cms.Path(process.prod)
process.outpath = cms.EndPath(process.out)
