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
F.register('mask',
           -1,
           F.multiplicity.singleton,
           F.varType.int,
           "Mask to be used. Accepted values: 3, 4 or 5.")
F.register('samples',
           '',
           F.multiplicity.singleton,
           F.varType.string,
           'Which samples to use. Inner ("inner"), outer ("outer"), or both ("all").')
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

indir1 = "/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/FlatRandomEGunProducer_bfontana_20190531/RECO/" #Inner radii
indir2 = "/eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/FlatRandomEGunProducer_bfontana_outer_20190605/RECO/" #Outer radii
glob1 = glob.glob(os.path.join(indir1,"*.root"))
glob2 = glob.glob(os.path.join(indir2,"*.root"))
if F.samples == 'all':
    glob_tot = glob1 + glob2
elif F.samples == 'inner':
    glob_tot = glob1
elif F.samples == 'outer':
    glob_tot = glob2
else:
    raise ValueError('Insert a valid "samples" option!')
fNames = ["file:" + it for it in glob_tot][F.fidx]

if isinstance(fNames,list):     
    process.source = cms.Source("PoolSource",
                        fileNames = cms.untracked.vstring(fNames),
                        duplicateCheckMode = cms.untracked.string("noDuplicateCheck"))
else:
    process.source = cms.Source("PoolSource",
                        fileNames = cms.untracked.vstring(fNames),
                        duplicateCheckMode = cms.untracked.string("noDuplicateCheck"))

process.prod = cms.EDProducer('HGCalMaskProd',
                              LayersAnalysed = cms.vuint32(1,2),
                              Mask = cms.uint32(F.mask))

process.an = cms.EDAnalyzer("HGCalMaskResolutionAna",
                            recHitCollection = cms.InputTag("HGCalRecHit","HGCEERecHits"),
                            dEdXWeights = dEdX.weights,
                            thicknessCorrection = cms.vdouble(1.132,1.092,1.084),
                            byClosest = cms.bool(False))
process.an_mask = process.an.clone(
    recHitCollection = cms.InputTag("RecHitsMasked","HGCEERecHits"))

pu_str = "pu" if F.pu else "nopu"
#outsubdir = F.outdir+'mask'+str(F.mask)+'_'+F.samples+'_weight/'
#f not os.path.isdir(outsubdir):
outsubdir = F.outdir
fileName = str(F.fidx)+"_mask"+str(F.mask)+"_"+F.samples+"_"+pu_str
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(fileName+".root"))
fileName = fileName + '_out'
process.out = cms.OutputModule("PoolOutputModule", 
                               fileName = cms.untracked.string(fileName+".root"))

process.p = cms.Path(process.prod * process.an_mask)
process.outpath = cms.EndPath(process.out)
