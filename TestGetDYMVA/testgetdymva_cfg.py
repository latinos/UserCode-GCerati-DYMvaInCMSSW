import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START50_V15::All"
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/smurf/cerati/Sync2012/00D848CB-F76F-E111-A0DE-0024E8769B39.root'
    )
)

process.demo = cms.EDAnalyzer('TestGetDYMVA')

process.p = cms.Path(process.demo)
