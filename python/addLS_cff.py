import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

leptonLessPFProducer = cms.EDProducer('LeptonLessPFProducer',
                                      vertices     = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                      pfParticles  = cms.InputTag('packedPFCandidates'),
                                      muons        = cms.InputTag('slimmedMuons'),
                                      mu_id        = cms.string("CutBasedIdMedium"),
                                      mu_dz        = cms.double(0.1),
                                      mu_d0        = cms.double(0.05),
                                      mu_sip3D     = cms.double(4),
                                      mu_pt        = cms.double(26),
                                      electrons    = cms.InputTag('slimmedElectrons'),
                                      e_id         = cms.string('mvaEleID-Fall17-noIso-V1-wp90'),
                                      e_isoMask    = cms.int32(-1),
                                      e_dz         = cms.double(0.1),
                                      e_d0         = cms.double(0.05),
                                      e_sip3D      = cms.double(4),
                                      e_pt         = cms.double(30)
)

from RecoJets.JetProducers.ak8PFJets_cfi  import ak8PFJetsPuppi,ak8PFJetsPuppiSoftDrop
from RecoJets.JetProducers.ak8PFJetsPuppi_groomingValueMaps_cfi import ak8PFJetsPuppiSoftDropMass
from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
from PhysicsTools.PatAlgos.tools.helpers  import getPatAlgosToolsTask, addToProcessAndTask
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

def _addProcessAndTask(proc, label, module):
    task = getPatAlgosToolsTask(proc)
    addToProcessAndTask(label, module, proc, task)

def addLSJets(process, pfCandLabel="customAK8ConstituentsTable:pfCandsNoLep", runOnMC=False):
    jetCollection = 'FatJetLSCollection'
    tagName = "FatJetLS"
    postfix = "Recluster"
    task = getPatAlgosToolsTask(process)
    addToProcessAndTask(jetCollection, ak8PFJetsPuppi.clone(src = cms.InputTag(pfCandLabel)),
                        process, task)
    getattr(process, jetCollection).jetAlgorithm = 'AntiKt'
    getattr(process, jetCollection).rParam = 0.8

    addJetCollection(process,
                     labelName          = "FatJetLS",
                     postfix            = postfix,
                     jetSource          = cms.InputTag(jetCollection),
                     pvSource           = cms.InputTag("offlineSlimmedPrimaryVertices"),
                     svSource           = cms.InputTag("slimmedSecondaryVertices"),
                     muSource           = cms.InputTag("slimmedMuons"),
                     elSource           = cms.InputTag("slimmedElectrons"),
                     algo               = "AK",
                     rParam             = 0.8,
                     jetCorrections     = ('AK8PFPuppi', ['L2Relative', 'L3Absolute'], 'None'),
                     pfCandidates       = cms.InputTag(pfCandLabel),
                     genJetCollection   = cms.InputTag("slimmedGenJetsAK8"),
                     genParticles       = cms.InputTag('prunedGenParticles'),
                     getJetMCFlavour    = False,
                     btagDiscriminators = ['pfDeepCSVJetTags:probb','pfDeepCSVJetTags:probbb','pfCombinedInclusiveSecondaryVertexV2BJetTags']
    )

    task = getPatAlgosToolsTask(process)
    addToProcessAndTask(jetCollection+'SoftDrop',
                        ak8PFJetsPuppiSoftDrop.clone( src =cms.InputTag(pfCandLabel),
                                                      rParam=0.8,
                                                      jetAlgorithm='AntiKt',
                                                      useExplicitGhosts=True,
                                                      R0= cms.double(0.8),
                                                      zcut=0.1,
                                                      beta=0,
                                                      doAreaFastjet = cms.bool(True),
                                                      writeCompound = cms.bool(True),
                                                      jetCollInstanceName=cms.string('SubJets') ),
                        process, task)

    task = getPatAlgosToolsTask(process)
    addToProcessAndTask(jetCollection+'SoftDropMass',
                        ak8PFJetsPuppiSoftDropMass.clone( src = cms.InputTag(jetCollection),
                                                          matched = cms.InputTag( jetCollection+'SoftDrop'),
                                                          distMax = cms.double(0.8)),
                        process, task)

    getattr(process,"patJetsFatJetLSRecluster").userData.userFloats.src = [jetCollection+'SoftDropMass']

    addJetCollection(process,
                     labelName          = "FatJetLSSubJets",
                     postfix            = postfix,
                     jetSource          = cms.InputTag(jetCollection+'SoftDrop','SubJets'),
                     pvSource           = cms.InputTag("offlineSlimmedPrimaryVertices"),
                     svSource           = cms.InputTag("slimmedSecondaryVertices"),
                     algo               = 'AK',
                     rParam             = 0.8,
                     btagDiscriminators = ['pfDeepCSVJetTags:probb', 'pfDeepCSVJetTags:probbb','pfCombinedInclusiveSecondaryVertexV2BJetTags'],
                     jetCorrections     = ('AK4PFPuppi', ['L2Relative', 'L3Absolute'], 'None'),
                     explicitJTA        = True,  # needed for subjet b tagging                                                                                                                             
                     svClustering       = True, # needed for subjet b tagging                                                                                                                              
                     genJetCollection   = cms.InputTag('slimmedGenJetsAK8SoftDropSubJets'),
                     genParticles       = cms.InputTag('prunedGenParticles'),
                     getJetMCFlavour    = False,
                     fatJets            = cms.InputTag(jetCollection), 
                     groomedFatJets     = cms.InputTag(jetCollection+'SoftDrop')
    )

    selectedPatJetCollection = "selectedPatJets{}{}".format(tagName,postfix)
    task = getPatAlgosToolsTask(process)
    addToProcessAndTask(selectedPatJetCollection+"Packed",
                        cms.EDProducer("BoostedJetMerger",
                                       jetSrc=cms.InputTag(selectedPatJetCollection),
                                       subjetSrc=cms.InputTag("patJetsFatJetLSSubJetsRecluster")),
                        process, task)

    task = getPatAlgosToolsTask(process)
    addToProcessAndTask('packedPatJets'+'FatJetLS',
                        cms.EDProducer("JetSubstructurePacker",
                                       jetSrc=cms.InputTag(selectedPatJetCollection),
                                       distMax = cms.double(0.8),
                                       fixDaughters = cms.bool(False),
                                       algoTags = cms.VInputTag(
                                           cms.InputTag(selectedPatJetCollection+"Packed")
                                       ),
                                       algoLabels =cms.vstring('SoftDrop')),
                        process, task)
    
    updateJetCollection(
        process,
        labelName          = "FatJetLS",
        postfix            = "Final",
        jetSource          = cms.InputTag('packedPatJets'+'FatJetLS'),
        rParam             = 0.8,
        jetCorrections     = ('AK8PFPuppi', ['L2Relative', 'L3Absolute'], 'None'),
        btagDiscriminators = ['pfDeepCSVJetTags:probb','pfDeepCSVJetTags:probbb','pfBoostedDoubleSecondaryVertexAK8BJetTags',
                              'pfCombinedInclusiveSecondaryVertexV2BJetTags',
                          ]
    )

    patJetFinalCollection="selectedUpdatedPatJets{}{}".format(tagName,"Final")

    process.customAK8LSTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                                              src= cms.InputTag("selectedUpdatedPatJetsFatJetLSFinal"),
                                              cut = cms.string(""),
                                              name = cms.string("FatJetLS"),
                                              doc = cms.string("Lepton subtracted fat jets"),
                                              singleton = cms.bool(False),
                                              extension = cms.bool(False),
                                              variables = cms.PSet(P4Vars,
                                                                   msoftdropraw = Var("userFloat('FatJetLSCollectionSoftDropMass')", float, doc="raw soft drop mass",precision=10),
                                                               )
                                          )

    process.customAK8LSsubjetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                                                    src = cms.InputTag(jetCollection+'SoftDrop',"SubJets"),
                                                    cut = cms.string(""),
                                                    name = cms.string("FatJetLSSubJet"),
                                                    doc = cms.string("Lepton subtracted fat jets sub jets"),
                                                    singleton = cms.bool(False),
                                                    extension = cms.bool(False),
                                                    variables = cms.PSet(P4Vars,
                                                                     )
                                                )

def addNoLep(process, runOnMC=False):
    process.customizedNoLepTask = cms.Task( )
    process.schedule.associate(process.customizedNoLepTask)

    # addLSJets(process, pfCandLabel="customAK8ConstituentsTable:pfCandsNoLep", runOnMC=False)                                                                                                                 # process.customizedNoLepTask.add(process.customAK8LSTable)                                                                                                                                           
    # process.customizedNoLepTask.add(process.customAK8LSsubjetTable)                                                                                                                                     

    def producePF(process) :
        from CommonTools.PileupAlgos.Puppi_cff import puppi
        puppi.useExistingWeights = True
        puppi.candName = "packedPFCandidates"
        puppi.vertexName = "offlineSlimmedPrimaryVertices"
        _addProcessAndTask(process,"leptonLessPFProducer",leptonLessPFProducer.clone())
        _addProcessAndTask(process,"leptonLesspuppi",puppi.clone(candName = cms.InputTag('leptonLessPFProducer')))

    def LSjetsequence(process,runOnMC):
        _btagDiscriminators = [
            'pfJetProbabilityBJetTags',
            'pfDeepCSVJetTags:probb',
            'pfDeepCSVJetTags:probc',
            'pfDeepCSVJetTags:probbb',
            'pfDeepCSVJetTags:probudsg',
            'pfDeepBoostedDiscriminatorsJetTags:WvsQCD',
            'pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:WvsQCD'
        ]
        from RecoBTag.ONNXRuntime.pfParticleNet_cff import _pfParticleNetJetTagsAll as pfParticleNetJetTagsAll
        _btagDiscriminators += pfParticleNetJetTagsAll
        subjetBTagDiscriminators = ['pfDeepCSVJetTags:probb','pfDeepCSVJetTags:probbb']
        JETCorrLevels = ['L1FastJet','L2Relative','L3Absolute']
        if not runOnMC: JETCorrLevels.append('L2L3Residual')
        
        from PhysicsTools.PFNano.jetToolbox_cff import jetToolbox
        jetToolbox(process, 'ak8',  'leptonSubtractedJetSequence','noOutput',PUMethod = 'Puppi', postFix='NoLep',
                   newPFCollection=True, nameNewPFCollection='leptonLesspuppi',
                   JETCorrPayload = 'AK8PFPuppi', JETCorrLevels = JETCorrLevels,
                   Cut='pt > 140.0 && abs(rapidity()) < 2.4', dataTier='miniAOD',
                   runOnMC=runOnMC, addSoftDrop=True, addSoftDropSubjets=True,
                   addNsub=True, maxTau=4,
                   GetSubjetMCFlavour=False,GetJetMCFlavour=False,
                   subJETCorrPayload='AK4PFPuppi',subJETCorrLevels = JETCorrLevels, 
                   bTagDiscriminators=['None'],subjetBTagDiscriminators=['None']
        )
        updateJetCollection(
            process,
            jetSource=cms.InputTag('packedPatJetsAK8PFPuppiNoLepSoftDrop'),
            pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
            svSource = cms.InputTag('slimmedSecondaryVertices'),
            rParam=0.8,
            jetCorrections = ('AK8PFPuppi', cms.vstring(JETCorrLevels), 'None'),
            btagDiscriminators = _btagDiscriminators,
            postfix='AK8NoLepWithPuppiDaughters',   # !!! postfix must contain "WithPuppiDaughter" !!!                                                                                                   
            printWarning = False
        )

    producePF(process)
    LSjetsequence(process,runOnMC)
    process.customAK8LSTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                                              src= cms.InputTag("selectedUpdatedPatJetsAK8NoLepWithPuppiDaughters"),
                                              cut = cms.string(""),
                                              name = cms.string("FatJetLS"),
                                              doc = cms.string("Lepton subtracted fat jets"),
                                              singleton = cms.bool(False),
                                              extension = cms.bool(False),
                                              variables = cms.PSet(P4Vars,
                                                        tau1 = Var("userFloat('NjettinessAK8PuppiNoLep:tau1')", float, doc="Nsubjettiness 1, after subtracting leptons from jet", precision=10),
                                                        tau2 = Var("userFloat('NjettinessAK8PuppiNoLep:tau2')", float, doc="Nsubjettiness 2, after subtracting leptons from jet", precision=10),
                                                        msoftdropraw = Var("userFloat('ak8PFJetsPuppiNoLepSoftDropMass')", float, doc="Raw soft drop", precision=10),
                                                        msoftdrop= Var("groomedMass()", float, doc="Corrected soft drop mass with PUPPI", precision=10),
                                                        deepTag_WvsQCD = Var("bDiscriminator('pfDeepBoostedDiscriminatorsJetTags:WvsQCD')",float,doc="DeepBoostedJet tagger W vs QCD discriminator",precision=10),
                                                        deepTagMD_WvsQCD = Var("bDiscriminator('pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:WvsQCD')",float,doc="Mass-decorrelated DeepBoostedJet tagger W vs QCD discriminator",precision=10),
                                                        particleNet_WvsQCD = Var("bDiscriminator('pfParticleNetDiscriminatorsJetTags:WvsQCD')",float,doc="ParticleNet tagger W vs QCD discriminator",precision=10),
                                                        particleNetMD_Xqq = Var("bDiscriminator('pfMassDecorrelatedParticleNetJetTags:probXqq')",float,doc="Mass-decorrelated ParticleNet tagger raw X->qq (uds) score. For X->qq vs QCD tagging, use Xqq/(Xqq+QCD). For W vs QCD tagging, use (Xcc+Xqq)/(Xcc+Xqq+QCD)",precision=10),
                                                        particleNetMD_Xcc = Var("bDiscriminator('pfMassDecorrelatedParticleNetJetTags:probXqq')",float,doc="Mass-decorrelated ParticleNet tagger raw X->cc score. For X->qq vs QCD tagging, use Xqq/(Xqq+QCD). For W vs QCD tagging, use (Xcc+Xqq)/(Xcc+Xqq+QCD)",precision=10),
                                                        particleNetMD_QCD = Var("bDiscriminator('pfMassDecorrelatedParticleNetJetTags:probQCDbb')+bDiscriminator('pfMassDecorrelatedParticleNetJetTags:probQCDcc')+bDiscriminator('pfMassDecorrelatedParticleNetJetTags:probQCDb')+bDiscriminator('pfMassDecorrelatedParticleNetJetTags:probQCDc')+bDiscriminator('pfMassDecorrelatedParticleNetJetTags:probQCDothers')",float,doc="Mass-decorrelated ParticleNet tagger raw QCD score",precision=10),
                                                               )
                                                            )
    process.customizedNoLepTask.add(process.customAK8LSTable)
    return process
