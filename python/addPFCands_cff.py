import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from RecoJets.JetProducers.ak8PFJets_cfi  import ak8PFJetsPuppi,ak8PFJetsPuppiSoftDrop
from RecoJets.JetProducers.ak8PFJetsPuppi_groomingValueMaps_cfi import ak8PFJetsPuppiSoftDropMass
from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
from PhysicsTools.PatAlgos.tools.helpers  import getPatAlgosToolsTask, addToProcessAndTask
from PhysicsTools.PatAlgos.tools.jetCollectionTools import RecoJetAdder
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

def ReclusterLSJets(proc, recoJA, runOnMC):
    recoJetInfo = recoJA.addRecoJetCollection(proc, 
                                              "ak8pfpuppi",
                                              inputCollection = "",
                                              genJetsCollection = cms.InputTag("slimmedGenJetsAK8"),
                                              bTagDiscriminators = ['pfDeepCSVJetTags:probb',
                                                                    'pfDeepCSVJetTags:probbb',
                                                                    'pfBoostedDoubleSecondaryVertexAK8BJetTags'],
                                          )

def _addProcessAndTask(proc, label, module):
    from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask, addToProcessAndTask
    task = getPatAlgosToolsTask(proc)
    addToProcessAndTask(label, module, proc, task)

def addPFCands(process, runOnMC=False, allPF = False, onlyAK4=False, onlyAK8=False, addLS=False):
    process.customizedPFCandsTask = cms.Task( )
    process.schedule.associate(process.customizedPFCandsTask)

    process.finalJetsAK8Constituents = cms.EDProducer("PatJetConstituentPtrSelector",
                                            src = cms.InputTag("finalJetsAK8"),
                                            cut = cms.string("")
                                            )
    process.finalJetsAK4Constituents = cms.EDProducer("PatJetConstituentPtrSelector",
                                            src = cms.InputTag("finalJets"),
                                            cut = cms.string("")
                                            )
    if allPF:
        candInput = cms.InputTag("packedPFCandidates")
    elif onlyAK4:
        candList = cms.VInputTag(cms.InputTag("finalJetsAK4Constituents", "constituents"))
        process.customizedPFCandsTask.add(process.finalJetsAK4Constituents)
        process.finalJetsConstituents = cms.EDProducer("PackedCandidatePtrMerger", src = candList, skipNulls = cms.bool(True), warnOnSkip = cms.bool(True))
        candInput = cms.InputTag("finalJetsConstituents")
    elif onlyAK8:
        candList = cms.VInputTag(cms.InputTag("finalJetsAK8Constituents", "constituents"))
        process.customizedPFCandsTask.add(process.finalJetsAK8Constituents)
        process.finalJetsConstituents = cms.EDProducer("PackedCandidatePtrMerger", src = candList, skipNulls = cms.bool(True), warnOnSkip = cms.bool(True))
        candInput = cms.InputTag("finalJetsConstituents")
    else:
        candList = cms.VInputTag(cms.InputTag("finalJetsAK4Constituents", "constituents"), cms.InputTag("finalJetsAK8Constituents", "constituents"))
        process.customizedPFCandsTask.add(process.finalJetsAK4Constituents)
        process.customizedPFCandsTask.add(process.finalJetsAK8Constituents)
        process.finalJetsConstituents = cms.EDProducer("PackedCandidatePtrMerger", src = candList, skipNulls = cms.bool(True), warnOnSkip = cms.bool(True))
        candInput = cms.InputTag("finalJetsConstituents")
    process.customConstituentsExtTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                                                        src = candInput,
                                                        cut = cms.string(""), #we should not filter after pruning
                                                        name = cms.string("PFCands"),
                                                        doc = cms.string("interesting particles from AK4 and AK8 jets"),
                                                        singleton = cms.bool(False), # the number of entries is variable
                                                        extension = cms.bool(False), # this is the extension table for the AK8 constituents
                                                        variables = cms.PSet(CandVars,
                                                            puppiWeight = Var("puppiWeight()", float, doc="Puppi weight",precision=10),
                                                            puppiWeightNoLep = Var("puppiWeightNoLep()", float, doc="Puppi weight removing leptons",precision=10),
                                                            vtxChi2 = Var("?hasTrackDetails()?vertexChi2():-1", float, doc="vertex chi2",precision=10),
                                                            trkChi2 = Var("?hasTrackDetails()?pseudoTrack().normalizedChi2():-1", float, doc="normalized trk chi2", precision=10),
                                                            dz = Var("?hasTrackDetails()?dz():-1", float, doc="pf dz", precision=10),
                                                            dzErr = Var("?hasTrackDetails()?dzError():-1", float, doc="pf dz err", precision=10),
                                                            d0 = Var("?hasTrackDetails()?dxy():-1", float, doc="pf d0", precision=10),
                                                            d0Err = Var("?hasTrackDetails()?dxyError():-1", float, doc="pf d0 err", precision=10),
                                                            pvAssocQuality = Var("pvAssociationQuality()", int, doc="primary vertex association quality"),
                                                            lostInnerHits = Var("lostInnerHits()", int, doc="lost inner hits"),
                                                            trkQuality = Var("?hasTrackDetails()?pseudoTrack().qualityMask():0", int, doc="track quality mask"),
                                                         )
                                    )
    process.customAK8ConstituentsTable = cms.EDProducer("PatJetConstituentTableProducer",
                                                        candidates = candInput,
                                                        jets = cms.InputTag("finalJetsAK8"),
                                                        jet_radius = cms.double(0.8),
                                                        name = cms.string("FatJetPFCands"),
                                                        idx_name = cms.string("pFCandsIdx"),
                                                        nameSV = cms.string("FatJetSVs"),
                                                        idx_nameSV = cms.string("sVIdx"),
                                                        )
    process.customAK4ConstituentsTable = cms.EDProducer("PatJetConstituentTableProducer",
                                                        #candidates = cms.InputTag("packedPFCandidates"),
                                                        candidates = candInput,
                                                        jets = cms.InputTag("finalJets"),
                                                        jet_radius = cms.double(0.4),
                                                        name = cms.string("JetPFCands"),
                                                        idx_name = cms.string("pFCandsIdx"),
                                                        nameSV = cms.string("JetSVs"),
                                                        idx_nameSV = cms.string("sVIdx"),
                                                        )


    jetCollection = 'FatJetLSCollection'
    tagName = "FatJetLS"
    postfix = "Recluster" 
    task = getPatAlgosToolsTask(process)
    addToProcessAndTask(jetCollection, ak8PFJetsPuppi.clone(src = cms.InputTag("customAK8ConstituentsTable:pfCandsNoLep")),
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
                     pfCandidates       = cms.InputTag("customAK8ConstituentsTable:pfCandsNoLep"),
                     genJetCollection   = cms.InputTag("slimmedGenJetsAK8"),
                     genParticles       = cms.InputTag('prunedGenParticles'),
                     getJetMCFlavour    = False,
                     btagDiscriminators = ['pfDeepCSVJetTags:probb','pfDeepCSVJetTags:probbb','pfCombinedInclusiveSecondaryVertexV2BJetTags']
    )

    task = getPatAlgosToolsTask(process)
    addToProcessAndTask(jetCollection+'SoftDrop',
                        ak8PFJetsPuppiSoftDrop.clone( src =cms.InputTag("customAK8ConstituentsTable:pfCandsNoLep"),
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
                     fatJets            = cms.InputTag(jetCollection),             # needed for subjet flavor clustering                                                                               
                     groomedFatJets     = cms.InputTag(jetCollection+'SoftDrop') # needed for subjet flavor clustering                                                                                    
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
        #jetSource          = cms.InputTag(selectedPatJetCollection),
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
                                              singleton = cms.bool(False), # the number of entries is variable                                                                                    
                                              extension = cms.bool(False), # this is the extension table for the AK8 constituents                                                                 
                                              variables = cms.PSet(P4Vars,
                                                                   nMuons = Var("?hasOverlaps('muons')?overlaps('muons').size():0", int, doc="number of muons in the jet"),
                                                                   nElectrons = Var("?hasOverlaps('electrons')?overlaps('electrons').size():0", int, doc="number of electrons in the jet"),
                                                                   msoftdropraw = Var("userFloat('FatJetLSCollectionSoftDropMass')", float, doc="raw soft drop mass",precision=10),
                                                                   #msoftdrop = Var("groomedMass('')",float, doc="Corrected soft drop mass with PUPPI",precision=10),
                                                               )
                                          )

    process.customAK8LSsubjetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                                                    src = cms.InputTag(jetCollection+'SoftDrop',"SubJets"),
                                                    cut = cms.string(""), #probably already applied in miniaod
                                                    name = cms.string("FatJetLSSubJet"),
                                                    doc = cms.string("Lepton subtracted fat jets sub jets"),
                                                    singleton = cms.bool(False), # the number of entries is variable
                                                    extension = cms.bool(False), # this is the main table for the jets
                                                    variables = cms.PSet(P4Vars,
                                                                         #btagDeepB = Var("bDiscriminator('pfDeepCSVJetTags:probb')+bDiscriminator('pfDeepCSVJetTags:probbb')",float,doc="DeepCSV b+bb tag discriminator",precision=10),
                                                                         #btagCSVV2 = Var("bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')",float,doc=" pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2)",precision=10),
     
                                                                     )
                                                )

    if not allPF:
        process.customizedPFCandsTask.add(process.finalJetsConstituents)
    process.customizedPFCandsTask.add(process.customConstituentsExtTable)
    process.customizedPFCandsTask.add(process.customAK8ConstituentsTable)
    process.customizedPFCandsTask.add(process.customAK4ConstituentsTable)
    process.customizedPFCandsTask.add(process.customAK8LSTable)
    process.customizedPFCandsTask.add(process.customAK8LSsubjetTable)

    if runOnMC:

        process.genJetsAK8Constituents = cms.EDProducer("GenJetPackedConstituentPtrSelector",
                                                    src = cms.InputTag("slimmedGenJetsAK8"),
                                                    cut = cms.string("pt > 80")
                                                    )

      
        process.genJetsAK4Constituents = process.genJetsAK8Constituents.clone(
                                                    src = cms.InputTag("slimmedGenJets"),
                                                    cut = cms.string("pt > 20")
                                                    )
        if allPF:
            genCandInput = cms.InputTag("packedGenParticles")
        elif onlyAK4:
            genCandList = cms.VInputTag(cms.InputTag("genJetsAK4Constituents", "constituents"))
            genCandInput =  cms.InputTag("genJetsConstituents")
            process.genJetsConstituents = cms.EDProducer("PackedGenParticlePtrMerger", src = genCandList, skipNulls = cms.bool(True), warnOnSkip = cms.bool(True))
        elif onlyAK8:
            genCandList = cms.VInputTag(cms.InputTag("genJetsAK8Constituents", "constituents"))
            genCandInput =  cms.InputTag("genJetsConstituents")
            process.genJetsConstituents = cms.EDProducer("PackedGenParticlePtrMerger", src = genCandList, skipNulls = cms.bool(True), warnOnSkip = cms.bool(True))
        else:
            genCandList = cms.VInputTag(cms.InputTag("genJetsAK4Constituents", "constituents"), cms.InputTag("genJetsAK8Constituents", "constituents"))
            genCandInput =  cms.InputTag("genJetsConstituents")
            process.genJetsConstituents = cms.EDProducer("PackedGenParticlePtrMerger", src = genCandList, skipNulls = cms.bool(True), warnOnSkip = cms.bool(True))
        process.genJetsParticleTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                                                         src = genCandInput,
                                                         cut = cms.string(""), #we should not filter after pruning
                                                         name= cms.string("GenCands"),
                                                         doc = cms.string("interesting gen particles from AK4 and AK8 jets"),
                                                         singleton = cms.bool(False), # the number of entries is variable
                                                         extension = cms.bool(False), # this is the main table for the AK8 constituents
                                                         variables = cms.PSet(CandVars
                                                                          )
                                                     )
        process.genAK8ConstituentsTable = cms.EDProducer("GenJetConstituentTableProducer",
                                                         candidates = genCandInput,
                                                         jets = cms.InputTag("genJetsAK8Constituents"), # Note: The name has "Constituents" in it, but these are the jets
                                                         name = cms.string("GenFatJetCands"),
                                                         nameSV = cms.string("GenFatJetSVs"),
                                                         idx_name = cms.string("pFCandsIdx"),
                                                         idx_nameSV = cms.string("sVIdx"),
                                                         readBtag = cms.bool(False))
        process.genAK4ConstituentsTable = cms.EDProducer("GenJetConstituentTableProducer",
                                                         candidates = genCandInput,
                                                         jets = cms.InputTag("genJetsAK4Constituents"), # Note: The name has "Constituents" in it, but these are the jets
                                                         name = cms.string("GenJetCands"),
                                                         nameSV = cms.string("GenJetSVs"),
                                                         idx_name = cms.string("pFCandsIdx"),
                                                         idx_nameSV = cms.string("sVIdx"),
                                                         readBtag = cms.bool(False))
        process.customizedPFCandsTask.add(process.genJetsAK4Constituents) #Note: For gen need to add jets to the process to keep pt cuts.
        process.customizedPFCandsTask.add(process.genJetsAK8Constituents)
        if not allPF:
            process.customizedPFCandsTask.add(process.genJetsConstituents)
        process.customizedPFCandsTask.add(process.genJetsParticleTable)
        process.customizedPFCandsTask.add(process.genAK8ConstituentsTable)
        process.customizedPFCandsTask.add(process.genAK4ConstituentsTable)
        
    return process
