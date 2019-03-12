import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import numpy as np

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest

#2016 MC: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco#Data_MC_Scale_Factors_period_dep
#2017 MC: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation94X

class LLObjectsProducer(Module):
    def __init__(self):
        self.metBranchName = "MET"
	self.p_tauminus = 15
	self.p_Z0       = 23
	self.p_Wplus    = 24

    def beginJob(self):
        pass
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
	self.out.branch("Stop0l_MtEleMET", "F",  lenVar="nElectron")
	self.out.branch("Stop0l_MtMuonMET", "F", lenVar="nMuon")
	self.out.branch("Stop0l_nElectron","I")
	self.out.branch("Stop0l_nMuon",    "I")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def SelEle(self, ele):
	#print "ele pt: %d, ele eta: %d", ele.pt, ele.eta
        if math.fabs(ele.eta) > 2.5 or ele.pt < 5:
            return False
        ## Veto ID electron
        if ele.cutBasedNoIso < 1:
            return False
        ## MiniIso < 0.1
        if ele.miniPFRelIso_all > 0.1:
            return False
        return True

    def SelMuon(self, mu):
        ## NanoAOD store loose ID Muon by default
        if math.fabs(mu.eta) > 2.4 or mu.pt < 5:
            return False
        ## MiniIso < 0.1
        if mu.miniPFRelIso_all > 0.2:
            return False
        return True

    def SelMtlepMET(self, ele, muon, met):
	mtele = []
	mtmuon = []
	for l in ele:
		mtele.append(math.sqrt( 2 * met.pt * l.pt * (1 - math.cos(met.phi-l.phi))))
	for l in muon:
		mtmuon.append(math.sqrt( 2 * met.pt * l.pt * (1 - math.cos(met.phi-l.phi))))
	return mtele, mtmuon

    def isA(self, particleID, p):
	return abs(p) == particleID

    def SelTau(self, genpart, pfc):
	for g in genpart:
		genPartMom = g.genPartIdxMother
		if self.isA(self.p_tauminus, g.pdgId) and deltaR(g.eta, g.phi, pfc.eta, pfc.phi) < 0.2:# and (self.isA(self.p_Z0, genPartMom) or self.isA(self.p_Wplus, genPartMom)):
			return True
	return False

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        ## Getting objects
	jets	  = Collection(event, "Jet")
        electrons = Collection(event, "Electron")
        muons     = Collection(event, "Muon")
        met       = Object(event, self.metBranchName)
	genpart   = Collection(event, "GenPart")

        ## Selecting objects
        self.Electron_Stop0l = map(self.SelEle, electrons)
        self.Muon_Stop0l     = map(self.SelMuon, muons)

        ## Jet variables
	MtEleMET, MtMuonMET = self.SelMtlepMET(electrons, muons, met)
	
        ### Store output
	self.out.fillBranch("Stop0l_MtEleMET",  MtEleMET)
	self.out.fillBranch("Stop0l_MtMuonMET", MtMuonMET)
	self.out.fillBranch("Stop0l_nElectron",sum(self.Electron_Stop0l))
	self.out.fillBranch("Stop0l_nMuon",    sum(self.Muon_Stop0l))
	return True


 # define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
