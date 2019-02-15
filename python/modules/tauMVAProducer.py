#!/usr/bin/env python
import os, sys
import ROOT
import math
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest
from PhysicsTools.NanoSUSYTools.modules.xgbHelper import XGBHelper

class tauMVAProducer(Module):
    def __init__(self):
	self.writeHistFile=True
	self.tauMVADisc = 0.56
	self.bdt_file = "/eos/uscms/store/user/mkilpatr/13TeV/tauMVA/tauDisc.root"
	self.bdt_vars = ["pt", "abseta", "chiso0p1", "chiso0p2", "chiso0p3", "chiso0p4", "totiso0p1", "totiso0p2", "totiso0p3", "totiso0p4", "neartrkdr", "contjetdr", "contjetcsv"]
	self.xgb = XGBHelper(self.bdt_file, self.bdt_vars)

    def beginJob(self,histFile=None,histDirName=None):
   	pass
    def endJob(self):
	pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.fillBranch("kinBDT", self.xgb.eval(inputs))
	self.out.branch("TauMVA_Stop0l", "I");

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def SelTauMVA(self, mva):
	if mva > self.tauMVADisc:
		return True
	else:
		return False

    def analyze(self, event):
	taus      = Collection(event, "PFcand")
	jet	  = Collection(event, "Jet")

	for tau in taus :
		jetmatch = tau.contJetIndex > -1 and jet[tau.contJetIndex].pt > 30.0 and math.fabs(jet[tau.contJetIndex].eta) < 2.4;
		contjetdr = deltaR(tau, jet[tau.contJetIndex]) if jetmatch else -1.0;
		contjetcsv = jet[tau.contJetIndex].btagDeepB if jetmatch else -1.0;
		mva = [tau.pt, tau.eta, tau.dz, tau.chiso0p1, tau.chiso0p2, tau.chiso0p3, tau.chiso0p4, tau.totiso0p1, tau.totiso0p2, tau.totiso0p3, tau.totiso0p4, tau.nearestTrkDR, contjetdr, contjetcsv]
		print "mva output: ", self.xgb.eval(mva)
		self.out.fillBranch("kinBDT", self.xgb.eval(mva))
		self.out.fillBranch("TauMVA_Stop0l", True if self.xgb.eval(mva) > self.tauMVADisc else False)

        return True
