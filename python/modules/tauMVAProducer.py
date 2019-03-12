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
	#self.bdt_file = "/eos/uscms/store/user/mkilpatr/13TeV/tauMVA/xgboost.xml"
	self.bdt_file = "/uscms_data/d3/mkilpatr/CMSSW_10_2_9/src/TauMVATraining/tauMVA-xgb.model"
	self.bdt_vars = ["Pfcand_pt", "Pfcand_eta", "Pfcand_dz", "Pfcand_chiso0p1", "Pfcand_chiso0p2", "Pfcand_chiso0p3", "Pfcand_chiso0p4", "Pfcand_totiso0p1", "Pfcand_totiso0p2", "Pfcand_totiso0p3", "Pfcand_totiso0p4", "Pfcand_nearestTrkDR", "Pfcand_contjetdr", "Pfcand_contjetcsv"]
	self.xgb = XGBHelper(self.bdt_file, self.bdt_vars)

    def beginJob(self,histFile=None,histDirName=None):
   	pass
    def endJob(self):
	pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("kinBDT", "F", lenVar="nPFcand")
	#self.out.branch("TauMVA_Stop0l", "O", lenVar="nPFcand");

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

	mva = {}
	mva_ = []
	for tau in taus :
		jetmatch = tau.contJetIndex > -1 and jet[tau.contJetIndex].pt > 30.0 and math.fabs(jet[tau.contJetIndex].eta) < 2.4;
		contjetdr = deltaR(tau, jet[tau.contJetIndex]) if jetmatch else -1.0;
		contjetcsv = jet[tau.contJetIndex].btagDeepB if jetmatch else -1.0;
		mva = {"Pfcand_pt": tau.pt, 
		       "Pfcand_eta": tau.eta, 
		       "Pfcand_dz": tau.dz, 
		       "Pfcand_chiso0p1": tau.chiso0p1, 
		       "Pfcand_chiso0p2": tau.chiso0p2, 
		       "Pfcand_chiso0p3": tau.chiso0p3, 
		       "Pfcand_chiso0p4": tau.chiso0p4, 
		       "Pfcand_totiso0p1": tau.totiso0p1, 
		       "Pfcand_totiso0p2": tau.totiso0p2, 
		       "Pfcand_totiso0p3": tau.totiso0p3, 
		       "Pfcand_totiso0p4": tau.totiso0p4, 
		       "Pfcand_nearestTrkDR": tau.nearestTrkDR, 
		       "Pfcand_contjetdr": contjetdr, 
		       "Pfcand_contjetcsv": contjetcsv}
		mva_.append(self.xgb.eval(mva))
	#print "mva output: ", mva_
	self.out.fillBranch("kinBDT", mva_)
	#self.out.fillBranch("TauMVA_Stop0l", True if mva_ > self.tauMVADisc else False)

        return True
