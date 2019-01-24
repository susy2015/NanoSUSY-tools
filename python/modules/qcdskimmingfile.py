#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest
#import "$CMSSW_BASE/src/AnalysisTools/QuickRefold/interface/TObjectContainer.h"

class qcdskimmingfile(Module): 
    def __init__(self):
	self.writeHistFile=True
	self.metBranchName="MET"
 
    def beginJob(self,histFile=None,histDirName=None):
   	pass
    def endJob(self):
	pass 

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("met",  "F");
        self.out.branch("met_phi", "F");
	self.out.branch("jet_pt", "F");
	self.out.branch("genjet_pt", "F");
	# self.out.branch("dphij1met", "F");
	# self.out.branch("dphij2met", "F");
	# self.out.branch("dphij3met", "F");
	# self.out.branch("dphij4met", "F");


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def jetResFunction(self, jets, genjets):
	res = jets.pt/genjets.pt
	return res

    def analyze(self, event):
	jets      = Collection(event, "Jet")
	genjets   = Collection(event, "GenJet")
	met       = Object(event,     self.metBranchName)


	# njets = 0	
	# for j in jets :
	# 	njets += 1
	# 	#dphij1met_ = deltaPhi(j, met)
	# 	for gj in genjets :
	# 		origRes_ = self.jetResFunction(j, gj)

	# if njets > 0 : self.out.fillBranch("dphij1met", deltaPhi(jets[0], met))
	# if njets > 1 : self.out.fillBranch("dphij2met", deltaPhi(jets[1], met))
	# if njets > 2 : self.out.fillBranch("dphij3met", deltaPhi(jets[2], met))
	# if njets > 3 : self.out.fillBranch("dphij4met", deltaPhi(jets[3], met))
	#self.out.fillBranch("origRes", origRes_)
        self.out.fillBranch("met", met.pt)
        self.out.fillBranch("met_phi", met.phi)
        return True
