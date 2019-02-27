#!/usr/bin/env python                                                                                                                                   
import os, sys
import ROOT
import math
import numpy as np
from array import array
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest
from PhysicsTools.NanoAODTools.postprocessing.framework.treeReaderArrayTools import *
from rootpy.tree import Tree, TreeModel, IntCol, FloatArrayCol


class JetResSkim(Module):
    def __init__(self):
	pass
    def beginJob(self):
        pass
    def endjob(self):
        pass
    def beginFile(self,inputFile,outputFile,inputTree,wrappedOutputTree, outputFileSkim, outputTreeSkim):
	self.out = wrappedOutputTree
	self.outSkim = outputTreeSkim
	self.outSkim.branch("weight"	,"F")
	self.outSkim.branch("genjetpt"	,"F")
	self.outSkim.branch("genjeteta"	,"F")
	self.outSkim.branch("recojetpt"	,"F")
	self.outSkim.branch("genjetrank"	,"I")
	self.outSkim.branch("flavour"	,"I")
	self.outSkim.branch("rechempt"	,"F")
	self.outSkim.branch("genhempt"	,"F")

    def analyze(self, event):
	jets      = Collection(event, "Jet")
	genjets   = Collection(event, "GenJet")
	weight    = event.genWeight
	eventNum  = event.event
	PassFilter  = event.Pass_EventFilter
	PassJetID   = event.Pass_JetID

	if PassFilter and PassJetID:
		for gJ in xrange(len(genjets)):
			gJet = genjets[gJ]
			if gJet.pt < 20: continue  
			
			rJet = 0
			for iR in xrange(len(jets)) :
				if jets[iR].genJetIdx != gJ: continue
				rJet = jets[iR]
				break
			
			self.outSkim.fillBranch("weight",	weight)
			self.outSkim.fillBranch("genjetpt", 	gJet.pt)
			self.outSkim.fillBranch("genjeteta", 	gJet.eta)
			self.outSkim.fillBranch("recojetpt", 	rJet.pt if rJet != 0 else 9.5)
			self.outSkim.fillBranch("genjetrank", 	min(gJ, 250))
			self.outSkim.fillBranch("flavour",	gJet.partonFlavour)
		
			if(gJet.eta > -2.8 and gJet.eta < -1.6 and gJet.phi >-1.37 and gJet.phi < -1.07):
				self.outSkim.fillBranch("genhempt", gJet.pt)
			else:
				self.outSkim.fillBranch("genhempt", 0)
			if rJet != 0:
				if(rJet.eta > -2.8 and rJet.eta < -1.6 and rJet.phi >-1.37 and rJet.phi < -1.07): 
					self.outSkim.fillBranch("rechempt",rJet.pt)
				else:
					self.outSkim.fillBranch("rechempt", 0)
			else:
				self.outSkim.fillBranch("rechempt", 0)

			self.outSkim.fill()	
	return True  
