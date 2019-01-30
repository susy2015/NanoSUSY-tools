#!/usr/bin/env python
import os, sys
import ROOT
import time
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest
from PhysicsTools.NanoAODTools.postprocessing.framework.treeReaderArrayTools import *
from rootpy.tree import Tree, TreeModel, IntCol, FloatArrayCol


class qcdSmearProducer(Module): 
    def __init__(self):
        self.writeHistFile=True
        self.metBranchName="MET"
        self.xBinWidth = 0.01
        self.minWindow = 0.01
        self.maxWindow = 0.5
        self.nSmears = 100
        self.nSmearJets = 2
        self.nBootstraps = 50
        self.LINEAR_GRANULATED=True
        self.winType = self.LINEAR_GRANULATED
        self.doFlatSampling = True
        self.respInputName = "JetResByFlav"
        self.respFileName = "file:/eos/uscms/store/user/ddash/qcd_smeared/resTailOut_combined_filtered_CHEF_puWeight_weight_WoH_NORMALIZED.root"
	self.outFileName = "Smear_tree_"

    def loadHisto(self,filename,hname):
	tf = ROOT.TFile.Open(filename)
	hist = tf.Get(hname)
	hist.SetDirectory(None)
	tf.Close()
	return hist

    def storeSmear(self, fname):
	inFile = ROOT.TFile.Open(fname)
	inTree = inFile.Get("Events")
	return inFile, inTree

    def storeSmearFile(self, outFileName, inTree, eventNum):
	outFile = ROOT.TFile.Open(outFileName, "recreate")
	outTree = inTree.CloneTree(inTree.GetReadEntry())
	inTree.GetEntry(eventNum)
	outTree.Fill()
	outTree.Write()
	
    def ptmapping(self,jets):
	ptrange = [0, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 700, 1000, 1500]	
        pt_index = -1
	for i in xrange(len(ptrange)):
		if jets.pt < float(ptrange[i]): 
			pt_index = i - 1
			break
	
	bname=["res_b_comp_14","res_b_comp_15","res_b_comp_16","res_b_comp_17","res_b_comp_18","res_b_comp_19","res_b_comp_20","res_b_comp_21","res_b_comp_22","res_b_comp_23","res_b_comp_24","res_b_comp_25","res_b_comp_26"]
        lgtname=["res_light_comp_1","res_light_comp_2","res_light_comp_3","res_light_comp_4","res_light_comp_5","res_light_comp_6","res_light_comp_7","res_light_comp_8","res_light_comp_9","res_light_comp_10","res_light_comp_11","res_light_comp_12","res_light_comp_13"]
        if jets.partonFlavour == 4 :
		return bname[pt_index]
        else :
		return lgtname[pt_index]

    def beginJob(self,histFile=None,histDirName=None):
        pass

    def endJob(self):
        pass 

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
	self.out.branch("Jet_pt", "F", lenVar="nJet")
	self.out.branch("MET_pt", "F")


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def jetResFunction(self, jets, genjets):
        res = jets.pt/genjets.pt
        return res

    def interpolateResToProb(self,cdf,resp):
        xbin = cdf.FindFixBin(resp)
	if xbin <= 0: return 0
	binwidth = cdf.GetBinWidth(xbin)
	b = ((cdf.GetBinWidth(xbin) + cdf.GetBinLowEdge(xbin)) * cdf.GetBinContent(xbin-1) - (cdf.GetBinLowEdge(xbin)) *cdf.GetBinContent(xbin)) / binwidth
	m = (cdf.GetBinContent(xbin) - cdf.GetBinContent(xbin-1))/binwidth
	return resp * m + b

    def interpolteProbToRes(self,cdf,probe):
        binAbove = cdf.FindFirstBinAbove(probe) 
        if binAbove <= 1 :return 0
        deltaProb = cdf.GetBinContent(binAbove) - cdf.GetBinContent(binAbove -1)
        newResValue = cdf.GEtBinCenter(binAbove)
        if deltaProb > 0 :
		b = (cdf.GetBinContent(binAbove) * cdf.GetBinLowEdge(binAbove)- cdf.GetBinContent(binAbove -1 ) * (cdf.GetBinWidth(binAbove) + cdf.GetBinLowEdge(binAbove)) )/deltaprob
		m = cdf.GetBinWidth(binAbove) / deltaProb
		newResValue = m* probe + b
        return newResValue
     
    def getScaledWindow(self,resp,minW,maxW):
        if resp < 1 :
		return (minW - maxW)*resp + maxW
        else :
		return -1 * (minW - maxW) * resp + 2 * minW - maxW 
        
    def getUpIntegratedScaledWindow(self,resp,minW,maxW):
        if resp < 1 - self.getScaledWindow(1,minW,maxW):
		return (resp + maxW)/(1-(minW - maxW))
        else :
		return (resp + 2*minW - maxW)/(1 +(minW - maxW))

    def getLowIntegratedScaledWindow(self,resp,minW,maxW):
        if resp > 1 + self.getScaledWindow(1,minW,maxW):
		return (resp - ( 2* minW - maxW))/(1- (minW - maxW))
        else :
		return (resp - maxW)/(1 + (minW- maxW))

    def getWindowProb(self,cdf,minRes,maxRes):
        minRes = max(0.0001,minRes)
        maxRes = min(1.9999,maxRes)
        if minRes >= maxRes : 
		minProb=0
		maxProb=0
		return
        else :
		minProb = self.interpolateResToProb(cdf,minRes)
		maxProb = self.interpolateResToProb(cdf,maxRes)
	return minProb, maxProb

    def getScaledWindowAndProb(self,cdf,resp,minWindow,maxWindow):
        window = self.getScaledWindow(resp,minWindow,maxWindow)
        minRes = resp - window
        maxRes = resp + window
	minProb, maxProb = self.getWindowProb(cdf,minRes,maxRes)
        return minProb, maxProb, minRes, maxRes

    def getContributionScaledWindowAndProb(self,cdf,resp,minWindow,maxWindow):
        minRes = self.getLowIntegratedScaledWindow(resp,minWindow,maxWindow)
        maxRes = self.getUpIntegratedScaledWindow(resp,minWindow,maxWindow)
        minProb, maxProb = self.getWindowProb(cdf,minRes,maxRes)
	return minProb, maxProb, minRes, maxRes

    def testMetCalc(self, obj1, obj2, obj3):
        tot = ROOT.TLorentzVector()
        v1 = ROOT.TLorentzVector()
        v2 = ROOT.TLorentzVector()
        v3 = ROOT.TLorentzVector()
        v1.SetPtEtaPhiM(obj1.pt, 0, obj1.phi, 0)
        v2.SetPtEtaPhiM(obj2.pt, 0, obj2.phi, 0)
        v3.SetPtEtaPhiM(obj3.pt, 0, obj3.phi, 0)
	tot = (v1 + (v2 - v3))
	return tot
    
    def addFourVector(self,obj1,obj2):
        tot = ROOT.TLorentzVector()
        v1 = ROOT.TLorentzVector()
        v2 = ROOT.TLorentzVector()
        v1.SetPtEtaPhiM(obj1.pt, 0, obj1.phi, 0)
        v2.SetPtEtaPhiM(obj2.pt, 0, obj2.phi, 0)
        tot = v1+v2
        return tot
    
    def addTLorentzVector(self,obj1,obj2):
        tot = ROOT.TLorentzVector()
        v1 = ROOT.TLorentzVector()
        v2 = ROOT.TLorentzVector()
        v1.SetPtEtaPhiM(obj1.Pt(), 0, obj1.Phi(), 0)
        v2.SetPtEtaPhiM(obj2.pt, 0, obj2.phi, 0)
        tot = v1+v2
        return tot
    
    def subFourVector(self,obj1,obj2):
        tot = ROOT.TLorentzVector()
        v1 = ROOT.TLorentzVector()
        v2 = ROOT.TLorentzVector()
        v1.SetPtEtaPhiM(obj1.pt, 0, obj1.phi, 0)
        v2.SetPtEtaPhiM(obj2.pt, 0, obj2.phi, 0)
        tot = v1-v2
        return tot
    
    def analyze(self, event):
        jets      = Collection(event, "Jet")
        genjets   = Collection(event, "GenJet")
        met       = Object(event,     self.metBranchName)
	weight    = event.genWeight
	eventNum  = event.event

	#Attempt to make new root file
	#Create a new file + a clone of old tree in new file
	inFile, inTree = self.storeSmear("/uscms_data/d3/lpcsusyhad/benwu/Moriond2019/TestNanoAOD/CMSSW_10_4_X_2018-12-11-2300/src/prod2017MC_NANO.root")
	
	#Need to initialize a random seed
	ROOT.gRandom.SetSeed(123456)
        
        #bootstrapping should be done here
        #the histogram can be accessed by doing self.targeth.{some root function to get the value}
        #xBinWidth = float(2/self.targeth.GetNbinsX())
        xBinWidth = 0.01

        #begin smearing
        smearWeight = 1
	SmearJets = []
        eventList = []
	eventList.append(event)
	for iJ in xrange(len(genjets)) :
		if iJ == self.nSmearJets: break
		gJ = genjets[iJ]
		rJI = -1
		if gJ.pt == 0: break
		for iR in xrange(len(jets)) :
			if jets[iR].genJetIdx != iJ:  continue
			rJI = iR
			break

		testMet = 0
		if rJI < 0:
			testMet = self.subFourVector(met, gJ).Pt()
		else:
			testMet = self.testMetCalc(met, jets[rJI], gJ).Pt()
		
		deltamet = testMet - met.pt
		if deltamet > met.pt + 100 and deltamet > 0.55 *gJ.pt: continue
		
		recoJet = jets[rJI]
		if rJI < 0 :
			rJI = len(jets)
			n1=ROOT.TLorentzVector()
			newjet = [n1.SetPtEtaPhiM(9.5,gJ.eta,gJ.phi,gJ.mass),-1,0,9.5,0,gJ.pt]
			recoJet = newjet
		
		origRes_ = self.jetResFunction(recoJet, gJ)
		if origRes_ < 0 or origRes_ > 2 : continue
		
		respHistoName = self.ptmapping(jets[rJI])
		targeth = self.loadHisto(self.respFileName,respHistoName)
		cdf = targeth.GetBinContent(int(origRes_/self.xBinWidth))
		print "CDF", cdf
		minProb, maxProb, minRes, maxRes = self.getScaledWindowAndProb(targeth,origRes_,self.minWindow,self.maxWindow)
		if minProb - maxProb == 0 : continue
		
		SmearJets_buff = [gJ,rJI,targeth,minProb,maxProb,minRes,maxRes] 
		SmearJets.append(SmearJets_buff)

        if len(SmearJets) == 0: return True

        originalRecoJets = jets
        originalMET = met
        originalWeight = weight
	canSmear = False
	SmearedJets = []

        for iS in xrange(self.nSmears) :
		recoJets = []
		originalRecoJets_pt = []
		for iJ in xrange(self.nSmearJets) :
			info = SmearJets[iJ]
			newResValue = 1
			if self.doFlatSampling :
				newResValue = ROOT.gRandom.Uniform(info[5], info[6]) 
			else :
				newResProb = ROOT.gRandom.Uniform(info[3], info[4])   
				newResValue=self.interpolateProbToRes(info[2], newResProb)

			minProb2, maxProb2, minRes2, maxRes2 = self.getContributionScaledWindowAndProb(info[2], newResValue, self.minWindow, self.maxWindow) 
			contribProb = maxProb2 - minProb2
			if contribProb == 0 : continue
			canSmear = True

			smearingCorr = 1
			if self.doFlatSampling:
				deltaMinRes = newResValue - 0.001
				deltaMaxRes = newResValue + 0.001 
				deltaMinProb, deltaMaxProb = self.getWindowProb(info[2], deltaMinRes, deltaMaxRes)
				flatProb = (deltaMaxRes - deltaMinRes)/ (info[6] - info[5])
				trueProb = deltaMaxProb - deltaMinProb
				smearingCorr = trueProb / flatProb
			else:
				smearingCorr = maxProb - minProb

			smearWeight *= smearingCorr / contribProb
			recoJet = jets[info[1]]		
	
			if iJ == 0 : met = self.addFourVector(met, recoJet)
			else:        met = self.addTLorentzVector(met, recoJet)
			newp4 = ROOT.TLorentzVector()
			newp4.SetPtEtaPhiM(newResValue * info[0].pt,recoJet.eta,recoJet.phi,recoJet.mass)
			recoJets.append(newp4.Pt())
			met -= newp4

		for j in xrange(len(originalRecoJets)):
			if j == SmearJets[0][1] or j == SmearJets[1][1] :
				continue
			else :
				recoJets.append(originalRecoJets[j].pt)
			originalRecoJets_pt.append(originalRecoJets[j].pt)

		if canSmear :
			recoJets.sort(key = lambda j : j, reverse = True)
			#print "recoJets: ", recoJets
			smearWeight /= float(self.nSmears)
			weight *= smearWeight
			self.out.fillBranch("Jet_pt",     recoJets)
			self.out.fillBranch("MET_pt",     met.Pt())
			self.storeSmearFile(self.outFileName + str(eventNum) + "_" + str(iS) + ".root", inTree, eventNum)
			#Here is where we need to push values to a new tree

		smearWeight = 1.0
		weight = originalWeight
		canSmear = False
		met = originalMET
		jets = originalRecoJets
		self.out.fillBranch("Jet_pt",     originalRecoJets_pt)
		self.out.fillBranch("MET_pt",     originalMET.pt)
	
        return True    
