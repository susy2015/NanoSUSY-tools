#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest

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
      #  self.respHistoName = sself.ptmapping(jets)
      #  self.targeth = self.loadHisto(self.respFileName,self.respHistoName)
        
    def loadHisto(self,filename,hname):
        tf = ROOT.TFile.Open(filename)
        hist = tf.Get(hname)
        hist.SetDirectory(None)
        tf.Close()
        return hist
    
    def ptmapping(self,jets):
        #self.jetspt      = self.analyze(event).jets
        #ptrange = [jets.pt.range(0,50),jets.pt.range(50,75),jets.pt.range(75,100),jets.pt.range(100,125),jets.pt.range(125,150),jets.pt.range(150,200),jets.pt.range(200,250),jets.pt.range(300,400),jets.pt.range(400,500),jets.pt.range(500,700),jets.pt.range(700,1000),jets.pt.range(1000,1500),jets.pt.range(1500,4000)]

	ptrange = [0, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 700, 1000, 1500]	
        pt_index = -1
	for i in xrange(len(ptrange)):
		if jets.pt < float(ptrange[i]): 
			pt_index = i - 1
			break
	
	bname=["res_b_comp_14","res_b_comp_15","res_b_comp_16","res_b_comp_17","res_b_comp_18","res_b_comp_19","res_b_comp_20","res_b_comp_21","res_b_comp_22","res_b_comp_23","res_b_comp_24","res_b_comp_25","res_b_comp_26"]
        lgtname=["res_light_comp_1","res_light_comp_2","res_light_comp_3","res_light_comp_4","res_light_comp_5","res_light_comp_6","res_light_comp_7","res_light_comp_8","res_light_comp_9","res_light_comp_10","res_light-comp_11","res_light_comp_12","res_light_comp_13"]
	print "index: ", pt_index
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
        self.out.branch("origRes", "F");
        self.out.branch("jetFlav", "F");

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def jetResFunction(self, jets, genjets):
        res = jets.pt/genjets.pt
        return res
#    def pushToTree()

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
        weight    = Object(event,     "genWeight")
 
# matching gen jet can be called by the index Jet_genJetIdx, jet.genJetIdx == matched GenJet
	ROOT.gRandom.SetSeed(123456)
        
        #bootstrapping should be done here
        #the histogram can be accessed by doing self.targeth.{some root function to get the value}
        #xBinWidth = float(2/self.targeth.GetNbinsX())
        xBinWidth = 0.01

        #begin smearing
        smearWeight = 1
	SmearJets = []
        for iJ in range(0, len(genjets)) :
		if iJ == self.nSmearJets: break
		gJ = genjets[iJ]
		rJI = -1
		if gJ.pt == 0: break
		for iR in range(0, len(jets)) :
			if jets[iR].genJetIdx != iJ:  continue
			rJI = iR
			print "recoJets: ", iR
			break

        #you know have a matching index to the reco jet
		testMet = 0
		if rJI < 0:
			testMet = self.subFourVector(met, gJ).Pt()
		else:
			testMet = self.testMetCalc(met, jets[rJI], gJ).Pt()
		
		deltamet = testMet - met.pt
		print "testmet pt value is : ",testMet
		if deltamet > met.pt + 100 and deltamet > 0.55 *gJ.pt: continue
		
		if rJI < 0 :
			rJI = len(jets)
			n1=ROOT.TLorentzVector()
			newjet = [n1.SetPtEtaPhiM(9.5,gJ.eta,gJ.phi,gJ.mass),-1,0,9.5,0,gJ.pt]
		
		print "iJ: ", iJ
		print "index: ", rJI
		
		origRes_ = self.jetResFunction(jets[rJI], gJ)
		print "origRes_: ", origRes_
		if origRes_ < 0 or origRes_ > 2 : continue
		
		respHistoName = self.ptmapping(jets[rJI])
		targeth = self.loadHisto(self.respFileName,respHistoName)
		cdf = targeth.GetBinContent(int(origRes_/self.xBinWidth))
		print "CDF", cdf
		minProb, maxProb, minRes, maxRes = self.getScaledWindowAndProb(targeth,origRes_,self.minWindow,self.maxWindow)
		if minProb - maxProb == 0 : continue
		
		SmearJets_buff = [gJ,rJI,targeth,minProb,maxProb,minRes,maxRes] #How is this storing multiple jets to smear
		SmearJets.append(SmearJets_buff)
		self.out.fillBranch("origRes", origRes_)

        if len(SmearJets) == 0: return True

        originalRecoJets = jets
        originalMET = met
        originalWeight = weight
	canSmear = False

	print "SmearJets Size: ", len(SmearJets)
	print "SmearJets: ", SmearJets

        for iS in range(0, self.nSmears) :
		for iJ in range(0, self.nSmearJets) :
			info = SmearJets[iJ]
			newResValue = 1
			if self.doFlatSampling :
				newResValue = ROOT.gRandom.Uniform(info[5], info[6]) 
			else :
				newResProb = ROOT.gRandom.Uniform(info[3], info[4])   
				newResValue=self.interpolateProbToRes(info[2], newResProb)

			minProb2, maxProb2, minRes2, maxRes2 = self.getContributionScaledWindowAndProb(info[2], newResValue, self.minWindow, self.maxWindow) #how is this getting the right values
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
			#this can't be working because met/jets aren't lorentz vectors
			recoJet = jets[info[1]]		
	
			if iJ == 0 : met = self.addFourVector(met, recoJet)
			else:        met = self.addTLorentzVector(met, recoJet)
			newp4 = ROOT.TLorentzVector()
			newp4.SetPtEtaPhiM(newResValue * info[0].pt,recoJet.eta,recoJet.phi,recoJet.mass)
			recoJet = newp4
			met -= recoJet
		if(canSmear):
			smearWeight /= float(self.nSmears)
			#weight *= smearWeight
			#Here is where we need to push values to a new tree

		smearWeight = 1.0
		weight = originalWeight
		canSmear = False
		met = originalMET
		jets = originalRecoJets

        return True    
