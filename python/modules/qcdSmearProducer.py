#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest

class qcdSmearProducer(Module): 
    def __init__(self)
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
        self.canSmear  =True
        self.respInputName = "JetResByFlav"
        self.respFileName = "file:/eos/uscms/store/user/ddash/qcd_smeared/resTailOut_combined_filtered_CHEF_puWeight_weight_WoH_NORMALIZED.root"
      #  self.respHistoName = self.ptmapping(jets)
      #  self.targeth = self.loadHisto(self.respFileName,self.respHistoName)
        
    def loadHisto(self,filename,hname):
        tf = ROOT.TFile.Open(filename)
        hist = tf.Get(hname)
        hist.SetDirectory(None)
        tf.Close()
        return hist
    
    def ptmapping(self,jets):
        #self.jetspt      = self.analyze(event).jets
        self.ptrange = [jets.pt.range(0,50),jets.pt.range(50,75),jets.pt.range(75,100),jets.pt.range(100,125),jets.pt.range(125,150),jets.pt.range(150,200),jets.pt.range(200,250),jets.pt.range(300,400),jets.pt.range(400,500),jets.pt.range(500,700),jets.pt.range(700,1000),jets.pt.range(1000,1500),jets.pt.range(1500,4000)]
        self.bname=["res_b_comp_14","res_b_comp_15","res_b_comp_16","res_b_comp_17","res_b_comp_18","res_b_comp_19","res_b_comp_20","res_b_comp_21","res_b_comp_22","res_b_comp_23","res_b_comp_24","res_b_comp_25","res_b_comp_26"]
        self.lgtname=["res_light_comp_1","res_light_comp_2","res_light_comp_3","res_light_comp_4","res_light_comp_5","res_light_comp_6","res_light_comp_7","res_light_comp_8","res_light_comp_9","res_light_comp_10","res_light-comp_11","res_light_comp_12","res_light_comp_13"]
        if jetFlavour == 4 :
           bhistmap = map(ptranges,bname)
           return bname
        else :
           lhistmap= map(ptranges,lgtname)
           return lgtname
    def beginJob(self,histFile=None,histDirName=None):
        pass
    def endJob(self):
        pass 

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("origRes", "F");
        self.out.branch("jetFlav", "F");
	self.writeHistFile=True
	self.metBranchName="MET"
	self.xBinWidth = 0.01
	self.minWindow = 0.01
	self.maxWindow = 0.5;
	self.nSmears = 100;
	self.nSmearJets = 2;
	self.nBootstraps = 50;
	self.LINEAR_GRANULATED=True
	self.winType = self.LINEAR_GRANULATED;
	self.doFlatSampling = True;
	self.respInputName = "JetResByFlav";
	self.respFileName = "file:/eos/uscms/store/user/mkilpatr/13TeV/qcd_smearing/resTailOut_combined_filtered_CHEF_puWeight_weight_WoH_NORMALIZED.root"
	self.respHistoName = "res_b_comp_14"
	self.targeth = self.loadHisto(self.respFileName,self.respHistoName)
 
    def loadHisto(self,filename,hname):
	tf = ROOT.TFile.Open(filename)
	hist = tf.Get(hname)
	hist.SetDirectory(None)
	tf.Close()
	return hist

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
        bin = self.cdf.FindFixBin(resp)
        if bin <= 0 :
            return 0
        else :
            binwidth = self.cdf.GetBinWidth(resp)
            b = ((self.cdf.GetBinWidth(bin) + self.cdf.GetBinLowEdge(bin)) * self.cdf.GetBinContent(bin-1) - (self.cdf.GetBinLowEdge(bin)) *self.cdf.GetBinContent(bin)) / bandwidth 
            m = (self.cdf.GetBinContent(bin) - self.cdf.GetBincontent(bin-1))/bandwidth 
            return resp * m + b
    def interpolteProbToRes(self,cdf,probe):
        binAbove = self.cdf.FindFirstBinAbove(probe) 
        if binAbove <= 1 :return 0
        deltaProb = cdf.GetBinContent(binAbove) - cdf.GetBinContent(binAbove -1)
        newResValue = cdf.GEtBinCenter(binAbove)
        if deltaProb > 0 :
           b = (cdf.GetBinContent(binAbove) * cdf.GetBinLowEdge(binAbove)- cdf.GetBinContent(binAbove -1 ) * (cdf.GetBinWidth(binAbove) + cdf.GetBinLowEdge(binAbove)) )/deltaprob
           m = cdf.GetBinWidth(binAbove) / deltaProb
           newResValue = m* probe + b
        return newResValue
     
    def getScaleWindow(self,resp,minW,maxW):
        if resp < 1 :
            return (minW -maxW)*resp + maxW
        else :
            return -1 * (minW -maxW) * resp + 2 * minW - maxW 
        
    def getUpIntegratedscaleWindow(self,resp,minW,maxW):
        if resp < 1 - self.getScaleWindow(1,minW,maxW):
            return (resp + maxW)/(1-(minW - maxW))
        else :
            return (resp + 2*minW - maxW)/(1 +(minW - maxW))
    def getLowIntegratedScaleWindow(self,resp,minW,maxW):
        if resp > 1 + self.getScaleWindow(1,minW,maxW):
            return (resp - ( 2* minW - maxW))/(1- (minW - maxW))
        else :
            return (resp - maxW)/(1 + (minW- maxW))
    def getWindowProb(self,cdf,minProb,maxProb,minRes,maxRes):
        minRes = max(0.0001,minRes)
        maxRes = min(1.9999,maxRes)
        if minRes >= maxRes : 
              minProb=0
              maxProb=0
              return
        else :
          minprob = interpolateresToProbe(cdf,minRes)
          maxprob = interpolateResToProbe(cdf,maxres)

    def getScalWindowAndProb(self,cdf,resp,minWindow,maxWindow,minprob,maxprob,minRes,maxRes):
        window = getScaledWindow(resp,minWindow,maxWindow)
        minRes = resp - window
        maxRes = resp + window
        return self.getWindowProb(cdf,minProb,maxProb,minRes,maxRes)
    def getContriWinAndProb (self,cdf,resp,miniwindow,maxwindow,minProb,maxProb,minRes,maxRes):
        minRes = getLowIntegratedScaledWindow(resp,miniwindow,maxwindow)
        maxRes = getUpIntegratedScaledWindow(resp,miniwindow,maxwindow)
        return self.getWindowProb(cdf,minprob,maxprob,minRes,maxRes)
    def getContributionScaleWindowAndProb(self,cdf,resp,minWindow,maxWindow,mnProb,maxProb,minRes,maxRes):
        minRes = getLowIntegrateScaledWindow(resp,minWindow,maxWindow)
        maxRes = getUpIntegrateScaledWindow(resp,minWindow,maxWindow)
        return self.genWindowProb(cdf,minProb,maxProb,minRes,maxRes)
    
    def addFourVector(self, obj1, obj2):
        tot = ROOT.TLorentzVector()
        v1 = ROOT.TLorentzVector()
        v2 = ROOT.TLorentzVector()
        v1.SetPtEtaPhiM(obj1.pt, 0, obj1.phi, 0)
        v2.SetPtEtaPhiM(obj2.pt, 0, obj2.phi, 0)
        tot = v1+v2
        return tot
    
    def subFourVector(self,obj1,obj2):
        tot = ROOT.TLorentzVector()
        v1 = ROOT.TLorentzVector()
        v2 = ROOT.TLorentzVector()
        v1.SetPtEtaPhiM(obj1.Pt(), 0, obj1.Phi(), 0)
        v2.SetPtEtaPhiM(obj2.pt, 0, obj2.phi, 0)
        tot = v1-v2
        return tot
    def subFourVector1(self,obj1,obj2):
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
        
        #bootstrapping should be done here
        #the histogram can be accessed by doing self.targeth.{some root function to get the value}
        #xBinWidth = float(2/self.targeth.GetNbinsX())
        xBinWidth = 0.01

        #begin smearing
        smearWeight = 1
        nj = 0
        ngjet=0
        for gj in genjets :
            if nj == self.nSmearJets:
                break
            else:
                nj+=1        
            for j in jets :
                rjI = j.genJetIdx
                print "rjI 1st is " ,rjI
                print "ngjet 1st is ",ngjet
                if rjI != ngjet:  continue
                recojet = j
        
                break
        #you know have a matching index to the reco jet
            ngjet+=1
            if rjI < 0 : 
                  #testmet = self.subFourVector1(met , genjets[rjI])
                print "NO MATCH FOUND, I AM SKIPPING ON BECAUSE I DON't KNOW WHAT TO DO"
                continue
            else :
                  testMet1 = self.addFourVector(met, j)
                  if rjI >= len(genjets):
                      print "PANIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                  testMet = self.subFourVector(testMet1,genjets[rjI])
            deltamet = self.subFourVector(testMet,met).Pt()
            print "testmet pt value is : ",testMet.Pt()
            if deltamet > met.pt + 100 and deltamet > 0.55 *genjets[rjI].pt:
                     continue
            if rjI < 0 :
                rjI = len(jets)
                n1=ROOT.TLorentzVector()
                newjet = [n1.SetPtEtaPhiM(9.5,genjets[rjI].eta,genjets[rjI].phi,genjets[rjI].mass),-1,0,9.5,0,genjets[rjI].pt]
            print "nj: ", nj
            print "index: ", rjI
            jetFlavour = j.partonFlavour
            self.out.fillBranch("jetFlav", jetFlavour)
           #This calculates the response with matched gen and reco jets
            origRes_ = self.jetResFunction(j, genjets[rjI])
            if origRes_ < 0 or origRes_ > 0 :
                   continue
            self.respHistoName = self.ptmapping(j)
            self.targeth = self.loadHisto(self.respFileName,self.respHistoName)
            CDF = self.targeth.GetBinContent(int(origRes_/self.xBinWidth))
            print "CDF", CDF
            SmearJets = []
            self.getScalWindowAndProb(cdf,origRes_,minWindow,maxWindow,minProb,maxProb,minRes,maxRes)
            if minProb - maxProb == 0 :
                  continue
            else :
                 SmearJets = [gj,rjI,CDF,minProb,maxProb,minres,maxres]
            self.out.fillBranch("origRes", origRes_)
            if len(SmearJets) == 0:
              pass
        origRecojet = jets
        origmet = met
        origweight = weight
        for js in nSmears :
           for jsj in nSmearJets :
             newResValue = 1
             if doFlatsampling :
              deltaMinRes = newResValue - 0.001
              deltaMaxEes = newResValue + 0.001 
              getWindowProb(cdf,deltaMinProb,deltaMaxProb,deltaMinRes,deltaMaxRes)
              flatProb = (deltaMaxRes - deltaMinRes)/ (maxRes - minRes)
              trueProb = deltaMaxProb - deltaMinProb
              smearingCorr =trueProb / flatProb
              newResValue = random.uniform(minProb,maxProb) 
             else :
              smearingCorr = maxProb - minProb
              newResProb =random.uniform(minProb,maxProb)   
              newResValue=self.interpolateProbToRes(cdf,newResProb)
             self.geContributionScaleWindowAndProb(cdf,newResValue,minWindow,maxWindow,minprob2,maxprob2)
             contriProb = maxProb2 - minProb2
             if contriprob == 0 : continue
             smearWeight *= smearingCorr / contriProb
             met.p4 += recojet.p4
             newp4 = ROOT.TLorentzVector()
             newp4.SetPtEtaPhiM(newResValue * genJets[rjI].pt,recojet.eta,recojet.phi,recojet.mass )
             recojet.SetP4(newp4)
             met.p4 -= recojet.p4
           if(canSmear):
               smearWeight /= self.nSmears
               weight *= smearWeight
           smearWeight =1.0
           weight = origweight
           canSmear = false
           met = origmet
           jets = origRecojet
       return True