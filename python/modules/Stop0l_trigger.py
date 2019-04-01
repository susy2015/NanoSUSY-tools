import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class Stop0l_trigger(Module):
    def __init__(self, era):
        self.era = era
	eff_file = "%s/src/PhysicsTools/NanoSUSYTools/data/trigger_eff/" % os.environ['CMSSW_BASE']
	eff_file = eff_file + self.era + "_trigger_eff.root"
	self.tf = ROOT.TFile.Open(eff_file)

    def beginJob(self):
        pass
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Pass_trigger_MET", "O")
        self.out.branch("Pass_trigger_muon", "O")
        self.out.branch("Pass_trigger_electron", "O")
        self.out.branch("Pass_trigger_photon", "O")

	self.out.branch("Stop0l_trigger_eff_MET_loose_baseline", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def mygetattr(self, my_obj, my_branch, default_bool):
	try: getattr(my_obj, my_branch)
	except RuntimeError:
	    #print my_branch, "not found"
	    return default_bool
	else: return getattr(my_obj, my_branch)

    def get_efficiency(self, trigger_name, kinematic):
	eff_hist = self.tf.Get(trigger_name)
	return eff_hist.GetBinContent(eff_hist.FindBinNumber(kinematic))

    def SelPhotons(self, photon):
        #if photon.pt < 200:
        #    return False
        abeta = math.fabs(photon.eta)
        if (abeta > 1.442 and abeta < 1.566) or (abeta > 2.5):
            return False
        ## cut-base ID, 2^0 loose ID
        cutbase =  photon.cutBasedBitmap  if self.era != "2016" else photon.cutBased
        if not cutbase & 0b1:
            return False
        return True

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        hlt       = Object(event, "HLT")
        met       = Object(event, "MET")
	electrons = Collection(event, "Electron")
	muons	  = Collection(event, "Muon")
	photons   = Collection(event, "Photon")

        Pass_trigger_MET = (
	self.mygetattr(hlt, 'PFMET100_PFMHT100_IDTight', False)
	or self.mygetattr(hlt, 'PFMET110_PFMHT110_IDTight', False)
	or self.mygetattr(hlt, 'PFMET120_PFMHT120_IDTight', False)
	or self.mygetattr(hlt, 'PFMET130_PFMHT130_IDTight', False)
	or self.mygetattr(hlt, 'PFMET140_PFMHT140_IDTight', False)
	or self.mygetattr(hlt, 'PFMETNoMu100_PFMHTNoMu100_IDTight', False)
	or self.mygetattr(hlt, 'PFMETNoMu110_PFMHTNoMu110_IDTight', False)
	or self.mygetattr(hlt, 'PFMETNoMu120_PFMHTNoMu120_IDTight', False)
	or self.mygetattr(hlt, 'PFMETNoMu130_PFMHTNoMu130_IDTight', False)
	or self.mygetattr(hlt, 'PFMETNoMu140_PFMHTNoMu140_IDTight', False)
	or self.mygetattr(hlt, 'PFMET100_PFMHT100_IDTight_PFHT60', False)
	or self.mygetattr(hlt, 'PFMET110_PFMHT110_IDTight_PFHT60', False)
	or self.mygetattr(hlt, 'PFMET120_PFMHT120_IDTight_PFHT60', False)
	or self.mygetattr(hlt, 'PFMET130_PFMHT130_IDTight_PFHT60', False)
	or self.mygetattr(hlt, 'PFMET140_PFMHT140_IDTight_PFHT60', False)
	or self.mygetattr(hlt, 'PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60', False)
	or self.mygetattr(hlt, 'PFMETNoMu110_PFMHTNoMu110_IDTight_PFHT60', False)
	or self.mygetattr(hlt, 'PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60', False)
	or self.mygetattr(hlt, 'PFMETNoMu130_PFMHTNoMu130_IDTight_PFHT60', False)
	or self.mygetattr(hlt, 'PFMETNoMu140_PFMHTNoMu140_IDTight_PFHT60', False)
	#or self.mygetattr(hlt, 'PFMET120_PFMHT120_IDTight_HFCleaned', False)
	#or self.mygetattr(hlt, 'PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned', False)
	#or self.mygetattr(hlt, 'PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned', False)
	)

        Pass_trigger_muon = (
	self.mygetattr(hlt, 'IsoMu20', False)
	or self.mygetattr(hlt, 'IsoMu22', False)
	or self.mygetattr(hlt, 'IsoMu24', False)
	or self.mygetattr(hlt, 'IsoMu27', False)
	or self.mygetattr(hlt, 'IsoMu22_eta2p1', False)
	or self.mygetattr(hlt, 'IsoMu24_eta2p1', False)
	or self.mygetattr(hlt, 'IsoTkMu22', False)
	or self.mygetattr(hlt, 'IsoTkMu24', False)
	or self.mygetattr(hlt, 'Mu50', False)
	or self.mygetattr(hlt, 'Mu55', False)
	)

        Pass_trigger_electron = (
	self.mygetattr(hlt, 'Ele105_CaloIdVT_GsfTrkIdT', False)
	or self.mygetattr(hlt, 'Ele115_CaloIdVT_GsfTrkIdT', False)
	or self.mygetattr(hlt, 'Ele135_CaloIdVT_GsfTrkIdT', False)
	or self.mygetattr(hlt, 'Ele145_CaloIdVT_GsfTrkIdT', False)
	or self.mygetattr(hlt, 'Ele25_eta2p1_WPTight_Gsf', False)
	or self.mygetattr(hlt, 'Ele20_eta2p1_WPLoose_Gsf', False)
	or self.mygetattr(hlt, 'Ele27_eta2p1_WPLoose_Gsf', False)
	or self.mygetattr(hlt, 'Ele27_WPTight_Gsf', False)
	or self.mygetattr(hlt, 'Ele35_WPTight_Gsf', False)
	or self.mygetattr(hlt, 'Ele20_WPLoose_Gsf', False)
	or self.mygetattr(hlt, 'Ele45_WPLoose_Gsf', False)
	or self.mygetattr(hlt, 'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL', False)
	or self.mygetattr(hlt, 'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ', False)
	or self.mygetattr(hlt, 'DoubleEle33_CaloIdL_GsfTrkIdVL', False)
	or self.mygetattr(hlt, 'DoubleEle33_CaloIdL_GsfTrkIdVL_MW', False)
	or self.mygetattr(hlt, 'DoubleEle25_CaloIdL_MW', False)
	or self.mygetattr(hlt, 'DoubleEle33_CaloIdL_MW', False)
	)

        Pass_trigger_photon = (
	self.mygetattr(hlt, 'Photon175', False)
	or self.mygetattr(hlt, 'Photon200', False)
	)

	ele_veto = []
	ele_mid = []
	ele_tight = []
	for ele in electrons:
		if (ele.pt > 5 and abs(ele.eta) < 2.5 and ele.cutBasedNoIso >= 1 and ele.miniPFRelIso_all < 0.1):
			ele_veto.append(ele)
			if (ele.cutBasedNoIso >= 3):
				ele_mid.append(ele)
			if (ele.cutBasedNoIso >= 4):
				ele_tight.append(ele)
	n_ele = len(ele_veto)
	n_ele_mid = len(ele_mid)
	n_ele_tight = len(ele_tight)

	zee_mid = []
	if (n_ele_mid == 2 and ele_mid[0].pt > 40 and ele_mid[1].pt > 20 and (ele_mid[0].charge + ele_mid[1].charge) == 0):
		zee_cand = ROOT.TLorentzVector()
		zee_cand = ele_mid[0].p4() + ele_mid[1].p4()
		#Z mass = 91 GeV
		#print "zee_cand mass = ", zee_cand.M()
		if (zee_cand.M() > 81 and zee_cand.M() < 101):
			zee_mid.append(zee_cand)
	n_zee = len(zee_mid)
		
	mu_loose = []
	mu_mid = []
	for mu in muons:
		if (mu.pt > 5 and abs(mu.eta) < 2.4 and mu.miniPFRelIso_all < 0.2):
			mu_loose.append(mu)
			if (mu.mediumId):
				mu_mid.append(mu)
	n_mu = len(mu_loose)
	n_mu_mid = len(mu_mid)

	zmumu_mid = []
	if (n_mu_mid == 2 and mu_mid[0].pt > 50 and mu_mid[1].pt > 20 and (mu_mid[0].charge + mu_mid[1].charge) == 0):
		zmumu_cand = ROOT.TLorentzVector()
		zmumu_cand = mu_mid[0].p4() + mu_mid[1].p4()
		#Z mass = 91 GeV
		#print "zmumu_cand mass = ", zmumu_cand.M()
		if (zmumu_cand.M() > 81 and zmumu_cand.M() < 101):
			zmumu_mid.append(zmumu_cand)
	n_zmumu = len(zmumu_mid)

	photon_loose = []
	for photon in photons:
		if (self.SelPhotons(photon)):
			photon_loose.append(photon)
	n_photon = len(photon_loose)

	MET_trigger_eff_loose_baseline = self.get_efficiency("MET_loose_baseline", met.pt)
	MET_trigger_eff_high_dm = self.get_efficiency("MET_high_dm", met.pt)
	MET_trigger_eff_low_dm = self.get_efficiency("MET_low_dm", met.pt)
	MET_trigger_eff_high_dm_QCD = self.get_efficiency("MET_high_dm_QCD", met.pt)
	MET_trigger_eff_low_dm_QCD = self.get_efficiency("MET_low_dm_QCD", met.pt)
	Electron_trigger_eff_pt = 0
	if (n_ele_mid >=1): Electron_trigger_eff_pt = self.get_efficiency("Electron_pt", ele_mid[0].pt)
	Electron_trigger_eff_eta = 0
	if (n_ele_mid >=1): Electron_trigger_eff_eta = self.get_efficiency("Electron_eta", ele_mid[0].eta)
	Muon_trigger_eff_pt = 0
	if (n_mu_mid >=1): Muon_trigger_eff_pt = self.get_efficiency("Muon_pt", mu_mid[0].pt)
	Muon_trigger_eff_eta = 0
	if (n_mu_mid >=1): Muon_trigger_eff_eta = self.get_efficiency("Muon_eta", mu_mid[0].eta)
	Photon_trigger_eff_pt = 0
	if (n_photon >=1): Photon_trigger_eff_pt = self.get_efficiency("Photon_pt", photon_loose[0].pt)
	Photon_trigger_eff_eta = 0
	if (n_photon >=1): Photon_trigger_eff_eta = self.get_efficiency("Photon_eta", photon_loose[0].eta)
	Zee_trigger_eff_pt = 0
	if (n_zee ==1): Zee_trigger_eff_pt = self.get_efficiency("Zee_pt", zee_mid[0].pt)
	Zmumu_trigger_eff_pt = 0
	if (n_zmumu ==1): Zmumu_trigger_eff_pt = self.get_efficiency("Zmumu_pt", zmumu_mid[0].pt)

        ### Store output
        self.out.fillBranch("Pass_trigger_MET", Pass_trigger_MET)
        self.out.fillBranch("Pass_trigger_muon", Pass_trigger_muon)
        self.out.fillBranch("Pass_trigger_electron", Pass_trigger_electron)
        self.out.fillBranch("Pass_trigger_photon", Pass_trigger_photon)

	self.out.fillBranch("Stop0l_trigger_eff_MET_loose_baseline", MET_trigger_eff_loose_baseline)

	self.tf.Close()
        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
# Stop0lBaseline = lambda : Stop0lBaselineProducer("2016", False)
