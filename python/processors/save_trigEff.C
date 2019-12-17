#include "TFile.h"
#include "TH1.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TString.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TStyle.h"

void save_trigEff() {

	TString year = "2016";
	year = "2017";
	year = "2018";

	TString out_file_name = year + "_trigger_eff.root";
	TFile out_file(out_file_name,"RECREATE");

	//std::vector<TString> histo_name_vec = {"MET_low_dm_QCD", "MET_low_dm_QCD_no_METSig"};
	std::vector<TString> histo_name_vec = {"MET_loose_baseline", "MET_high_dm", "MET_low_dm", "MET_loose_baseline_QCD", "MET_high_dm_QCD", "MET_low_dm_QCD", "MET_low_dm_QCD_METSig", "Electron_pt", "Electron_eta", "Muon_pt", "Muon_eta", "Photon_pt", "Photon_eta", "Zmumu_pt", "Zee_pt"};

	for(int i=0; i < histo_name_vec.size(); i++)
	{
		TString dataset = "SingleElectron";
		TString refTrg = "h_met_all";
		TString sigTrg = "h_met_passtrig";
		TString postfix = "_loose_baseline";
		TString histo_name = histo_name_vec.at(i);
		std::cout<<"\n"<<histo_name<<std::endl;

		if (histo_name == "MET_loose_baseline")
		{
			dataset = "SingleElectron";
			refTrg = "h_met_all_mid";
			sigTrg = "h_met_passtrig_mid";
			postfix = "_loose_baseline";
		}

		if (histo_name == "MET_high_dm")
		{
			dataset = "SingleElectron";
			refTrg = "h_met_all_mid";
			sigTrg = "h_met_passtrig_mid";
			postfix = "_highdm";
		}

		if (histo_name == "MET_low_dm")
		{
			dataset = "SingleElectron";
			refTrg = "h_met_all_mid";
			sigTrg = "h_met_passtrig_mid";
			postfix = "_lowdm";
		}

		if (histo_name == "MET_loose_baseline_QCD")
		{
			dataset = "JetHT_QCD";
			refTrg = "h_met_all_mid";
			sigTrg = "h_met_passtrig_mid";
			postfix = "_loose_baseline";
		}

		if (histo_name == "MET_high_dm_QCD")
		{
			dataset = "JetHT_QCD";
			refTrg = "h_met_all_mid";
			sigTrg = "h_met_passtrig_mid";
			postfix = "_highdm";
		}

		if (histo_name == "MET_low_dm_QCD_METSig")
		{
			dataset = "JetHT_QCD";
			refTrg = "h_met_all_mid";
			sigTrg = "h_met_passtrig_mid";
			postfix = "_lowdm";
		}

		if (histo_name == "MET_low_dm_QCD")
		{
			dataset = "JetHT_QCD_no_METSig";
			refTrg = "h_met_all_mid";
			sigTrg = "h_met_passtrig_mid";
			postfix = "_lowdm";
		}

		if (histo_name == "Electron_pt")
		{
			dataset = "MET";
			refTrg = "h_ele_all_mid";
			sigTrg = "h_ele_passtrig_mid";
			postfix = "_loose_baseline";
		}

		if (histo_name == "Electron_eta")
		{
			dataset = "MET";
			refTrg = "h_ele_all_eta_mid";
			sigTrg = "h_ele_passtrig_eta_mid";
			postfix = "_loose_baseline";
		}

		if (histo_name == "Muon_pt")
		{
			dataset = "MET";
			refTrg = "h_mu_all_mid";
			sigTrg = "h_mu_passtrig_mid";
			postfix = "_loose_baseline";
		}

		if (histo_name == "Muon_eta")
		{
			dataset = "MET";
			refTrg = "h_mu_all_eta_mid";
			sigTrg = "h_mu_passtrig_eta_mid";
			postfix = "_loose_baseline";
		}

		if (histo_name == "Photon_pt")
		{
			dataset = "JetHT";
			refTrg = "h_photon_all_mid";
			sigTrg = "h_photon_passtrig_mid";
			postfix = "_loose_baseline";
		}

		if (histo_name == "Photon_eta")
		{
			dataset = "JetHT";
			refTrg = "h_photon_all_eta_mid";
			sigTrg = "h_photon_passtrig_eta_mid";
			postfix = "_loose_baseline";
		}

		if (histo_name == "Zmumu_pt")
		{
			dataset = "MET";
			refTrg = "h_zmumu_all_mid";
			sigTrg = "h_zmumu_passtrig_mid";
			postfix = "_loose_baseline";
		}

		if (histo_name == "Zee_pt")
		{
			dataset = "MET";
			refTrg = "h_zee_all_mid";
			sigTrg = "h_zee_passtrig_mid";
			postfix = "_loose_baseline";
		}

		TString infile = "results" + postfix + "/" + year + "_" + + dataset + ".root";
		std::cout <<infile<< std::endl;
		TFile* f_in = new TFile(infile);
		TH1F* h_met_denom = (TH1F*) f_in->Get("TrigAnalyzerMiniAOD/" + refTrg);
		TH1F* h_met_num = (TH1F*) f_in->Get("TrigAnalyzerMiniAOD/" + sigTrg);

		TEfficiency* h_met_TEff = new TEfficiency(*h_met_num, *h_met_denom);

		h_met_TEff->Draw();
		gPad->Update();
		auto h_temp = h_met_TEff->GetPaintedGraph();
		h_temp->SetName(histo_name);

		TString title = "Ref Trigger: " + dataset;
		h_temp->SetTitle(title);

		out_file.cd();
		h_temp->Write();
		std::cout << typeid(h_temp).name() << std::endl;
	}

	out_file.Close();

	return;
}
