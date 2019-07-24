#include "TFile.h"
#include "TH1.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TString.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TStyle.h"

void plot_trigEff() {

	TString year = "2016";
	//year = "2017";
	//year = "2018";
	TString postfix = "_loose_baseline";
	//postfix = "_highdm";
	//postfix = "_lowdm";

	bool plot_MET = false;
	bool plot_MET_mid = false;
	bool plot_MET_QCD = false;
	bool plot_MET_QCD_compare = false;
	bool plot_MET_mid_compare = false;
	bool plot_Ele = false;
	bool plot_Ele_eta = false;
	bool plot_Ele_mid = false;
	bool plot_Ele_eta_mid = false;
	bool plot_Ele_tight = false;
	bool plot_Ele_eta_tight = false;
	bool plot_Mu = false;
	bool plot_Mu_eta = false;
	bool plot_Mu_mid = false;
	bool plot_Mu_eta_mid = false;
	bool plot_photon = false;
	bool plot_photon_eta = false;
	bool plot_photon_mid = false;
	bool plot_photon_eta_mid = true;
	bool plot_zee_mid = false;
	bool plot_zmumu_mid = false;

	TString out_name = "MET";
	bool use_3_dataset = true;
	TString dataset = "SingleElectron";
	TString dataset2 = "SingleMuon";
	TString dataset3 = "JetHT";
	TString refTrg = "h_met_all";
	TString sigTrg = "h_met_passtrig";
	//TString title = "E_{T}^{miss} [GeV]";

	if (plot_MET)
	{
	out_name = "MET";
	use_3_dataset = true;
	dataset = "SingleElectron";
	dataset2 = "SingleMuon";
	dataset3 = "JetHT";
	refTrg = "h_met_all";
	sigTrg = "h_met_passtrig";
	//title = "E_{T}^{miss} [GeV]";
	}

	if (plot_MET_mid)
	{
	out_name = "MET_mid";
	use_3_dataset = true;
	dataset = "SingleElectron";
	dataset2 = "SingleMuon";
	dataset3 = "JetHT";
	//dataset3 = "SinglePhoton";
	refTrg = "h_met_all_mid";
	sigTrg = "h_met_passtrig_mid";
	//title = "E_{T}^{miss} [GeV]";
	}

	if (plot_MET_QCD)
	{
	out_name = "MET_QCD";
	use_3_dataset = false;
	dataset = "JetHT_QCD";
	dataset2 = "SinglePhoton_QCD";
	refTrg = "h_met_all";
	sigTrg = "h_met_passtrig";
	//title = "E_{T}^{miss} [GeV]";
	}

	if (plot_MET_QCD_compare)
	{
	out_name = "MET_QCD_compare";
	use_3_dataset = true;
	dataset = "JetHT_QCD";
	dataset2 = "JetHT_QCD_highdm";
	dataset3 = "JetHT_QCD_lowdm";
	//dataset = "SinglePhoton_QCD";
	//dataset2 = "SinglePhoton_QCD_highdm";
	//dataset3 = "SinglePhoton_QCD_lowdm";
	refTrg = "h_met_all";
	sigTrg = "h_met_passtrig";
	//title = "E_{T}^{miss} [GeV]";
	}

	if (plot_MET_mid_compare)
	{
	out_name = "MET_mid_compare";
	use_3_dataset = true;
	dataset = "SingleMuon";
	dataset2 = "SingleMuon_plus_jets";
	dataset3 = "SingleMuon_plus_b";
	refTrg = "h_met_all_mid";
	sigTrg = "h_met_passtrig_mid";
	//title = "E_{T}^{miss} [GeV]";
	}

	if (plot_Ele)
	{
	out_name = "Ele";
	use_3_dataset = false;
	dataset = "SingleMuon";
	dataset2 = "JetHT";
	refTrg = "h_ele_all";
	sigTrg = "h_ele_passtrig";
	//title = "Electron p_{T} [GeV]";
	}

	if (plot_Ele_eta)
	{
	out_name = "Ele_eta";
	use_3_dataset = false;
	dataset = "SingleMuon";
	dataset2 = "JetHT";
	refTrg = "h_ele_all_eta";
	sigTrg = "h_ele_passtrig_eta";
	//title = "Electron |Eta|, p_{T} > 30 [GeV]";
	}

	if (plot_Ele_mid)
	{
	out_name = "Ele_mid";
	use_3_dataset = true;
	dataset = "MET";
	dataset2 = "JetHT";
	dataset3 = "SingleMuon";
	refTrg = "h_ele_all_mid";
	sigTrg = "h_ele_passtrig_mid";
	//title = "Electron p_{T} [GeV]";
	}

	if (plot_Ele_eta_mid)
	{
	out_name = "Ele_eta_mid";
	use_3_dataset = true;
	dataset = "MET";
	dataset2 = "JetHT";
	dataset3 = "SingleMuon";
	refTrg = "h_ele_all_eta_mid";
	sigTrg = "h_ele_passtrig_eta_mid";
	//title = "Electron |Eta|, p_{T} > 30 [GeV]";
	}

	if (plot_Ele_tight)
	{
	out_name = "Ele_tight";
	use_3_dataset = true;
	dataset = "SingleMuon";
	dataset2 = "JetHT";
	dataset3 = "MET";
	refTrg = "h_ele_all_tight";
	sigTrg = "h_ele_passtrig_tight";
	//title = "Electron p_{T} [GeV]";
	}

	if (plot_Ele_eta_tight)
	{
	out_name = "Ele_eta_tight";
	use_3_dataset = true;
	dataset = "SingleMuon";
	dataset2 = "JetHT";
	dataset3 = "MET";
	refTrg = "h_ele_all_eta_tight";
	sigTrg = "h_ele_passtrig_eta_tight";
	//title = "Electron |Eta|, p_{T} > 30 [GeV]";
	}

	if (plot_Mu)
	{
	out_name = "Mu";
	use_3_dataset = false;
	dataset = "SingleElectron";
	dataset2 = "JetHT";
	refTrg = "h_mu_all";
	sigTrg = "h_mu_passtrig";
	//title = "Electron p_{T} [GeV]";
	}

	if (plot_Mu_eta)
	{
	out_name = "Mu_eta";
	use_3_dataset = false;
	dataset = "SingleElectron";
	dataset2 = "JetHT";
	refTrg = "h_mu_all_eta";
	sigTrg = "h_mu_passtrig_eta";
	//title = "Electron |Eta|, p_{T} > 30 [GeV]";
	}

	if (plot_Mu_mid)
	{
	out_name = "Mu_mid";
	use_3_dataset = true;
	dataset = "MET";
	dataset2 = "JetHT";
	dataset3 = "SingleElectron";
	refTrg = "h_mu_all_mid";
	sigTrg = "h_mu_passtrig_mid";
	//title = "Electron p_{T} [GeV]";
	}

	if (plot_Mu_eta_mid)
	{
	out_name = "Mu_eta_mid";
	use_3_dataset = true;
	dataset = "MET";
	dataset2 = "JetHT";
	dataset3 = "SingleElectron";
	refTrg = "h_mu_all_eta_mid";
	sigTrg = "h_mu_passtrig_eta_mid";
	//title = "Electron p_{T} [GeV]";
	}

	if (plot_photon)
	{
	out_name = "photon";
	use_3_dataset = true;
	dataset = "SingleMuon";
	dataset2 = "JetHT";
	dataset3 = "MET";
	refTrg = "h_photon_all";
	sigTrg = "h_photon_passtrig";
	//title = "E_{T}^{miss} [GeV]";
	}

	if (plot_photon_eta)
	{
	out_name = "photon_eta";
	use_3_dataset = true;
	dataset = "SingleElectron";
	dataset2 = "SingleMuon";
	dataset3 = "JetHT";
	refTrg = "h_photon_all_eta";
	sigTrg = "h_photon_passtrig_eta";
	//title = "E_{T}^{miss} [GeV]";
	}

	if (plot_photon_mid)
	{
	out_name = "photon_mid";
	use_3_dataset = true;
	dataset = "MET";
	dataset2 = "JetHT";
	dataset3 = "SingleMuon";
	refTrg = "h_photon_all_mid";
	sigTrg = "h_photon_passtrig_mid";
	//title = "E_{T}^{miss} [GeV]";
	}

	if (plot_photon_eta_mid)
	{
	out_name = "photon_eta_mid";
	use_3_dataset = true;
	dataset = "MET";
	dataset2 = "JetHT";
	dataset3 = "SingleMuon";
	refTrg = "h_photon_all_eta_mid";
	sigTrg = "h_photon_passtrig_eta_mid";
	//title = "E_{T}^{miss} [GeV]";
	}

	if (plot_zee_mid)
	{
	out_name = "zee_mid";
	use_3_dataset = true;
	dataset = "MET";
	dataset2 = "JetHT";
	dataset3 = "SingleMuon";
	refTrg = "h_zee_all_mid";
	sigTrg = "h_zee_passtrig_mid";
	//title = "Electron p_{T} [GeV]";
	}

	if (plot_zmumu_mid)
	{
	out_name = "zmumu_mid";
	use_3_dataset = true;
	dataset = "MET";
	dataset2 = "JetHT";
	dataset3 = "SingleElectron";
	refTrg = "h_zmumu_all_mid";
	sigTrg = "h_zmumu_passtrig_mid";
	//title = "Electron p_{T} [GeV]";
	}

	//gStyle->SetPadTopMargin(0.08);
	gStyle->SetPadBottomMargin(0.2);
	gStyle->SetPadLeftMargin(0.12);
	//gStyle->SetPadRightMargin(0.05);
	gStyle->SetOptStat(0);
	TH1::SetDefaultSumw2();

	TCanvas* myCanvas = new TCanvas("myCanvas","myCanvas", 600, 600);

	TPad *padup = new TPad("padup", "padup", 0, 0.3, 1, 1.0);
	padup -> SetBottomMargin(0);
	padup -> SetGrid();
	padup -> Draw();
	padup -> cd();

        TLegend* leg = new TLegend(0.3,0,0.7,0.3);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.06);

	TMultiGraph *h_met_MG = new TMultiGraph();

	TString infile = "results" + postfix + "/" + year + "_" + dataset + ".root";
	TFile* f_in = new TFile(infile);
	TH1F* h_met_denom = (TH1F*) f_in->Get("TrigAnalyzerMiniAOD/" + refTrg);
	TH1F* h_met_num = (TH1F*) f_in->Get("TrigAnalyzerMiniAOD/" + sigTrg);

	TH1F* h_met_eff = (TH1F*) h_met_num->Clone();
	h_met_eff->Sumw2();
	h_met_eff->Divide(h_met_denom);
	h_met_eff->SetLineColor(kBlack);
	h_met_eff->SetMarkerStyle(20);
	h_met_eff->SetMarkerColor(kBlack);


	//h_met_eff->Draw("pe");
	//h_met_eff->Draw("axis");

	TEfficiency* h_met_TEff = new TEfficiency(*h_met_num, *h_met_denom);
	h_met_TEff->SetLineColor(kBlack);
	h_met_TEff->SetMarkerStyle(20);
	h_met_TEff->SetMarkerColor(kBlack);
	h_met_TEff->Draw();
	gPad->Update();
	auto h_temp = h_met_TEff->GetPaintedGraph();
	h_met_MG->Add(h_temp);
	
	leg->AddEntry(h_met_TEff,year + " " + dataset,"lep");

	TString infile2 = "results" + postfix + "/" + year + "_" + dataset2 + ".root";
	TFile* f_in2 = new TFile(infile2);
	TH1F* h_met_denom2 = (TH1F*) f_in2->Get("TrigAnalyzerMiniAOD/" + refTrg);
	TH1F* h_met_num2 = (TH1F*) f_in2->Get("TrigAnalyzerMiniAOD/" + sigTrg);

	TH1F* h_met_eff2 = (TH1F*) h_met_num2->Clone();
	h_met_eff2->Sumw2();
	h_met_eff2->Divide(h_met_denom2);
	h_met_eff2->SetLineColor(kBlue);
	h_met_eff2->SetMarkerStyle(21);
	h_met_eff2->SetMarkerColor(kBlue);

	//h_met_eff2->Draw("pe same");

	TEfficiency* h_met_TEff2 = new TEfficiency(*h_met_num2, *h_met_denom2);
	h_met_TEff2->SetLineColor(kBlue);
	h_met_TEff2->SetMarkerStyle(21);
	h_met_TEff2->SetMarkerColor(kBlue);
	h_met_TEff2->Draw();
	gPad->Update();
	auto h_temp2 = h_met_TEff2->GetPaintedGraph();
	h_met_MG->Add(h_temp2);

	leg->AddEntry(h_met_TEff2,year + " " + dataset2,"lep");

	TH1F* h_met_eff3 = NULL;

	if(use_3_dataset)
	{
	TString infile3 = "results" + postfix + "/" + year + "_" + dataset3 + ".root";
	TFile* f_in3 = new TFile(infile3);
	TH1F* h_met_denom3 = (TH1F*) f_in3->Get("TrigAnalyzerMiniAOD/" + refTrg);
	TH1F* h_met_num3 = (TH1F*) f_in3->Get("TrigAnalyzerMiniAOD/" + sigTrg);

	h_met_eff3 = (TH1F*) h_met_num3->Clone();
	h_met_eff3->Sumw2();
	h_met_eff3->Divide(h_met_denom3);
	h_met_eff3->SetLineColor(kRed);
	h_met_eff3->SetMarkerStyle(22);
	h_met_eff3->SetMarkerColor(kRed);

	//h_met_eff3->Draw("pe same");

	TEfficiency* h_met_TEff3 = new TEfficiency(*h_met_num3, *h_met_denom3);
	h_met_TEff3->SetLineColor(kRed);
	h_met_TEff3->SetMarkerStyle(22);
	h_met_TEff3->SetMarkerColor(kRed);
	h_met_TEff3->Draw();
	gPad->Update();
	auto h_temp3 = h_met_TEff3->GetPaintedGraph();
	h_met_MG->Add(h_temp3);

	leg->AddEntry(h_met_TEff3,year + " " + dataset3,"lep");
	}

	h_met_MG->Draw("ape");

	//h_met_MG->GetXaxis()->SetTitle(title);
	//std::cout<< h_met_eff->GetXaxis()->GetXmin() << std::endl;
	//h_met_MG->GetXaxis()->SetRangeUser(h_met_eff->GetXaxis()->GetXmin(), h_met_eff->GetXaxis()->GetXmax());
	h_met_MG->GetXaxis()->SetLimits(h_met_eff->GetXaxis()->GetXmin(), h_met_eff->GetXaxis()->GetXmax());
	h_met_MG->GetYaxis()->SetTitle("Efficiency");
	h_met_MG->GetYaxis()->SetTitleSize(0.05);
	h_met_MG->GetYaxis()->SetLabelSize(0.05);
	h_met_MG->GetYaxis()->SetTitleOffset(0.9);
	h_met_MG->GetYaxis()->SetRangeUser(0,1);

	leg->Draw("same");

	int lumi;
	if(year == "2016") lumi = 36;
	if(year == "2017") lumi = 42;
	if(year == "2018") lumi = 59;

        TLatex latex;
        latex.SetTextSize(0.04);
        latex.SetNDC();
        //latex.SetTextAlign(13);  //align at top
        //latex.DrawLatex(0.5,ymax+0.4,"#bf{CMS} Preliminary, 2017 data");
        latex.DrawLatex(0.12,0.91,"CMS #bf{Preliminary}");
        TString lumi_and_energy = "#bf{" + std::to_string(lumi) + " fb^{-1} (13TeV)}";
        latex.DrawLatex(0.74,0.91,lumi_and_energy);

	myCanvas->cd();
	TPad *paddown = new TPad("paddown", "paddown", 0, 0, 1, 0.3);
	paddown -> SetTopMargin(0);
	paddown -> SetBottomMargin(0.3);
	paddown -> SetGrid();
	paddown -> Draw();
	paddown -> cd();

	TH1F* h_ratio2 = (TH1F*)h_met_eff2->Clone();
	h_ratio2->Divide(h_met_eff);

	h_ratio2->GetXaxis()->SetLabelSize(0.12);
	h_ratio2->GetXaxis()->SetTitleSize(0.12);
	//h_ratio2->GetXaxis()->SetTitleOffset(0.9);
	h_ratio2->GetYaxis()->SetLabelSize(0.08);
	h_ratio2->GetYaxis()->SetTitleSize(0.1);
	h_ratio2->GetYaxis()->SetTitleOffset(0.35);
	h_ratio2->GetYaxis()->SetTitle("Ratio ");
	h_ratio2->GetYaxis()->SetRangeUser(0,2);
	//h_ratio2->GetYaxis()->SetNdivisions(5);

	h_ratio2->Draw();

	if(use_3_dataset)
	{
	TH1F* h_ratio3 = (TH1F*)h_met_eff3->Clone();
	h_ratio3->Divide(h_met_eff);
	h_ratio3->Draw("same");
	}

	myCanvas->SaveAs("plots/" + year + "_HLT_" + out_name + "_Eff" + postfix + ".png"); 

	return;
}
