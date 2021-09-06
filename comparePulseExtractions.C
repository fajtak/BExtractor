#include <iostream>

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TRatioPlot.h"
#include "TROOT.h"
#include "TStyle.h"

#include "BExtractedImpulseTel.h"
#include "BExtractedHeader.h"

TH1F* h_nHitsOld = new TH1F("h_nHitsOld","old technique;N_{hits} [#]; NoE [#]",100,0,100);
TH1F* h_nHitsNew = new TH1F("h_nHitsNew","new technique;N_{hits} [#]; NoE [#]",100,0,100);
TH1F* h_nHitsDiff = new TH1F("h_nHitsDiff","N_{hits}^{old} - N_{hits}^{new}; N_{hits}^{old} - N_{hits}^{new} [#]; NoE [#]",100,-50,50);
TH1F* h_nHitsDiffPos = new TH1F("h_nHitsDiffPos","N_{hits}^{old} (positive pulses) - N_{hits}^{new}; N_{hits}^{old} - N_{hits}^{new} [#]; NoE [#]",100,-50,50);
TH2F* h_nHitsScatter = new TH2F("h_nHitsScatter","N_{hits}^{old} vs. N_{hits}^{new}; N_{hits}^{old} [#]; N_{hits}^{new} [#]",100,0,100,100,0,100);
TH2F* h_nHitsScatterPos = new TH2F("h_nHitsScatterPos","N_{hits}^{old} (positive pulses) vs. N_{hits}^{new}; N_{hits}^{old} [#]; N_{hits}^{new} [#]",100,0,100,100,0,100);

TH1F* h_chargeOld = new TH1F("h_chargeOld","old technique;Q [FADC channels]",1000,0,10000);
TH1F* h_chargeNew = new TH1F("h_chargeNew","new technique;Q [FADC channels]",1000,0,10000);
TH1F* h_ampOld = new TH1F("h_ampOld","old technique; A [FADC channels]",1000,0,2000);
TH1F* h_ampNew = new TH1F("h_ampNew","new technique; A [FADC channels]",1000,0,2000);
TH1F* h_timeOld = new TH1F("h_timeOld","old technique; T [FADC channels]",1024,0,1024);
TH1F* h_timeNew = new TH1F("h_timeNew","new technique; T [FADC channels]",1024,0,1024);
TH1F* h_NchOld = new TH1F("h_NchOld","old technique; OM ID [#}",300,0,300);
TH1F* h_NchNew = new TH1F("h_NchNew","new technique; OM ID [#]",300,0,300);
TH1F* h_noiseChargeOld = new TH1F("h_noiseChargeOld","old technique;Q [FADC channels]; NoE [#]",1000,0,1000);
TH1F* h_noiseChargeNew = new TH1F("h_noiseChargeNew","new technique;Q [FADC channels]; NoE [#]",1000,0,1000);

int comparePulseExtractions()
{
	TChain* oldTree = new TChain("Events");
	oldTree->Add("/home/fajtak/work/data/baikalData/exp16_barsv051/cluster0/0070/f0070.extr.events.*.sorted.old.root");

	// TFile* oldFile = new TFile("/home/fajtak/work/data/baikalData/exp16_barsv051/cluster0/0070/f0070.extr.events.212.sorted.old.root");
	// TTree* oldTree = (TTree*)oldFile->Get("Events");


	BExtractedImpulseTel* impulseTelOld = NULL;
	BExtractedHeader* headerOld = NULL;
	oldTree->SetBranchAddress("BExtractedImpulseTel.",&impulseTelOld);
	oldTree->SetBranchAddress("BExtractedHeader.",&headerOld);

	TChain* newTree = new TChain("Events");
	newTree->Add("/home/fajtak/work/data/baikalData/exp16_barsv051/cluster0/0070/f0070.extr.events.*.sorted.dp.root");

	// TFile* newFile = new TFile("/home/fajtak/work/data/baikalData/exp16_barsv051/cluster0/0070/f0070.extr.events.212.sorted.dp.root");
	// TTree* newTree = (TTree*)newFile->Get("Events");

	BExtractedImpulseTel* impulseTelNew = NULL;
	BExtractedHeader* headerNew = NULL;
	newTree->SetBranchAddress("BExtractedImpulseTel.",&impulseTelNew);
	newTree->SetBranchAddress("BExtractedHeader.",&headerNew);

	int dummy;

	int nHitsOld = 0;
	int nHitsPosOld = 0;
	int nHitsNew = 0;
	int nHitsPosNew = 0;

	cout << "Old Entries: " << oldTree->GetEntries() << endl;
	cout << "New Entries: " << newTree->GetEntries() << endl;

	for (int i = 0; i < oldTree->GetEntries(); ++i)
	{
		if (i%(oldTree->GetEntries()/10) == 0)
		{
			cout << round((double)(i)/oldTree->GetEntries()*100) << "% ";
			cout << std::flush;
		}

		oldTree->GetEntry(i);
		newTree->GetEntry(i);

		nHitsOld += impulseTelOld->GetNimpulse();
		nHitsPosOld += impulseTelOld->GetNpos();
		nHitsNew += impulseTelNew->GetNimpulse();
		nHitsPosNew += impulseTelNew->GetNpos();

		h_nHitsOld->Fill(impulseTelOld->GetNimpulse());
		h_nHitsNew->Fill(impulseTelNew->GetNimpulse());
		h_nHitsDiff->Fill(impulseTelOld->GetNimpulse()-impulseTelNew->GetNimpulse());
		h_nHitsDiffPos->Fill(impulseTelOld->GetNpos()-impulseTelNew->GetNimpulse());
		h_nHitsScatter->Fill(impulseTelOld->GetNimpulse(),impulseTelNew->GetNimpulse());
		h_nHitsScatterPos->Fill(impulseTelOld->GetNpos(),impulseTelNew->GetNimpulse());

		if(impulseTelOld->GetNimpulse()-impulseTelNew->GetNimpulse() < -5)
		{
			// impulseTelOld->Print();
			// impulseTelNew->Print();
			// std::cin >> dummy;
		}

		for (int j = 0; j < impulseTelOld->GetNimpulse(); ++j)
		{
			h_chargeOld->Fill(impulseTelOld->GetQ(j));
			h_ampOld->Fill(impulseTelOld->GetA(j));
			h_timeOld->Fill(impulseTelOld->GetT(j));
			h_NchOld->Fill(impulseTelOld->GetNch(j));
			if (impulseTelOld->GetT(j) < 400)
				h_noiseChargeOld->Fill(impulseTelOld->GetQ(j));
		}
		for (int j = 0; j < impulseTelNew->GetNimpulse(); ++j)
		{
			h_chargeNew->Fill(impulseTelNew->GetQ(j));
			h_ampNew->Fill(impulseTelNew->GetA(j));
			h_timeNew->Fill(impulseTelNew->GetT(j));
			h_NchNew->Fill(impulseTelNew->GetNch(j));
			if (impulseTelNew->GetT(j) < 400)
				h_noiseChargeNew->Fill(impulseTelNew->GetQ(j));
		}
	}

	std::cout << "NHitsOld: " << nHitsOld << std::endl;
	std::cout << "NHitsOldPos: " << nHitsPosOld << std::endl;
	std::cout << "NHitsNew: " << nHitsNew << std::endl;
	std::cout << "NHitsNewPos: " << nHitsPosNew << std::endl;

	TCanvas* c_nHits = new TCanvas("c_nHits","NHits",800,600);
	gPad->SetGrid();

	h_nHitsOld->Sumw2();
	h_nHitsOld->SetLineColor(kRed);
	h_nHitsNew->Sumw2();
	TRatioPlot* rp_nHits = new TRatioPlot(h_nHitsNew,h_nHitsOld,"pois");
	rp_nHits->Draw("nogrid");
	rp_nHits->GetUpperRefYaxis()->SetTitle("NoE [#]");
	rp_nHits->GetLowerRefYaxis()->SetTitle("N_{hits}^{new}/N_{hits}^{old}");
	rp_nHits->GetUpperPad()->BuildLegend(0.75,0.75,0.95,0.95,"");

	c_nHits->Update();

	TCanvas* c_nHitsDiff = new TCanvas("c_nHitsDiff","NHitsDiff");
	h_nHitsDiff->Draw();
	h_nHitsDiffPos->SetLineColor(kRed);
	h_nHitsDiffPos->Draw("same");

	TCanvas* c_nHitsScatter = new TCanvas("c_nHitsScatter","NHitsScatter");
	h_nHitsScatter->Draw("colz");

	TCanvas* c_nHitsScatterPos = new TCanvas("c_nHitsScatterPos","NHitsScatterPos");
	h_nHitsScatterPos->Draw("colz");

	TCanvas* c_charge = new TCanvas("c_charge","Charge",800,600);
	gPad->SetGrid();

	h_chargeOld->Sumw2();
	h_chargeOld->SetLineColor(kRed);
	h_chargeNew->Sumw2();
	TRatioPlot* rp_charge = new TRatioPlot(h_chargeNew,h_chargeOld,"pois");
	rp_charge->Draw("nogrid");
	rp_charge->GetUpperRefYaxis()->SetTitle("NoE [#]");
	rp_charge->GetLowerRefYaxis()->SetTitle("Q^{new}/Q^{old}");
	rp_charge->GetUpperPad()->BuildLegend(0.75,0.75,0.95,0.95,"");

	TCanvas* c_amp = new TCanvas("c_amp","Amplitude",800,600);
	gPad->SetGrid();

	h_ampOld->Sumw2();
	h_ampOld->SetLineColor(kRed);
	h_ampNew->Sumw2();
	TRatioPlot* rp_amplitude = new TRatioPlot(h_ampNew,h_ampOld,"pois");
	rp_amplitude->Draw("nogrid");
	rp_amplitude->GetUpperRefYaxis()->SetTitle("NoE [#]");
	rp_amplitude->GetLowerRefYaxis()->SetTitle("A^{new}/A^{old}");
	rp_amplitude->GetUpperPad()->BuildLegend(0.75,0.75,0.95,0.95,"");

	TCanvas* c_time = new TCanvas("c_time","Time",800,600);
	gPad->SetGrid();

	h_timeOld->Sumw2();
	h_timeOld->SetLineColor(kRed);
	h_timeNew->Sumw2();
	TRatioPlot* rp_time = new TRatioPlot(h_timeNew,h_timeOld,"pois");
	rp_time->Draw("nogrid");
	rp_time->GetUpperRefYaxis()->SetTitle("NoE [#]");
	rp_time->GetLowerRefYaxis()->SetTitle("T^{new}/T^{old}");
	rp_time->GetUpperPad()->BuildLegend(0.75,0.75,0.95,0.95,"");

	TCanvas* c_Nch = new TCanvas("c_Nch","Nch",800,600);
	gPad->SetGrid();

	h_NchOld->Sumw2();
	h_NchOld->SetLineColor(kRed);
	h_NchNew->Sumw2();
	TRatioPlot* rp_Nch = new TRatioPlot(h_NchNew,h_NchOld,"pois");
	rp_Nch->Draw("nogrid");
	rp_Nch->GetUpperRefYaxis()->SetTitle("NoE [#]");
	rp_Nch->GetLowerRefYaxis()->SetTitle("Nch^{new}/Nch^{old}");
	rp_Nch->GetUpperPad()->BuildLegend(0.75,0.75,0.95,0.95,"");

	TCanvas* c_noiseCharge = new TCanvas("c_noiseCharge","NoiseCharge",800,600);
	h_noiseChargeOld->Draw();
	h_noiseChargeOld->SetLineColor(kRed);
	h_noiseChargeNew->Draw("same");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");


	return 0;
}