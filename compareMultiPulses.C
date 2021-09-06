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

int compareMultiPulses()
{
	TFile* run70 = new TFile("../results/doublePulses/BExtractor/dp_s2016_c0_r0070.root","READ");
	TH1F* ratioRun70  = (TH1F*)run70->Get("h_nHitsRatio");
	ratioRun70->SetTitle("50 kHz (run70)");
	ratioRun70->SetLineColor(kRed);
	ratioRun70->SetMarkerColor(kRed);
	ratioRun70->SetMarkerStyle(21);
	TH1F* noiseRatioRun70  = (TH1F*)run70->Get("h_nHitsNoiseRatio");
	noiseRatioRun70->SetTitle("50 kHz (run70)");
	noiseRatioRun70->SetLineColor(kRed);
	noiseRatioRun70->SetMarkerColor(kRed);
	noiseRatioRun70->SetMarkerStyle(21);

	TFile* run140 = new TFile("../results/doublePulses/BExtractor/dp_s2016_c0_r0140.root","READ");
	TH1F* ratioRun140  = (TH1F*)run140->Get("h_nHitsRatio");
	ratioRun140->SetTitle("30 kHz (run140)");
	ratioRun140->SetLineColor(kGreen);
	ratioRun140->SetMarkerColor(kGreen);
	ratioRun140->SetMarkerStyle(22);
	TH1F* noiseRatioRun140  = (TH1F*)run140->Get("h_nHitsNoiseRatio");
	noiseRatioRun140->SetTitle("30 kHz (run140)");
	noiseRatioRun140->SetLineColor(kGreen);
	noiseRatioRun140->SetMarkerColor(kGreen);
	noiseRatioRun140->SetMarkerStyle(22);

	TFile* run271 = new TFile("../results/doublePulses/BExtractor/dp_s2016_c0_r0271.root","READ");
	TH1F* ratioRun271  = (TH1F*)run271->Get("h_nHitsRatio");
	ratioRun271->SetTitle("100 kHz (run271)");
	ratioRun271->SetLineColor(kBlue);
	ratioRun271->SetMarkerColor(kBlue);
	ratioRun271->SetMarkerStyle(23);
	TH1F* noiseRatioRun271  = (TH1F*)run271->Get("h_nHitsNoiseRatio");
	noiseRatioRun271->SetTitle("100 kHz (run271)");
	noiseRatioRun271->SetLineColor(kBlue);
	noiseRatioRun271->SetMarkerColor(kBlue);
	noiseRatioRun271->SetMarkerStyle(23);


	THStack* s_nHitsRatio = new THStack("s_nHitsRatio","All pulses; OM ID [#];N_{DP}/N_{Hits}");
	s_nHitsRatio->Add(ratioRun70,"");
	s_nHitsRatio->Add(ratioRun140,"");
	s_nHitsRatio->Add(ratioRun271,"");

	TCanvas* c_nHitsRatio = new TCanvas("c_nHitsRatio","NHitsRatio",800,600);
	gPad->SetGrid();
	s_nHitsRatio->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	THStack* s_nHitsNoiseRatio = new THStack("s_nHitsNoiseRatio","Only hits from \"noise\" region (fT < 400); OM ID [#];N_{DP}/N_{Hits} (fT < 400)");
	s_nHitsNoiseRatio->Add(noiseRatioRun70,"");
	s_nHitsNoiseRatio->Add(noiseRatioRun140,"");
	s_nHitsNoiseRatio->Add(noiseRatioRun271,"");

	TCanvas* c_nHitsNoiseRatio = new TCanvas("c_nHitsNoiseRatio","NHitsNoiseRatio",800,600);
	gPad->SetGrid();
	s_nHitsNoiseRatio->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	// mcNHits->SetTitle("Muon group MC");
	// mcNHits->SetLineColor(kRed);
	// mcNHits->Rebin(5);
	// mcNHits->Scale(1/(mcDataDuration*24*3600));
	// TH1F* mcNHitsTFilter  = (TH1F*)mcData->Get("h_nHitsAfterTFilter");
	// mcNHitsTFilter->SetTitle("Muon group MC");
	// mcNHitsTFilter->SetLineColor(kRed);
	// mcNHitsTFilter->Rebin(2);
	// mcNHitsTFilter->Scale(1/(mcDataDuration*24*3600));
	// TH1F* mcNStringsTFilter  = (TH1F*)mcData->Get("h_nStringsAfterTFilter");
	// mcNStringsTFilter->SetTitle("Muon group MC");
	// mcNStringsTFilter->SetLineColor(kRed);

	return 0;
}