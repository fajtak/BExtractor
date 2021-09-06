#include <iostream>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TProfile.h"
#include "TRatioPlot.h"

#include "BExtractedImpulseTel.h"
#include "BExtractedImpulse.h"
#include "BExtractedHeader.h"

TH1F* h_pulseID = new TH1F("h_pulseID", "Pulse ID (0 = SP, 1,2,3,... = order in Multipulse); PulseID [#]; NoE [#]",50,0,50);
TH1F* h_nDPsOMDist = new TH1F("h_nDPsOMDist","DP distribution over OMs;OM ID [#]; Number of DPs [#]",300,0,300);
TH1F* h_nDPsOMDistNoise = new TH1F("h_nDPsOMDistNoise","DP distribution over OMs (fT < 400 FADC);OM ID [#]; Number of DPs [#]",300,0,300);
TH1F* h_nHighDPsOMDist = new TH1F("h_nHighDPsOMDist","Distribution of DPs of higher orders (6 >) over OMs;OM ID [#]; Number of DPs [#]",300,0,300);
TH1F* h_nHitsOMDist = new TH1F("h_nHitsOMDist","Hits distribution over OMs;OM ID [#]; Number of Hits [#]",300,0,300);
TH1F* h_nHitsOMDistNoise = new TH1F("h_nHitsOMDistNoise","Hits distribution over OMs (fT < 400 FADC);OM ID [#]; Number of Hits [#]",300,0,300);

TH1F* h_nHitsRatio = new TH1F("h_nHitsRatio","DP/SP Ratio;OM ID [#]; Number of DPs [#]",300,0,300);
TH1F* h_nHitsNoiseRatio = new TH1F("h_nHitsNoiseRatio","DP/SP Ratio (fT < 400 FADC);OM ID [#]; Number of DPs [#]",300,0,300);

TH1F* h_time = new TH1F("h_time","Time distribution of all pulses;T [FADC channels]; NoE [#]",1024,0,1024);
TH1F* h_timeDP = new TH1F("h_timeDP","Time distribution of DPs;T [FADC channels]; NoE [#]",1024,0,1024);

TH2F* h_timeScatter = new TH2F("h_timeScatter","Time vs. GumpTime for SPs; T [FADC channels]; GumpT [FADC channels]",256,0,1024,256,0,1024);
TH1F* h_timeDiff = new TH1F("h_timeDiff", "Time - TimeGump; T-T_{Gump} [ns]; NoE [#]",200,-100,100);
TH2F* h_chargeScatter = new TH2F("h_chargeScatter","Q vs. GumpQ for SPs; Q [FADC channels]; GumpQ [FADC channels]",1000,0,100000,1000,0,100000);
TH2F* h_amplitudeScatter = new TH2F("h_amplitudeScatter","A vs. GumpA for SPs; A [FADC channels]; GumpA [FADC channels]",300,0,3000,300,0,3000);
TH2F* h_gumpBetaOMDist = new TH2F("h_gumpBetaOMDist","GumpBeta vs. OMID; OM ID [#]; GumpBeta [FADC channels]",300,0,300,500,0,5);
TH2F* h_timeDiffOMDist = new TH2F("h_timeDiffOMDist","Time - GumpTime vs. OMID; OM ID [#]; Time - TimeGump [ns]",300,0,300,2000,-100,100);

TGraph* g_dpRatioPerString = new TGraph(36);
TRatioPlot* rp_nHits;
TRatioPlot* rp_nHitsNoise;
TRatioPlot* rp_time;

void SaveResults(int runID)
{
	TString outputFileName = Form("../results/doublePulses/BExtractor/dp_s2016_c0_r%04d.root",runID);
	TFile* outputFile = new TFile(outputFileName,"RECREATE");

	// h_nDPsOMDistNoise->Divide(h_nHitsOMDistNoise);
	// h_nDPsOMDist->Divide(h_nHitsOMDist);

	h_nHitsRatio->Divide(h_nDPsOMDist,h_nHitsOMDist);
	h_nHitsNoiseRatio->Divide(h_nDPsOMDistNoise,h_nHitsOMDistNoise);
	h_nHitsRatio->Write();
	h_nHitsNoiseRatio->Write();

	h_pulseID->Write();
	h_nDPsOMDist->Write();
	h_nDPsOMDistNoise->Write();
	h_nHighDPsOMDist->Write();
	h_nHitsOMDist->Write();
	h_nHitsOMDistNoise->Write();

	h_time->Write();
	h_timeDP->Write();
	h_timeScatter->Write();
	h_timeDiff->Write();
	h_chargeScatter->Write();
	h_amplitudeScatter->Write();
	h_gumpBetaOMDist->Write();
	h_timeDiffOMDist->Write();
	g_dpRatioPerString->Write();

	// rp_nHits->Write();
	// rp_nHitsNoise->Write();
	// rp_time->Write();
}

void DrawResults()
{
	TCanvas* c_pulseID = new TCanvas("c_pulseID","PulseID",800,600);
	h_pulseID->Draw();

	TCanvas* c_dpOMDist = new TCanvas("c_dpOMDist","DPOMDist",800,600);
	h_nDPsOMDist->Draw();
	h_nHighDPsOMDist->Draw("same");
	h_nHighDPsOMDist->SetLineColor(kGreen);

	TCanvas* c_nHitsOMDist = new TCanvas("c_nHitsOMDist","NHitsOMDist",800,600);
	h_nHitsOMDist->Draw();

	TCanvas* c_dpRatio = new TCanvas("c_dpRatio","DPRatio",800,600);
	gPad->SetGrid();

	// h_nDPsOMDist->Sumw2();
	h_nDPsOMDist->SetLineColor(kRed);
	// h_nHitsOMDist->Sumw2();
	rp_nHits = new TRatioPlot(h_nDPsOMDist,h_nHitsOMDist,"");
	rp_nHits->SetH2DrawOpt("Hist");
	rp_nHits->Draw("nogrid");
	rp_nHits->GetUpperRefYaxis()->SetTitle("NoE [#]");
	rp_nHits->GetUpperRefYaxis()->SetRangeUser(0,h_nHitsOMDist->GetMaximum()*1.1);
	rp_nHits->GetLowerRefYaxis()->SetTitle("N_{DPs}/N^{hits}");
	rp_nHits->GetUpperPad()->BuildLegend(0.75,0.75,0.95,0.95,"");

	TCanvas* c_dpRatioNoise = new TCanvas("c_dpRatioNoise","DPRatioNoise",800,600);
	gPad->SetGrid();

	h_nDPsOMDistNoise->SetLineColor(kRed);
	rp_nHitsNoise = new TRatioPlot(h_nDPsOMDistNoise,h_nHitsOMDistNoise,"");
	rp_nHitsNoise->SetH2DrawOpt("Hist");
	rp_nHitsNoise->Draw("nogrid");
	rp_nHitsNoise->GetUpperRefYaxis()->SetTitle("NoE [#]");
	rp_nHitsNoise->GetUpperRefYaxis()->SetRangeUser(0,h_nHitsOMDistNoise->GetMaximum()*1.1);
	rp_nHitsNoise->GetLowerRefYaxis()->SetTitle("N_{DPs}/N^{hits}");
	rp_nHitsNoise->GetUpperPad()->BuildLegend(0.75,0.75,0.95,0.95,"");

	TGraph* g_nHits = rp_nHits->GetLowerRefGraph();
	for (int i = 0; i < g_nHits->GetN(); ++i)
	{
		cout << i << " " << g_nHits->GetPointY(i) << endl;
		g_dpRatioPerString->SetPoint(i,(((int)g_nHits->GetPointX(i))%36),g_nHits->GetPointY(i)*100);
	}

	TCanvas* c_dpRatioPerString = new TCanvas("c_dpRatioPerString","DPRatioPerString",800,600);
	g_dpRatioPerString->Draw("AP");
	g_dpRatioPerString->SetTitle(";OM ID on the string [#];Percentage of DPs [%]");
	g_dpRatioPerString->SetMarkerStyle(7);

	THStack* s_time = new THStack("s_time","; T [FADC channels];NoE [#]");
	s_time->Add(h_time);
	h_time->SetLineColor(kRed);
	s_time->Add(h_timeDP);

	TCanvas* c_time = new TCanvas("c_time","Time",800,600);
	gPad->SetGrid();
	s_time->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	TCanvas* c_timeRatio = new TCanvas("c_timeRatio","TimeRatio",800,600);
	gPad->SetGrid();

	// h_timeDP->Sumw2();
	// h_timeDP->SetLineColor(kRed);
	// h_time->Sumw2();
	rp_time = new TRatioPlot(h_timeDP,h_time,"pois");
	rp_time->SetH2DrawOpt("Hist");
	rp_time->Draw("nogrid");
	rp_time->GetUpperRefYaxis()->SetTitle("NoE [#]");
	rp_time->GetUpperRefYaxis()->SetRangeUser(0,h_time->GetMaximum()*1.1);
	rp_time->GetLowerRefYaxis()->SetTitle("Ratio");
	rp_time->GetUpperPad()->BuildLegend(0.75,0.75,0.95,0.95,"");


	TCanvas* c_timeScatter = new TCanvas("c_timeScatter","TimeScatter",800,600);
	h_timeScatter->Draw("colz");

	TCanvas* c_timeDiff = new TCanvas("c_timeDiff","TimeDiff",800,600);
	h_timeDiff->Draw();

	TCanvas* c_chargeScatter = new TCanvas("c_chargeScatter","ChargeScatter",800,600);
	h_chargeScatter->Draw("colz");

	TCanvas* c_amplitudeScatter = new TCanvas("c_amplitudeScatter","AmplitudeScatter",800,600);
	h_amplitudeScatter->Draw("colz");

	TCanvas* c_amplitudeProfile = new TCanvas("c_amplitudeProfile","AmplitudeProfile",800,600);
	h_amplitudeScatter->ProfileX()->Draw("");

	TCanvas* c_gumpBetaOMDist = new TCanvas("c_gumpBetaOMDist","GumpBetaOMDist",800,600);
	h_gumpBetaOMDist->Draw("colz");

	TCanvas* c_timeDiffOMDist = new TCanvas("c_timeDiffOMDist","TimeDiffOMDist",800,600);
	h_timeDiffOMDist->Draw("colz");
}

void studyMultiPulses(int runID)
{
	TString storage = "/home/fajtak/work/data/baikalData/exp16_barsv051/cluster0/";

	TChain* tree = new TChain("Events");
	tree->Add(Form("%s/%04d/f%04d.extr.events.*.sorted.dp.root",storage.Data(),runID,runID));
	// tree->Add(v);

	// TFile* file = new TFile("/home/fajtak/work/data/baikalData/exp16_barsv051/cluster0/0070/f0070.extr.events.212.sorted.dp.root");
	// TTree* tree = (TTree*)file->Get("Events");

	BExtractedImpulseTel* impulseTel = NULL;
	BExtractedHeader* header = NULL;
	tree->SetBranchAddress("BExtractedImpulseTel.",&impulseTel);
	tree->SetBranchAddress("BExtractedHeader.",&header);

	BExtractedImpulse* imp = NULL;

	for (int i = 0; i < tree->GetEntries(); ++i)
	{
		if (i%(tree->GetEntries()/10) == 0)
		{
			std::cout << round((double)(i)/tree->GetEntries()*100) << "% ";
			std::cout << std::flush;
		}

		tree->GetEntry(i);

		for (int j = 0; j < impulseTel->GetNimpulse(); ++j)
		{
			imp = impulseTel->At(j);

			if (imp->GetIdDP() == 0)
			{
				// cout << imp->GetT() << " " << imp->GetTGump() << " " << imp->GetBetaGump() << endl;
				h_nHitsOMDist->Fill(imp->GetNch());
				h_time->Fill(imp->GetTGump());
				h_timeScatter->Fill(imp->GetT(),imp->GetTGump());
				h_chargeScatter->Fill(imp->GetQ(),imp->GetQGump());
				h_amplitudeScatter->Fill(imp->GetA(),imp->GetAGump()*0.37);
				h_timeDiff->Fill((imp->GetT()-imp->GetTGump())*5);
				h_gumpBetaOMDist->Fill(imp->GetNch(),imp->GetBetaGump());
				h_timeDiffOMDist->Fill(imp->GetNch(),(imp->GetT()-imp->GetTGump())*5);
				if (imp->GetTGump() < 400)
					h_nHitsOMDistNoise->Fill(imp->GetNch());
			}else{
				;
			}

			h_pulseID->Fill(imp->GetIdDP());

			if (imp->GetIdDP() == 1)
			{
				h_time->Fill(imp->GetTGump());
				h_timeDP->Fill(imp->GetTGump());
				h_nHitsOMDist->Fill(imp->GetNch());
				h_nDPsOMDist->Fill(imp->GetNch());
				if (imp->GetTGump() < 400)
					h_nDPsOMDistNoise->Fill(imp->GetNch());
			}

			if (imp->GetIdDP() == 10)
				h_nHighDPsOMDist->Fill(imp->GetNch());
		}
	}

	DrawResults();
	SaveResults(runID);
}