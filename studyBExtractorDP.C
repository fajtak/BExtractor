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

// TH1F* h_nHitsOld = new TH1F("h_nHitsOld","old technique;N_{hits} [#]; NoE [#]",100,0,100);
// TH1F* h_nHitsNew = new TH1F("h_nHitsNew","new technique;N_{hits} [#]; NoE [#]",100,0,100);
// TH1F* h_nHitsDiff = new TH1F("h_nHitsDiff","N_{hits}^{old} - N_{hits}^{new}; N_{hits}^{old} - N_{hits}^{new} [#]; NoE [#]",100,-50,50);
// TH1F* h_nHitsDiffPos = new TH1F("h_nHitsDiffPos","N_{hits}^{old} (positive pulses) - N_{hits}^{new}; N_{hits}^{old} - N_{hits}^{new} [#]; NoE [#]",100,-50,50);
// TH2F* h_nHitsScatter = new TH2F("h_nHitsScatter","N_{hits}^{old} vs. N_{hits}^{new}; N_{hits}^{old} [#]; N_{hits}^{new} [#]",100,0,100,100,0,100);
// TH2F* h_nHitsScatterPos = new TH2F("h_nHitsScatterPos","N_{hits}^{old} (positive pulses) vs. N_{hits}^{new}; N_{hits}^{old} [#]; N_{hits}^{new} [#]",100,0,100,100,0,100);

// TH1F* h_chargeOld = new TH1F("h_chargeOld","old technique;Q [FADC channels]",1000,0,10000);
// TH1F* h_chargeNew = new TH1F("h_chargeNew","new technique;Q [FADC channels]",1000,0,10000);
// TH1F* h_ampOld = new TH1F("h_ampOld","old technique; A [FADC channels]",1000,0,2000);
// TH1F* h_ampNew = new TH1F("h_ampNew","new technique; A [FADC channels]",1000,0,2000);
// TH1F* h_timeOld = new TH1F("h_timeOld","old technique; T [FADC channels]",1024,0,1024);
// TH1F* h_timeNew = new TH1F("h_timeNew","new technique; T [FADC channels]",1024,0,1024);
// TH1F* h_NchOld = new TH1F("h_NchOld","old technique; OM ID [#}",300,0,300);
// TH1F* h_NchNew = new TH1F("h_NchNew","new technique; OM ID [#]",300,0,300);
// TH1F* h_noiseChargeOld = new TH1F("h_noiseChargeOld","old technique;Q [FADC channels]; NoE [#]",1000,0,1000);
// TH1F* h_noiseChargeNew = new TH1F("h_noiseChargeNew","new technique;Q [FADC channels]; NoE [#]",1000,0,1000);

TH2F* h_amplitudeChi2 = new TH2F("h_amplitudeChi2","Amplitude vs. Chi^{2};Amplitude [FADC channels]; chi^{2}",2000,0,2000,1000,0,100000);
TH2F* h_amplitudeTOT = new TH2F("h_amplitudeTOT","Amplitude vs. TOT;Amplitude [FADC channels]; TOT [FADC channels]",2000,0,2000,1000,0,1000);
TH2F* h_amplitudeQA = new TH2F("h_amplitudeQA","Amplitude vs. QA Ratio;Amplitude [FADC channels]; Q/A [FADC channels/FADC channels]",2000,0,2000,1000,0,1000);
TH2F* h_amplitudeTFWHM = new TH2F("h_amplitudeTFWHM","Amplitude vs. TFWHM;Amplitude [FADC channels]; TFWHM [FADC channels]",2000,0,2000,1000,0,1000);

TH2F* h_amplitudeScatter = new TH2F("h_amplitudeScatter","A vs. GumpA for SPs; A [FADC channels]; GumpA [FADC channels]",300,0,3000,300,0,3000);
TH2F* h_chargeScatter = new TH2F("h_chargeScatter","Q vs. GumpQ for SPs; Q [FADC channels]; GumpQ [FADC channels]",1000,0,100000,1000,0,100000);

TH1F* h_timeDiff = new TH1F("h_timeDiff", "Time - TimeGump; T-T_{Gump} [ns]; NoE [#]",200,-100,100);

int studyBExtractorDP()
{
	TChain* tree = new TChain("Events");
	// tree->Add("/home/fajtak/work/data/baikalData/exp16_barsv051/cluster0/0070/f0070.extr.events.*.sorted.dp.root");
	// tree->Add("/home/fajtak/work/data/baikalData/exp16_barsv051/cluster0/0070/f0070.extr.events.*.sorted.root");
	tree->Add("/home/fajtak/work/data/baikalData/exp16_barsv051/cluster0/0070/f0070.extr.events.*.sorted.new.root");

	BExtractedImpulseTel* impulseTel = NULL;
	BExtractedHeader* header = NULL;
	tree->SetBranchAddress("BExtractedImpulseTel.",&impulseTel);
	tree->SetBranchAddress("BExtractedHeader.",&header);

	BExtractedImpulse* imp = NULL;

	for (int i = 0; i < tree->GetEntries(); ++i)
	{
		if (i%(tree->GetEntries()/10) == 0)
		{
			cout << round((double)(i)/tree->GetEntries()*100) << "% ";
			cout << std::flush;
		}

		tree->GetEntry(i);

		for (int j = 0; j < impulseTel->GetNimpulse(); ++j)
		{
			imp = impulseTel->At(j);
			if (imp->GetIdDP() == 0)
			{
				h_amplitudeChi2->Fill(imp->GetA(),imp->GetChi2()/imp->GetNDF());
				h_amplitudeTOT->Fill(imp->GetA(),imp->GetTOT());
				h_amplitudeQA->Fill(imp->GetA(),imp->GetQ()/imp->GetA());
				h_amplitudeTFWHM->Fill(imp->GetA(),imp->GetTFWHM());
				h_amplitudeScatter->Fill(imp->GetA(),imp->GetAGump());
				h_chargeScatter->Fill(imp->GetQ(),imp->GetQGump());
				// if (!imp->GetIsSat())
					h_timeDiff->Fill((imp->GetT()-imp->GetTGump())*5);
			}
			// h_chargeOld->Fill(impulseTelOld->GetQ(j));
			// h_ampOld->Fill(impulseTelOld->GetA(j));
			// h_timeOld->Fill(impulseTelOld->GetT(j));
			// h_NchOld->Fill(impulseTelOld->GetNch(j));
			// if (impulseTelOld->GetT(j) < 400)
			// 	h_noiseChargeOld->Fill(impulseTelOld->GetQ(j));
		}
	}

	TCanvas* c_amplitudeChi2 = new TCanvas("c_amplitudeChi2","AmplitudeChi2",800,600);
	h_amplitudeChi2->Draw("colz");

	TCanvas* c_amplitudeTOT = new TCanvas("c_amplitudeTOT","AmplitudeTOT",800,600);
	h_amplitudeTOT->Draw("colz");

	TCanvas* c_amplitudeQA = new TCanvas("c_amplitudeQA","AmplitudeQA",800,600);
	h_amplitudeQA->Draw("colz");

	TCanvas* c_amplitudeTFWHM = new TCanvas("c_amplitudeTFWHM","AmplitudeTFWHM",800,600);
	h_amplitudeTFWHM->Draw("colz");

	TCanvas* c_amplitudeScatter = new TCanvas("c_amplitudeScatter","AmplitudeScatter",800,600);
	h_amplitudeScatter->Draw("colz");

	TCanvas* c_chargeScatter = new TCanvas("c_chargeScatter","ChargeScatter",800,600);
	h_chargeScatter->Draw("colz");

	TCanvas* c_timeDiff = new TCanvas("c_timeDiff","TimeDiff",800,600);
	h_timeDiff->Draw();

	return 0;
}