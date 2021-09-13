#include "BRawMasterData.h"
#include "BRawMasterHeader.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraph.h"

TH1F* h_nWaveforms = new TH1F("h_nWaveforms","Number of waveforms per event; NoWaveforms [#]; NoE [#]",100,0,100);
TH1F* h_nSamples = new TH1F("h_nSamples","Number of samples per waveform; NoSamples [#]; NoE [#]",1025,0,1025);
TH1F* h_offset = new TH1F("h_offset","Time offset of waveform; Offset [FADC]; NoE [#]",1025,0,1025);

TCanvas* myCan = new TCanvas("myCan","Results",1024,720);

// Draw the results saved in the global histograms
int DrawResults()
{
	TCanvas* c_nWaveforms = new TCanvas("c_nWaveforms","NWaveforms",800,600);
	h_nWaveforms->Draw();
	TCanvas* c_nSamples = new TCanvas("c_nSamples","NSamples",800,600);
	h_nSamples->Draw();
	TCanvas* c_offset = new TCanvas("c_offset","Offset",800,600);
	h_offset->Draw();

	return 0;
}

int DrawWaveform(Int_t nSamples, Short_t* data, Short_t offset)
{
	myCan->cd();
    TGraph* waveform = new TGraph(nSamples);
    for (int i = 0; i < nSamples; ++i)
    {
        waveform->SetPoint(i,i+offset,data[i]);
    }
	waveform->SetTitle("");
	waveform->GetXaxis()->SetTitle("Time [FADC channels]");
	waveform->GetYaxis()->SetTitle("Amplitude [FADC channels]");
	waveform->Draw("APL*");
	// outputFile->cd();
	// waveform->Write();

    myCan->Modified();
    myCan->Update();
    int dummy;
    cin >> dummy;
    delete waveform;
    return dummy;
}

int AnalyzeFADCSample(BRawFADCSample* waveform)
{
	//Set parameters which are waveform shape independent and the same for all pulses in the waveform
    Int_t nSamples = waveform->GetNbins(); // nSamples in the waveform
    Short_t* data = waveform->GetData(); //array of sampled values of waveform
    Short_t offset = waveform->GetOffset(); // offset position in the 1024 FADC time window

    DrawWaveform(nSamples,data,offset);

	return 0;
}

// The script should mimic the behaviour of the new BExtractorDP
// It means to read raw.events, process the pulses and export them to ext.events
int BExtractorDP(TString fileName)
{
	//open the right file and set the proper variables in the TTree
	TFile *file = new TFile(fileName.Data(), "READ");
	TTree *tree = (TTree*)file->Get("Events");

	BRawMasterData *nsdmaster = 0;
	tree->SetBranchAddress("BRawMasterData.", &nsdmaster);
	BRawMasterHeader *nsdheader = 0;
	tree->SetBranchAddress("BRawMasterHeader.", &nsdheader);

	// loop over all the events
    for (int i = 0; i < tree->GetEntries(); ++i)
    {
        tree->GetEntry(i);
        h_nWaveforms->Fill(nsdmaster->GetNumSamples());
        // loop over all the waveforms in given events
        for (int j = 0; j < nsdmaster->GetNumSamples(); ++j)
        {
        	BRawFADCSample* waveform = nsdmaster->GetFADCSample(j);
        	h_nSamples->Fill(waveform->GetNbins());
        	h_offset->Fill(waveform->GetOffset());
            AnalyzeFADCSample(nsdmaster->GetFADCSample(j));
        }
    }

    DrawResults();

	return 0;
}