#include "BRawMasterData.h"
#include "BRawMasterHeader.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"
#include "TLine.h"

TH1F* h_nWaveforms = new TH1F("h_nWaveforms","Number of waveforms per event; NoWaveforms [#]; NoE [#]",100,0,100);
TH1F* h_nSamples = new TH1F("h_nSamples","Number of samples per waveform; NoSamples [#]; NoE [#]",1025,0,1025);
TH1F* h_offset = new TH1F("h_offset","Time offset of waveform; Offset [FADC]; NoE [#]",1025,0,1025);

TH1F* h_initPed = new TH1F("h_initPed","Initial pedestal value; Pedestal [FADC channels]; NoE [#]",10000,-50000,50000);
TH1F* h_initPedRMS = new TH1F("h_initPedRMS","Initial RMS pedestal value; Pedestal RMS [FADC channels]; NoE [#]",1000,0,1000);

TCanvas* myCan = new TCanvas("myCan","Results",1024,720);

Int_t fSizeAvgWindowB = 8; // number of samples that are used to calculate pedestal value at the beginning of the waveform
Int_t fSizeAvgWindowE = 5; // number of samples that are used to calculate pedestal value at the end of the waveform
Int_t fNSigmas = 4; // number of sigmas of pedestal that has to be reached to start pulse processing

// Draw the results saved in the global histograms
int DrawResults()
{
	TCanvas* c_nWaveforms = new TCanvas("c_nWaveforms","NWaveforms",800,600);
	h_nWaveforms->Draw();
	TCanvas* c_nSamples = new TCanvas("c_nSamples","NSamples",800,600);
	h_nSamples->Draw();
	TCanvas* c_offset = new TCanvas("c_offset","Offset",800,600);
	h_offset->Draw();
	TCanvas* c_initPed = new TCanvas("c_initPed","InitPed",800,600);
	h_initPed->Draw();
	TCanvas* c_initPedRMS = new TCanvas("c_initPedRMS","InitPedRMS",800,600);
	h_initPedRMS->Draw();

	return 0;
}

int DrawWaveform(Int_t nSamples, Short_t* data, Short_t offset, Double_t initPed)
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
	waveform->Draw("APL");
	// outputFile->cd();
	// waveform->Write();

	TLine* pedestalLine = new TLine(offset,initPed,nSamples+offset,initPed);
    pedestalLine->SetLineColor(kRed);
    pedestalLine->SetLineWidth(3);
    pedestalLine->SetLineStyle(6);
    pedestalLine->Draw();

    myCan->Modified();
    myCan->Update();
    int dummy;
    cin >> dummy;
    delete waveform;
    return dummy;
}

// It calculates average pedestal and pedestalRMS from given bins [firstBin,lastBin)
// Can be used for pedestal estimation at the beginning as well as end of the pulse
void InitialPed(BRawFADCSample* waveforms, Double_t& ped, Double_t& pedrms, Int_t firstBin, Int_t lastBin)
{
  	Short_t* data = waveforms->GetData();

    //initialize desired variables
    ped = 0;
    pedrms = 0;
    Double_t pedSquared = 0;
    Int_t nBins = lastBin-firstBin;

    // pedestal and pedestalRMS calculations
  for(int i = firstBin; i < lastBin; i++) {
    ped += data[i];
    pedSquared += data[i]*data[i];
  }
  pedrms = TMath::Sqrt((nBins*pedSquared - ped*ped)/(nBins*(nBins-1)));
  ped /= nBins;
}

int AnalyzeFADCSample(BRawFADCSample* waveform)
{
	//Set parameters which are waveform shape independent and the same for all pulses in the waveform
    Int_t nSamples = waveform->GetNbins(); // nSamples in the waveform
    Short_t* data = waveform->GetData(); //array of sampled values of waveform
    Short_t offset = waveform->GetOffset(); // offset position in the 1024 FADC time window

    //Calculate initial values for pedestal and pedestalRMS
    Double_t initPed;
    Double_t initPedRMS;
    InitialPed(waveform,initPed,initPedRMS,0,fSizeAvgWindowB);
    h_initPed->Fill(initPed);
    h_initPedRMS->Fill(initPedRMS);
    // h_begEndDiff->Fill(data[0]-data[nbins-1]);
    // if (data[0]-data[nbins-1] > 2 && offset == 0)
    if (initPedRMS < 10 && initPedRMS > 5 && offset != 0)
	    DrawWaveform(nSamples,data,offset,initPed);
    	// DrawWaveform(data,nbins,initPed,offset,0,0,0,0,0,0);

    // if (nSamples > 150 and nSamples < 1000)
    // if (nSamples > 1000)
	    // DrawWaveform(nSamples,data,offset);

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
        	if (waveform->GetNbins() < 15)
        		continue;
            AnalyzeFADCSample(nsdmaster->GetFADCSample(j));
        }
    }

    DrawResults();

	return 0;
}