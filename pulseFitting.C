#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TLegend.h"

#include "BRawMasterData.h"
#include "BRawMasterHeader.h"
#include "BExtractor.h"
#include "BExtractedImpulseTel.h"
#include "BGeomTel.h"
#include "BExtractedImpulse.h"
#include "BExtractedCrossTalkTel.h"
#include "BExtractedCrossTalk.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include <string>
#include <ctime>

const int c_nPed = 5;

class BExtractorMy : public BExtractor
{
public:
	void SetImpulseTel(BExtractedImpulseTel* pExtractedImpulseTel){fImpulses = pExtractedImpulseTel;};
	void SetGeomTel(BGeomTel* pGeomTel){fGeomTel = pGeomTel;};
	void SetRawMasterHeader(BRawMasterHeader* header){fRawMasterHeader = header;};
	void SetCrossTalk(BExtractedCrossTalkTel* crossTalk){fCrossTalks = crossTalk;};
	void Analyze(BRawFADCSample* ch){AnalyzeFADCSample(ch);};

private:
	int i;
};

class BGeomTelMy : public BGeomTel
{
public:
	int GetNchGeom(int a, int b){return 1;};
	BGeomTelMy(){;};
};

class BRawMasterHeaderMy : public BRawMasterHeader
{
public:
	BRawMasterHeaderMy(){;};
	int GetSdc(){return 1;};
};

class BExtractedCrossTalkTelMy: public BExtractedCrossTalkTel
{
public:
	BExtractedCrossTalkTelMy(){;};
	BExtractedCrossTalk* Add(int a){return new BExtractedCrossTalk();};
};

class Pulse
{
public:
	Pulse(){;};
	Pulse(double pedestal, double amplitude, double charge, double fitPed, double fitAmpl, double fitSig, double fitPos, double fitChiSq, int ndf);
	~Pulse();

	double GetAmplitude(){return fAmplitude;};
	double GetCharge(){return fCharge;};
	double GetFitAmplitude(){return fFitAmp;};
	double GetChiSq(){return fFitChiSq;};
	int GetNDF(){return fNDF;};
	double GetChiSqNDF(){return fFitChiSq/fNDF;};
	double GetFitSig(){return fFitSig;}

private:
	double fPedestal;
	double fAmplitude;
	double fCharge;
	double fFitPed;
	double fFitAmp;
	double fFitSig;
	double fFitPos;
	double fFitChiSq;
	int    fNDF;
};

Pulse::Pulse(double pedestal, double amplitude, double charge, double fitPed, double fitAmpl, double fitSig, double fitPos, double fitChiSq, int ndf)
{
	fPedestal = pedestal;
	fAmplitude = amplitude;
	fCharge = charge;
	fFitPed = fitPed;
	fFitAmp = fitAmpl;
	fFitSig = fitSig;
	fFitPos = fitPos;
	fFitChiSq = fitChiSq;
	fNDF = ndf;
}

// Show the dependence of the fitted amplitude and fitted sigma on the charge of the pulse
// plus some additional histograms
int pulseFitting(void)
{
	TString fileName = "/Data/BaikalData/2017/cluster-1/g0020_1.root";

	TFile* file = new TFile(fileName,"READ");

	TTree* tree = (TTree*)file->Get("Events");
	BRawMasterData* rawMasterData = NULL;
    tree->SetBranchAddress("BRawMasterData.", &rawMasterData);
    BRawMasterHeader* rawMasterHeader = NULL;
    tree->SetBranchAddress("BRawMasterHeader.", &rawMasterHeader);

    std::vector<Pulse*> pulses;

    TH1F* h_nbins = new TH1F("h_nbins","Number of bins in waveform",1024,0,1024);

    TF1* fitFunc = new TF1("fitFunc","[0]+[1]*exp(-((x-[2])/[3] + exp(-(x-[2])/[3])))",0,200);

    TCanvas* myCan = new TCanvas("myCan","Results",800,600);

    TGraph* graph = new TGraph();

    for (int i = 0; i < 10000000; ++i)
    {
    	tree->GetEntry(i);
    	// std::cout << rawMasterData->GetNumSamples() << std::endl;
    	for (int j = 0; j < rawMasterData->GetNumSamples(); ++j)
    	{
    		BRawFADCSample* sample = rawMasterData->GetFADCSample(j);
    		h_nbins->Fill(sample->GetNbins());
    		if (sample->GetNch() == 0 && sample->GetNbins() < 90 && sample->GetNbins() != 0 && rawMasterHeader->GetSdc() == 201)
    		{
				Int_t nbins = sample->GetNbins();
				Short_t *data = sample->GetData();
				double pedestal = 0;
				int amplitude = 0;
				int charge = 0;
				int amplitudeBin = 0;
				for (int k = 0; k < c_nPed; ++k)
				{
					pedestal += data[k];
				}
				pedestal /= c_nPed;
				for (int k = 0; k < nbins; ++k)
				{
					if (data[k]-pedestal > amplitude)
					{
						amplitude = data[k]-pedestal;
						amplitudeBin = k;
					}
					charge += data[k]-pedestal;
				}
				// if (amplitude-pedestal > 70	|| amplitude-pedestal < 65)
				if (charge < 0)
					continue;

				graph->Set(nbins);
				for(Int_t n = 0; n < nbins; n++) {
					graph->SetPoint(n, n, data[n]-pedestal);
				}	
				fitFunc->SetParameters(0,amplitude/0.37,amplitudeBin,2.2);
				graph->Fit(fitFunc,"QW");
				// cout << "amplitude = " << amplitude << " charge = " << charge << endl;
				Pulse* currPulse = new Pulse(pedestal,amplitude,charge,fitFunc->GetParameter(0),fitFunc->GetParameter(1),fitFunc->GetParameter(3),fitFunc->GetParameter(2),fitFunc->GetChisquare(),fitFunc->GetNDF());
				pulses.push_back(currPulse);
				// if (pow(amplitude/30,2) + 5 - fitFunc->GetChisquare()/fitFunc->GetNDF() > 0 && fitFunc->GetChisquare()/fitFunc->GetNDF() > 350)
				// {
				// 	cout << "ChiSq: " << fitFunc->GetChisquare() << " " << fitFunc->GetNDF() << " " << pulses.back()->GetChiSqNDF() << " " << pulses.back()->GetAmplitude() << endl;
				// 	graph->Draw();
				// 	myCan->Update();
				// 	cin >> pedestal;
				// 	graph->Clear();
				// }
    		}
    	}
    }

    h_nbins->Draw();
    //graph->Draw();

    TH2F* AvsQ = new TH2F("AvsQ","Amplitude vs. Charge",100,0,500,100,0,2500);
    TH2F* AvsFA = new TH2F("AvsFA","Amplitude vs. Fitted Amplitude",100,0,500,150,0,1500);
    TH2F* QvsFA = new TH2F("QvsFA","Charge of the pulse vs. parameter a;Charge [FADC channels];a [FADC channels]",100,0,2500,150,0,1500);
    TH2F* QvsFS = new TH2F("QvsFS","Charge of the pulse vs. #beta;Charge [FADC channels];#beta [FADC channels]",100,0,2500,100,0,20);
    TH2F* AvsChiSqu = new TH2F("AvsChiSqu","Amplitude vs. Chi2/NDF",200,0,1000,1000,0,2000);
    TH1F* chiSqu = new TH1F("chiSqu","Chi2",1000,0,1000);

    for (unsigned int i = 0; i < pulses.size(); ++i)
    {
    	if (pow(pulses[i]->GetAmplitude()/30,2) + 5 - pulses[i]->GetChiSqNDF() < 0)
    		continue;
    	AvsQ->Fill(pulses[i]->GetAmplitude(),pulses[i]->GetCharge());
    	AvsFA->Fill(pulses[i]->GetAmplitude(),pulses[i]->GetFitAmplitude());
    	QvsFA->Fill(pulses[i]->GetCharge(),pulses[i]->GetFitAmplitude());
    	QvsFS->Fill(pulses[i]->GetCharge(),pulses[i]->GetFitSig());
    	AvsChiSqu->Fill(pulses[i]->GetAmplitude(),pulses[i]->GetChiSqNDF());
    	chiSqu->Fill(pulses[i]->GetChiSqNDF());
    }

    TCanvas* myCan2 = new TCanvas("myCan2","Results",800,600);
    myCan2->Divide(2,3);
    myCan2->cd(1);
    AvsQ->Draw("colz");
    myCan2->cd(2);
    AvsQ->ProfileX()->Draw("colz");
    myCan2->cd(3);
    QvsFA->Draw("colz");
    myCan2->cd(4);
    QvsFA->ProfileX()->Draw("colz");
    myCan2->cd(5);
    AvsFA->Draw("colz");
    myCan2->cd(6);
    AvsFA->ProfileX()->Draw("colz");


    TCanvas* myCan3 = new TCanvas("myCan3","Results",800,600);
    myCan3->Divide(2,2);
    myCan3->cd(1);
    chiSqu->Draw();
    myCan3->cd(2);
    AvsChiSqu->Draw("colz");
    myCan3->cd(3);
    QvsFS->Draw("colz");
    myCan3->cd(4);
    QvsFS->ProfileX()->Draw("colz");

    TCanvas* myCan4 = new TCanvas("myCan4","Results",800,600);
    QvsFA->ProfileX()->Draw("colz");

    TCanvas* myCan5 = new TCanvas("myCan5","Results",800,600);
    QvsFS->ProfileX()->Draw("colz");

    return 0;

}

// Draw a big number of pulses in one TGraph
int pulseFitting(int pulse)
{
	TString fileName = "/Data/BaikalData/2017/cluster-1/g0020_1.root";

	TFile* file = new TFile(fileName,"READ");

	TTree* tree = (TTree*)file->Get("Events");
	BRawMasterData* rawMasterData = NULL;
    tree->SetBranchAddress("BRawMasterData.", &rawMasterData);
    BRawMasterHeader* rawMasterHeader = NULL;
    tree->SetBranchAddress("BRawMasterHeader.", &rawMasterHeader);

    std::vector<Pulse> pulses;

    TH1F* h_nbins = new TH1F("h_nbins","Number of bins in waveform",1024,0,1024);

	TGraph *graph = new TGraph();
	graph->SetName("graph");
	graph->SetLineColor(4);				
	graph->SetMarkerColor(4);				
	graph->SetMarkerStyle(21);				
	graph->SetMarkerSize(0.5);	

    for (int i = 0; i < 10000; ++i)
    {
    	tree->GetEntry(i);
    	// std::cout << rawMasterData->GetNumSamples() << std::endl;
    	for (int j = 0; j < rawMasterData->GetNumSamples(); ++j)
    	{
    		BRawFADCSample* sample = rawMasterData->GetFADCSample(j);
    		h_nbins->Fill(sample->GetNbins());
    		if (sample->GetNch() == 0 && sample->GetNbins() < 90 && sample->GetNbins() != 0 && rawMasterHeader->GetSdc() == 201)
    		{
				Int_t nbins = sample->GetNbins();
				Int_t nfirst = graph->GetN();
				Short_t *data = sample->GetData();
				double pedestal = 0;
				int amplitude = 0;
				int charge = 0;
				for (int k = 0; k < c_nPed; ++k)
				{
					pedestal += data[k];
				}
				pedestal /= c_nPed;
				for (int k = 0; k < nbins; ++k)
				{
					if (data[k]-pedestal > amplitude)
						amplitude = data[k]-pedestal;
					charge += data[k]-pedestal;
				}
				if ((amplitude < 15	|| amplitude > 20) && charge > 0)
				// if (charge < 0)
					continue;
				graph->Set(nfirst + nbins);
				for(Int_t n = 0; n < nbins; n++) {
					graph->SetPoint(nfirst + n, n, data[n]-pedestal);
				}	
				cout << "amplitude = " << amplitude << " charge = " << charge << endl;
				// pulses.push_back(Pulse(pdestal,amplitude))		
    		}
    	}
    }

    h_nbins->Draw();
    graph->Draw();

    return 0;
}

// Draw all the waveforms from the one chosen OM when the test pulse was detected
int pulseFitting(string string)
{
	TString fileName = "/Data/BaikalData/2017/cluster-1/0441/g0441.raw.events.195.sorted.root";

	TFile* file = new TFile(fileName,"READ");

	TTree* tree = (TTree*)file->Get("Events");
	BRawMasterData* rawMasterData = NULL;
    tree->SetBranchAddress("BRawMasterData.", &rawMasterData);
    BRawMasterHeader* rawMasterHeader = NULL;
    tree->SetBranchAddress("BRawMasterHeader.", &rawMasterHeader);

    std::vector<Pulse> pulses;

    TH1F* h_nbins = new TH1F("h_nbins","Number of bins in waveform",1024,0,1024);

	TGraph *graph = new TGraph();
	graph->SetName("graph");
	graph->SetLineColor(4);				
	graph->SetMarkerColor(4);				
	graph->SetMarkerStyle(21);				
	graph->SetMarkerSize(0.5);	

	cout << "Number of events: " << tree->GetEntries() << endl;

	int nEvents = 0;

	TH1F* nSamples = new TH1F("nSamples","Number of samples",50,0,50);
	TH1F* nSamplesFirst = new TH1F("nSamplesFirst","Number of samples on first",50,0,50);
	TH1F* offsetPosition = new TH1F("offsetPosition","Offset position",1024,0,1024);
	TH1F* hNBins = new TH1F("hNBins","Number of bins per sample",200,0,200);
	TH2F* OffsetvsNBins = new TH2F("OffsetvsNBins","Offset vs NBins",1024,0,1024,200,0,200);		

	TCanvas* myCan2 = new TCanvas("myCan2","Results",1200,800);

	int channel = 0;

    for (int i = 0; i < tree->GetEntries(); ++i)
    {
    	// if (i%10 != 0)
    	// 	continue;
    	tree->GetEntry(i);

    	nSamples->Fill(rawMasterData->GetNumSamples());
    	
    	int nSamplesOnFirst = 0;
    	bool testPulse = false;
    	for (int j = 0; j < rawMasterData->GetNumSamples(); ++j)
    	{
    		// cout << "Sample: " << j << " Nch: " << rawMasterData->GetFADCSample(j)->GetNch() << endl;
    		if (rawMasterData->GetFADCSample(j)->GetNch() == channel && rawMasterData->GetFADCSample(j)->GetNbins() != 0)
    		{
    			if (rawMasterData->GetFADCSample(j)->GetOffset() > 706 && rawMasterData->GetFADCSample(j)->GetOffset() < 712)
    				testPulse = true;
    			nSamplesOnFirst++;
    			offsetPosition->Fill(rawMasterData->GetFADCSample(j)->GetOffset());
    			hNBins->Fill(rawMasterData->GetFADCSample(j)->GetNbins());
    			OffsetvsNBins->Fill(rawMasterData->GetFADCSample(j)->GetOffset(),rawMasterData->GetFADCSample(j)->GetNbins());
    		}
    	}
    	if (!testPulse)
    		continue;
    	nSamplesFirst->Fill(nSamplesOnFirst);
    	if (nSamplesOnFirst < 2)
    		continue;

    	//if (!(rawMasterData->GetNumSamples() == 2 && (rawMasterData->GetFADCSample(0)->GetOffset() > 708 && rawMasterData->GetFADCSample(0)->GetOffset() < 710)) )
    		// continue;
    	for (int j = 0; j < rawMasterData->GetNumSamples(); ++j)
    	{
    		BRawFADCSample* sample = rawMasterData->GetFADCSample(j);
    		h_nbins->Fill(sample->GetNbins());
    		if (sample->GetNch() == channel && sample->GetNbins() != 0)
    		{
				Int_t nbins = sample->GetNbins();
				Int_t nfirst = graph->GetN();
				Int_t offset = sample->GetOffset();
				Short_t *data = sample->GetData();
				double pedestal = 0;
				int amplitude = 0;
				int charge = 0;
				for (int k = 0; k < c_nPed; ++k)
				{
					pedestal += data[k];
				}
				pedestal /= c_nPed;
				for (int k = 0; k < nbins; ++k)
				{
					if (data[k]-pedestal > amplitude)
						amplitude = data[k]-pedestal;
					charge += data[k]-pedestal;
				}
				// if ((amplitude < 15	|| amplitude > 20) && charge > 0)
				// if (charge < 0)
					// continue;
				// graph->Set(nfirst + nbins);
				// for(Int_t n = 0; n < nbins; n++) {
				// 	graph->SetPoint(nfirst + n, n+offset, data[n]-pedestal);
				// }	
				// cout << "ID: " << i << " amplitude = " << amplitude << " charge = " << charge << endl;
				// pulses.push_back(Pulse(pdestal,amplitude))	
				nEvents++;	
				
    		}
    	}
  //   	graph->Draw();
  //   	graph->GetXaxis()->SetTitle("Time in CeM FADC window [bins]");
  //   	graph->GetYaxis()->SetTitle("Pulse Amplitude [bins]");
		// myCan2->Update();
		// cin >> nSamplesOnFirst;
		// graph->Clear();
		// graph->Set(0);
    }

    cout << "Number of test pulse & LED pulses: " << nEvents << endl;
    // h_nbins->Draw();
    graph->Draw();

    TCanvas* myCan = new TCanvas("myCan","Results",800,600);
    myCan->Divide(2,2);
    myCan->cd(1);
    nSamples->Draw();
    myCan->cd(2);
    offsetPosition->Draw();
    myCan->cd(3);
    nSamplesFirst->Draw();
    myCan->cd(4);
    OffsetvsNBins->Draw("colz");

    TCanvas* myCan4 = new TCanvas("myCan4","Offset",800,600);
    offsetPosition->Draw();

    TCanvas* myCan3 = new TCanvas("myCan3","OffsetvsNBins",800,600);
    OffsetvsNBins->Draw("colz");

    return 0;
}

int AnalyzeEvent()
{

	TString fileName = "/Data/BaikalData/2017/cluster-1/0411/g0411.raw.events.195.sorted.root";
	// TString fileName = "/Data/BaikalData/2017/cluster-1/g0411.root";

	TFile* file = new TFile(fileName,"READ");

	TTree* tree = (TTree*)file->Get("Events");
	BRawMasterData* rawMasterData = NULL;
    tree->SetBranchAddress("BRawMasterData.", &rawMasterData);
    BRawMasterHeader* rawMasterHeader = NULL;
    tree->SetBranchAddress("BRawMasterHeader.", &rawMasterHeader);

	BExtractedImpulseTel* extImpTel = new BExtractedImpulseTel();
	BGeomTelMy* geom = new BGeomTelMy();
	BExtractorMy* extractor = new BExtractorMy();
	BRawMasterHeaderMy* masterHeader = new BRawMasterHeaderMy();
	BExtractedCrossTalkTel* crossTalkTel = new BExtractedCrossTalkTel();
	extractor->SetImpulseTel(extImpTel);
	extractor->SetGeomTel(geom);
	extractor->SetRawMasterHeader(masterHeader);
	extractor->SetCrossTalk(crossTalkTel);

	TGraph *graph = new TGraph();
	graph->SetName("graph");
	graph->SetLineColor(4);				
	graph->SetMarkerColor(4);				
	graph->SetMarkerStyle(21);				
	graph->SetMarkerSize(0.5);	

	TGraph* timesExtr = new TGraph();
	timesExtr->SetName("timesExtr");
	timesExtr->SetLineColor(2);				
	timesExtr->SetMarkerColor(2);				
	timesExtr->SetMarkerStyle(21);				
	timesExtr->SetMarkerSize(0.5);	

	TH1F* TimeDist = new TH1F("TimeDist","Time distribution of pulses with respect to test pulse; Pulse extracted Time - Test Pulse Time [bins]",1000,-500,500);
	TH1F* TimeDistUpper = new TH1F("TimeDistUpper","Time distribution of pulses with respect to test pulse",1000,-500,500);
	TH2F* timeVsAmpl = new TH2F("timeVsAmpl","Time of detection vs. Amplitude at half width",200,50,100,1000,0,1000);
	TH2F* timeVsCharge = new TH2F("timeVsCharge","Time of detection vs. Charge",100,50,100,1000,0,100000);
	TH2F* timeVsAmplUpper = new TH2F("timeVsAmplUpper","Time of detection vs. Amplitude at half width on UPPER",200,100,200,1000,0,1000);
	TH2F* timeVsChargeUpper = new TH2F("timeVsChargeUpper","Time of detection vs. Charge on UPPER",200,100,200,1000,0,10000);

	cout << "Number of events: " << tree->GetEntries() << endl;

	int channel = 4;

	tree->GetEntry(0);
	int startTime = 1510301938;
	// cout << "Start time: " << startTime << endl;

	int nTestPulses = 0;

	TCanvas* myCan2 = new TCanvas("myCan2","Results",1200,800);

	for (int i = 0; i < tree->GetEntries(); ++i)
	{
		tree->GetEntry(i);
		// if (rawMasterHeader->GetTimePC() - startTime < 800 || rawMasterHeader->GetTimePC() - startTime > 902)
		// 	continue;

		int nSamplesOnFirst = 0;
    	bool testPulse = false;
    	for (int j = 0; j < rawMasterData->GetNumSamples(); ++j)
    	{
    		if (rawMasterData->GetFADCSample(j)->GetNch() == channel && rawMasterData->GetFADCSample(j)->GetNbins() != 0)
    		{
    			if (rawMasterData->GetFADCSample(j)->GetOffset() > 706 && rawMasterData->GetFADCSample(j)->GetOffset() < 712)
    			{
    				testPulse = true;
    			}
    			nSamplesOnFirst++;
    		}
    	}
    	if (!testPulse)
    		continue;
    	// if (nSamplesOnFirst < 2)
    	// 	continue;
    	//nTestPulses++;
    	int nSamples = 0;
    	double testPulseTime = -1;

    	for (int j = 0; j < rawMasterData->GetNumSamples(); ++j)
    	{
			BRawFADCSample* sample = rawMasterData->GetFADCSample(j);
    		if (sample->GetNch() == channel && sample->GetNbins() != 0)
    		{
    			nSamples++;
    			extractor->Analyze(sample);
    		}
    		if (sample->GetNch() == channel+1 && sample->GetNbins() != 0)
    		{
    			nSamples++;
    			extractor->Analyze(sample);
    		}
    		/*if (sample->GetNch() == 1 && sample->GetNbins() != 0)
    		{
    			nSamples++;
    			extractor->Analyze(sample);
    			if (nSamples == extImpTel->GetNimpulse() && testPulseTime != -1 && abs(extImpTel->GetT(nSamples-1) - testPulseTime - 127) < 10)
    			{
    				//TimeDist->Fill(extImpTel->GetT(nSamples-1)-testPulseTime);
    				Int_t nbins = sample->GetNbins();
    				Int_t nfirst = graph->GetN();
    				Int_t offset = sample->GetOffset();
    				Short_t *data = sample->GetData();
    				Int_t npoints = timesExtr->GetN();
    				// graph->Set(nfirst + nbins);
    				timesExtr->Set(npoints + 1);
    				timesExtr->SetPoint(npoints,extImpTel->GetT(nSamples-1)-testPulseTime,extImpTel->GetA(nSamples-1)/2);
    				timeVsAmplUpper->Fill(extImpTel->GetT(nSamples-1)-testPulseTime,extImpTel->GetA(nSamples-1)/2);
    				timeVsChargeUpper->Fill(extImpTel->GetT(nSamples-1)-testPulseTime,extImpTel->GetQ(nSamples-1));
    				for(Int_t n = 0; n < nbins; n++) {
    					graph->SetPoint(nfirst + n, n+offset-testPulseTime, data[n]-extImpTel->At(nSamples-1)->GetPed());
					}
    			}
    		}*/
    	}

    	for (int j = 0; j < extImpTel->GetNimpulse(); ++j)
    	{
    		if (extImpTel->GetNch(j)%12 == channel && extImpTel->GetA(j) > 200 && extImpTel->GetA(j) < 500 && extImpTel->At(j)->GetTFWHM() > 17.8 && extImpTel->At(j)->GetTFWHM() < 18.8 && extImpTel->At(j)->GetTfront() > 8.5 && extImpTel->At(j)->GetTfront() < 11)
    		{
    			testPulseTime = extImpTel->GetT(j);
    			nTestPulses++;
    		}
    	}
    	for (int j = 0; j < extImpTel->GetNimpulse(); ++j)
    	{
    		if (extImpTel->GetNch(j)%12 == channel &&  testPulseTime != -1 && extImpTel->GetQ(j) > 0 && abs(extImpTel->GetT(j) - testPulseTime - 115) < 5)
    		{
				timeVsAmpl->Fill((extImpTel->GetT(j)-testPulseTime-100)*5,extImpTel->GetA(j)/2);
				timeVsCharge->Fill((extImpTel->GetT(j)-testPulseTime-100)*5,extImpTel->GetQ(j));    			
    		}
    	}
    	for (int j = 0; j < extImpTel->GetNimpulse(); ++j)
    	{
    		if (extImpTel->GetNch(j)%12 == channel + 1 && testPulseTime != -1 && extImpTel->GetQ(j) > 0 && abs(extImpTel->GetT(j) - testPulseTime - 128) < 4)
    		{
				timeVsAmplUpper->Fill((extImpTel->GetT(j)-testPulseTime-100)*5,extImpTel->GetA(j)/2);
				timeVsChargeUpper->Fill((extImpTel->GetT(j)-testPulseTime-100)*5,extImpTel->GetQ(j));    			
    		}
    	}


  //   	cout << "Nsamples: " << nSamples << " " << extImpTel->GetNimpulse() << endl;
  //   	for (int j = 0; j < extImpTel->GetNimpulse(); ++j)
  //   	{
  //   		cout << "Imp " << j << ": " << extImpTel->GetT(j) << " " << extImpTel->GetQ(j) << " " << extImpTel->GetA(j) << endl;
  //   	}
  //   	graph->Draw();
		// myCan2->Update();
		// cin >> nSamples;
		// graph->Clear();
		// graph->Set(0);

    	// cout << extImpTel->GetNimpulse() << endl;

    	extImpTel->Clear();
	}

	cout << "Number of test pulses: " << nTestPulses << endl;

	//extImpTel->At(j)->GetT()
	// TimeDist->Draw();

	TCanvas* c_waveforms = new TCanvas("c_waveforms","Waveforms",800,600);
	graph->Draw("APL");
	timesExtr->Draw("SAMEP");
	TCanvas* c_timeVsAmpl = new TCanvas("c_timeVsAmpl","Time vs. Amplitudes",800,600);
	timeVsAmpl->Draw("colz");
	TCanvas* c_timeVsAmplProfile = new TCanvas("c_timeVsAmplProfile","Time vs. Amplitudes - profile",800,600);
	timeVsAmpl->ProfileY()->Draw();
	TCanvas* c_timeVsCharge = new TCanvas("c_timeVsCharge","Time vs. Charge",800,600);
	timeVsCharge->Draw("colz");
	TCanvas* c_timeVsChargeProfile = new TCanvas("c_timeVsChargeProfile","Time vs. Charge - profile",800,600);
	timeVsCharge->ProfileY()->Draw();
	TCanvas* c_timeDist = new TCanvas("c_timeDist","TimeDist",800,600);
	TimeDist->Draw();

	TCanvas* c_timeVsAmplU = new TCanvas("c_timeVsAmplU","Time vs. Amplitudes - upper",800,600);
	timeVsAmplUpper->Draw("colz");
	TCanvas* c_timeVsAmplProfileU = new TCanvas("c_timeVsAmplProfileU","Time vs. Amplitudes - profile  - upper",800,600);
	timeVsAmplUpper->ProfileY()->Draw();
	TCanvas* c_timeVsChargeU = new TCanvas("c_timeVsChargeU","Time vs. Charge  - upper",800,600);
	timeVsChargeUpper->Draw("colz");
	TCanvas* c_timeVsChargeProfileU = new TCanvas("c_timeVsChargeProfileU","Time vs. Charge - profile  - upper",800,600);
	timeVsChargeUpper->ProfileY()->Draw();

	TCanvas* c_timeSpecSample = new TCanvas("c_timeSpecSample","Time hist samples",800,600);
	c_timeSpecSample->Divide(2,2);
	c_timeSpecSample->cd(1);
	timeVsCharge->ProjectionX("_a",1,50)->Draw();
	c_timeSpecSample->cd(2);
	timeVsCharge->ProjectionX("_b",2,2)->Draw();
	c_timeSpecSample->cd(3);
	timeVsCharge->ProjectionX("_c",3,3)->Draw();
	c_timeSpecSample->cd(4);
	timeVsCharge->ProjectionX("_d",4,4)->Draw();


	return 0;
}

const int nValues = 40;
const int Minimum = -15;
const int nSamples = 1000;


// Study the effect of final sampling with help of Gumpbel function
int SamplingStudy(void)
{
	TF1* singlePE = new TF1("singlePE","53.86*exp(-((x)/2.1978 + exp(-(x)/2.1978)))",-50,50);
	singlePE->SetNpx(1000);
	TCanvas* c_singlePE = new TCanvas("c_singlePE","Single p.e. pulse",1600,1200);
	singlePE->Draw();

	TRandom* randomGenerator = new TRandom(time(0));

	TGraph* g_values = new TGraph(nValues);
	g_values->SetLineColor(4);				
	g_values->SetMarkerColor(4);				
	g_values->SetMarkerStyle(8);				
	g_values->SetMarkerSize(1.5);	

	BRawFADCSample* fadcSample = new BRawFADCSample();
	fadcSample->SetNbins(nValues);

	TH1F* timeDist = new TH1F("timeDist","Time distribution of the extracted pulse times for different sampling",200,-3,-1);

	double realAmplitude = singlePE->GetMaximum();
	// std:;cout << "Real amplitude: " << realAmplitude << std::endl;
	Short_t values[nValues] {0};

	BExtractedImpulseTel* extImpTel = new BExtractedImpulseTel();
	BGeomTelMy* geom = new BGeomTelMy();
	BExtractorMy* extractor = new BExtractorMy();
	BRawMasterHeaderMy* masterHeader = new BRawMasterHeaderMy();
	BExtractedCrossTalkTel* crossTalkTel = new BExtractedCrossTalkTel();
	extractor->SetImpulseTel(extImpTel);
	extractor->SetGeomTel(geom);
	extractor->SetRawMasterHeader(masterHeader);
	extractor->SetCrossTalk(crossTalkTel);

	double minTime = 500;
	double maxTime = -500;

	double realTime = (singlePE->GetX(singlePE->GetMaximum()/2));

	for (int j = 0; j < nSamples+1; ++j)
	{
		// std::cout << "J: " << j << endl;
		for (int i = 0; i < nValues; ++i)
		{
			// std::cout << static_cast<Short_t>(singlePE->Eval(i+Minimum+j*0.1)) << std::endl;
			// std::cout << randomGenerator->Gaus(0,0.52) << std::endl;
			values[i] = static_cast<Short_t>(round(singlePE->Eval(i+Minimum+j*(1.0/nSamples))+randomGenerator->Gaus(0,0.52) + 50));
			// values[i] = static_cast<Short_t>(round (singlePE->Eval(i+Minimum+j*(1.0/nSamples)) + 50));
			// g_values->SetPoint(i,i+Minimum+j*(1.0/nSamples),values[i]-50);
		}
		fadcSample->SetData(values);
		extractor->Analyze(fadcSample);
		// std::cout << extImpTel->GetNimpulse() << std::endl;
		// std::cout << extImpTel->At(j)->GetT()+Minimum+j*(1.0/nSamples) << std::endl;
		if (extImpTel->At(j)->GetT()+Minimum+j*(1.0/nSamples) < minTime)
			minTime = extImpTel->At(j)->GetT()+Minimum+j*(1.0/nSamples);
		if (extImpTel->At(j)->GetT()+Minimum+j*(1.0/nSamples) > maxTime)
			maxTime	= extImpTel->At(j)->GetT()+Minimum+j*(1.0/nSamples);
		timeDist->Fill(extImpTel->At(j)->GetT()+Minimum+j*(1.0/nSamples));
		
		double ampMax = 0;

		for (int i = 0; i < nValues; ++i)
		{
			if (values[i] > ampMax)
			{
				ampMax = values[i];
			}
		}

		ampMax /= 2;
		// std::cout << ampMax << std::endl;
		double time = 0;

		for (int i = 0; i < nValues; ++i)
		{
			if (values[i] > ampMax)
			{
				time = (i-1) + (ampMax-values[i-1])/(values[i]-values[i-1]) + Minimum + j*0.1;
				break;
			}
		}
		// g_values->Draw("samePL");
		// c_singlePE->Update();
		// cin >> ampMax;
		// g_values->Clear();

		// std::cout << "Time: " << time << std::endl;
	}	

	timeDist->Draw();

	std::cout << "Real time [bins]: " << realTime << std::endl;
	std::cout << "Minimal time [bins]: " << minTime << std::endl;
	std::cout << "Maximal time [bins]: " << maxTime << std::endl;
	std::cout << "Maximal diference [ns]: " << (maxTime - minTime)*5 << std::endl;
	//std::cout << extImpTel->GetNimpulse() << std::endl;

	return 0;
}

const int nEventsTWE = 5000;
const int nValuesTWE = 60;
const int MinimumTWE = 15;
const double TTS = 4/5.0;

int TWEStudy(void)
{
	TRandom* randomGenerator = new TRandom();

	TF1* singlePE;

	BRawFADCSample* fadcSample = new BRawFADCSample();
	fadcSample->SetNbins(nValuesTWE);

	Short_t values[nValuesTWE] {0};

	BExtractedImpulseTel* extImpTel = new BExtractedImpulseTel();
	BGeomTelMy* geom = new BGeomTelMy();
	BExtractorMy* extractor = new BExtractorMy();
	BRawMasterHeaderMy* masterHeader = new BRawMasterHeaderMy();
	BExtractedCrossTalkTel* crossTalkTel = new BExtractedCrossTalkTel();
	extractor->SetImpulseTel(extImpTel);
	extractor->SetGeomTel(geom);
	extractor->SetRawMasterHeader(masterHeader);
	extractor->SetCrossTalk(crossTalkTel);

	TH1F* timeDist = new TH1F("timeDist","Time distribution one single p.e. pulse",200,-10,10);

	for (int i = 0; i < nEventsTWE; ++i)
	{
		double timeShift = randomGenerator->Gaus(0,TTS);
		std::string shift = "";
		if (timeShift >= 0)
			shift = "(x-" + std::to_string(timeShift) + ")/2.1978";
		else
			shift = "(x+" + std::to_string(-1*timeShift) + ")/2.1978";
		std::string function = "53.86*exp(-(" + shift + " + exp(-" + shift + ")))";
		// std::cout << function << std::endl;
		singlePE = new TF1("singlePE",function.c_str(),-50,50);
		for (int j = 0 ; j < nValuesTWE ; j++)
		{
			values[j] = static_cast<Short_t>(singlePE->Eval(j+MinimumTWE)+randomGenerator->Gaus(0,0.52) + 50);
			// values[j] = static_cast<Short_t>(singlePE->Eval(j+MinimumTWE) + 50);
			// std::cout << j << " " << values[j] << std::endl;
		}

		fadcSample->SetData(values);
		extractor->Analyze(fadcSample);
		timeDist->Fill(extImpTel->At(i)->GetT()+MinimumTWE);
		// std::cout << extImpTel->GetNimpulse() << std::endl;
		// std::cout << extImpTel->At(i)->GetT()+MinimumTWE << std::endl;

		delete singlePE;
	}

	timeDist->Draw();

	return 0;
}

int TWEStudy2(void)
{
	TRandom* randomGenerator = new TRandom();

	TF1* singlePE;

	BRawFADCSample* fadcSample = new BRawFADCSample();
	fadcSample->SetNbins(nValuesTWE);

	Short_t values[nValuesTWE] {0};

	BExtractedImpulseTel* extImpTel = new BExtractedImpulseTel();
	BGeomTelMy* geom = new BGeomTelMy();
	BExtractorMy* extractor = new BExtractorMy();
	BRawMasterHeaderMy* masterHeader = new BRawMasterHeaderMy();
	BExtractedCrossTalkTel* crossTalkTel = new BExtractedCrossTalkTel();
	extractor->SetImpulseTel(extImpTel);
	extractor->SetGeomTel(geom);
	extractor->SetRawMasterHeader(masterHeader);
	extractor->SetCrossTalk(crossTalkTel);

	TH1F* timeDist = new TH1F("timeDist","Time distribution two single p.e. pulses",200,-10,10);

	for (int i = 0; i < nEventsTWE; ++i)
	{
		double timeShift1 = randomGenerator->Gaus(0,TTS);
		double timeShift2 = randomGenerator->Gaus(0,TTS);
		std::string shift1 = "";
		std::string shift2 = "";
		if (timeShift1 >= 0)
			shift1 = "(x-" + std::to_string(timeShift1) + ")/2.1978";
		else
			shift1 = "(x+" + std::to_string(-1*timeShift1) + ")/2.1978";
		if (timeShift2 >= 0)
			shift2 = "(x-" + std::to_string(timeShift2) + ")/2.1978";
		else
			shift2 = "(x+" + std::to_string(-1*timeShift2) + ")/2.1978";

		std::string function = "53.86*exp(-(" + shift1 + " + exp(-" + shift1 + ")))" + " + 53.86*exp(-(" + shift2 + " + exp(-" + shift2 + ")))";
		// std::cout << function << std::endl;
		singlePE = new TF1("singlePE",function.c_str(),-50,50);
		for (int j = 0 ; j < nValuesTWE ; j++)
		{
			values[j] = static_cast<Short_t>(singlePE->Eval(j+MinimumTWE)+randomGenerator->Gaus(0,0.52) + 50);
			// values[j] = static_cast<Short_t>(singlePE->Eval(j+MinimumTWE) + 50);
			// std::cout << j << " " << values[j] << std::endl;
		}

		fadcSample->SetData(values);
		extractor->Analyze(fadcSample);
		timeDist->Fill(extImpTel->At(i)->GetT()+MinimumTWE);
		// std::cout << extImpTel->GetNimpulse() << std::endl;
		// std::cout << extImpTel->At(i)->GetT()+MinimumTWE << std::endl;

		delete singlePE;
	}

	timeDist->Draw();

	return 0;
}

int GumpVis(void)
{
	//TString fileName = "/Data/BaikalData/2017/cluster-1/g0020_1.root";
	TString fileName = "/Data/BaikalData/2018/cluster-2/0005/h0005.raw.events.196.sorted.root";

	TFile* file = new TFile(fileName,"READ");

	TTree* tree = (TTree*)file->Get("Events");
	BRawMasterData* rawMasterData = NULL;
    tree->SetBranchAddress("BRawMasterData.", &rawMasterData);
    BRawMasterHeader* rawMasterHeader = NULL;
    tree->SetBranchAddress("BRawMasterHeader.", &rawMasterHeader);

    TF1* fitFunc = new TF1("fitFunc","[0]+[1]*exp(-((x-[2])/[3] + exp(-(x-[2])/[3])))",0,200);
    fitFunc->SetLineColor(kRed);

    TF1* fitFunc2 = new TF1("fitFunc2","gaus",0,200);
    fitFunc2->SetLineColor(7);

    TF1* fitFunc3 = new TF1("fitFunc3","landau",0,200);
    fitFunc3->SetLineColor(3);


    TCanvas* myCan = new TCanvas("myCan","Results",800,600);

    TGraph* graph = new TGraph();
    graph->SetName("graph");
	graph->SetLineColor(4);				
	graph->SetMarkerColor(4);				
	graph->SetMarkerStyle(21);				
	graph->SetMarkerSize(0.9);

    cout << tree->GetEntries() << endl;

	tree->GetEntry(35);

	cout << rawMasterData->GetNumSamples() << endl;

	BRawFADCSample* sample = rawMasterData->GetFADCSample(0);

	cout << sample->GetNbins() << endl;

	Int_t nbins = sample->GetNbins();
	Short_t *data = sample->GetData();
	double pedestal = 0;
	int amplitude = 0;
	int charge = 0;
	int amplitudeBin = 0;
	for (int k = 0; k < c_nPed; ++k)
	{
		pedestal += data[k];
	}
	pedestal /= c_nPed;
	for (int k = 0; k < nbins; ++k)
	{
		if (data[k]-pedestal > amplitude)
		{
			amplitude = data[k]-pedestal;
			amplitudeBin = k;
		}
		charge += data[k]-pedestal;
	}
	// if (amplitude-pedestal > 70	|| amplitude-pedestal < 65)
	if (charge < 0)
		cout << "NEGATIVE CHARGE!" << endl;

	graph->Set(nbins);
	for(Int_t n = 0; n < nbins; n++) {
		graph->SetPoint(n, n, data[n]-pedestal);
	}	
	fitFunc->SetParameters(0,amplitude/0.37,amplitudeBin,2.2);
	graph->Fit(fitFunc,"W");
	graph->Fit(fitFunc2,"W+");
	graph->Fit(fitFunc3,"W+");
	// cout << "amplitude = " << amplitude << " charge = " << charge << endl;
	graph->Draw();
	graph->GetXaxis()->SetTitle("Time [FADC channels]");	
	graph->GetYaxis()->SetTitle("Amplitude [FADC channels]");	
	graph->SetTitle("Comparison of different fitting functions");
	myCan->Update();

	TLegend* myLegend = new TLegend(0.6,0.7,0.9,0.9);
	myLegend->AddEntry(fitFunc,"Gumpbel function, #chi^{2}/NDF = 0.32","l");
	myLegend->AddEntry(fitFunc2,"Gaussian function, #chi^{2}/NDF = 2.86","l");
	myLegend->AddEntry(fitFunc3,"Landau function, #chi^{2}/NDF = 2.19","l");
	myLegend->Draw();

    return 0;
}

int GumpVis2(void)
{
	TF1* singlePE = new TF1("singlePE","53.86*exp(-((x)/2.1978 + exp(-(x)/2.1978)))",-15,25);
	singlePE->SetNpx(1000);
	singlePE->SetMinimum(-2);
	singlePE->GetHistogram()->GetYaxis()->SetTitle("Amplitude [FADC channels]");
   	singlePE->GetHistogram()->GetXaxis()->SetTitle("Time [FADC samples]");
	TCanvas* c_singlePE = new TCanvas("c_singlePE","Single p.e. pulse",800,600);
	singlePE->Draw();

	TRandom* randomGenerator = new TRandom(time(0));

	TGraph* g_values = new TGraph(nValues);
	g_values->SetLineColor(4);				
	g_values->SetMarkerColor(4);				
	g_values->SetMarkerStyle(8);				
	g_values->SetMarkerSize(1.5);	

	double realAmplitude = singlePE->GetMaximum();
	// std:;cout << "Real amplitude: " << realAmplitude << std::endl;
	Short_t values[nValues] {0};

	for (int i = 0; i < nValues; ++i)
	{
		// std::cout << static_cast<Short_t>(singlePE->Eval(i+Minimum+j*0.1)) << std::endl;
		// std::cout << randomGenerator->Gaus(0,0.52) << std::endl;
		values[i] = static_cast<Short_t>(round(singlePE->Eval(i-15)));//+randomGenerator->Gaus(0,0.52)));
		// values[i] = static_cast<Short_t>(round (singlePE->Eval(i+Minimum+j*(1.0/nSamples)) + 50));
		g_values->SetPoint(i,i-15,values[i]);
	}
	g_values->Draw("samePL");
	g_values->SetMinimum(-2);
	c_singlePE->Update();
	// cin >> ampMax;
	// g_values->Clear();


	return 0;
}

int GumpVis3(void)
{
	TString fileName = "/Data/BaikalData/muon/0061/0061/f0061.raw.events.214.sorted.root";

	TFile* file = new TFile(fileName,"READ");

	TTree* tree = (TTree*)file->Get("Events");
	BRawMasterData* rawMasterData = NULL;
    tree->SetBranchAddress("BRawMasterData.", &rawMasterData);
    BRawMasterHeader* rawMasterHeader = NULL;
    tree->SetBranchAddress("BRawMasterHeader.", &rawMasterHeader);

    TF1* fitFunc = new TF1("fitFunc","[0]+[1]*exp(-((x-[2])/[3] + exp(-(x-[2])/[3]))) + [4]*exp(-((x-[5])/[6] + exp(-(x-[5])/[6])))",0,200);
    TF1* fitFunc2 = new TF1("fitFunc","[0]+[1]*exp(-((x-[2])/[3] + exp(-(x-[2])/[3])))",0,200);
    TF1* fitFunc3 = new TF1("fitFunc","[0]+[1]*exp(-((x-[2])/[3] + exp(-(x-[2])/[3])))",0,200);


    TCanvas* myCan = new TCanvas("myCan","Results",800,600);

    TGraph* graph = new TGraph();
    graph->SetName("graph");
	graph->SetLineColor(4);				
	graph->SetMarkerColor(4);				
	graph->SetMarkerStyle(21);				
	graph->SetMarkerSize(0.9);

    cout << tree->GetEntries() << endl;

	tree->GetEntry(370);

	cout << rawMasterData->GetNumSamples() << endl;

	BRawFADCSample* sample = rawMasterData->GetFADCSample(2);

	cout << sample->GetNbins() << endl;

	Int_t nbins = sample->GetNbins();
	Short_t *data = sample->GetData();
	double pedestal = 0;
	int amplitude = 0;
	int charge = 0;
	int amplitudeBin = 0;
	for (int k = 0; k < c_nPed; ++k)
	{
		pedestal += data[k];
	}
	pedestal /= c_nPed;
	for (int k = 0; k < nbins; ++k)
	{
		if (data[k]-pedestal > amplitude)
		{
			amplitude = data[k]-pedestal;
			amplitudeBin = k;
		}
		charge += data[k]-pedestal;
	}
	// if (amplitude-pedestal > 70	|| amplitude-pedestal < 65)
	if (charge < 0)
		cout << "NEGATIVE CHARGE!" << endl;

	graph->Set(nbins);
	for(Int_t n = 0; n < nbins; n++) {
		graph->SetPoint(n, n, data[n]-pedestal);
	}	
	// fitFunc->SetParameters(0,amplitude/0.37,amplitudeBin,2.2);
	// fitFunc->SetParameters(0,1315,13,2.2,300,25,2.2);
	fitFunc->SetParameters(0,702,13,2.2,135,22,2.2);
	graph->Fit(fitFunc,"W");
	fitFunc2->SetParameters(fitFunc->GetParameter(0),fitFunc->GetParameter(1),fitFunc->GetParameter(2),fitFunc->GetParameter(3));
	fitFunc2->SetLineColor(4);
	fitFunc2->SetLineStyle(4);
	fitFunc3->SetParameters(fitFunc->GetParameter(0),fitFunc->GetParameter(4),fitFunc->GetParameter(5),fitFunc->GetParameter(6));
	fitFunc3->SetLineColor(4);
	fitFunc3->SetLineStyle(4);
	// cout << "amplitude = " << amplitude << " charge = " << charge << endl;
	graph->DrawClone();
	graph->GetXaxis()->SetTitle("Time [FADC channels]");	
	graph->GetYaxis()->SetTitle("Amplitude [FADC channels]");	
	graph->SetTitle("Pulse separation - step 1");
	graph->DrawClone();

	fitFunc2->Draw("Same");
	fitFunc3->Draw("same");

	myCan->Update();

	for(Int_t n = 0; n < nbins; n++) {
		graph->SetPoint(n, n, data[n]-pedestal-fitFunc2->Eval(n));
	}	

	graph->SetTitle("Pulse separation - step 2");
    TCanvas* myCan2 = new TCanvas("myCan2","Results",800,600);
    graph->Draw();

    return 0;
}

TCanvas* myCan = new TCanvas("myCan","Results");

Double_t GetRealThres(Double_t recentThres, Int_t nPulses, Int_t smoothing, Double_t &width)
{
	Double_t beta = 1.813;
	Double_t thresFADC = (2.605 + 61.36)*recentThres/2.7183;

	const int nSamples = 20;
	const int minSample = -5;
	const int nBins = 2000;

	// myCan->cd();
	TH1I* hist = new TH1I("hist","Pulses detected with 1.5 p.e. threshold, smoothing = 4;Pulse amplitude [p.e.];Simulated and Detected [#]",nBins,0,20);
	TH1I* hist2 = new TH1I("hist2","Hist",nBins,0,20);
	TH1F* ratio = new TH1F("ratio","Ratio",nBins,0,20);
	TRandom2* randomGenerator = new TRandom2(time(0));
	TF1* singlePE;

	Short_t values[nSamples] {0};

	for (int i = 0; i < nPulses; ++i)
	{
		if (i % (nPulses/10) == 0)
		{
			cout << (double)i/nPulses*100 << "% ";
			cout << flush;
		}
		Double_t nPE = randomGenerator->Uniform(recentThres,2*recentThres); 
		Double_t sampleShift = randomGenerator->Uniform(1);
		hist->Fill(nPE);
		Double_t a = 2.605 + 61.36*nPE;
		// string function = std::to_string(a) + "*exp(-((x-" + std::to_string(sampleShift) + ")/" + std::to_string(beta) + " + exp(-(x-" + std::to_string(sampleShift) + ")/" + std::to_string(beta) +")))";
		string function = "[0]*exp(-((x-[1])/[2] + exp(-(x-[1])/[2])))";
		// string function = std::to_string(a) + "*exp(-((x)/" + std::to_string(beta) + " + exp(-(x)/" + std::to_string(beta) +")))";
		singlePE = new TF1("singlePE",function.c_str(),-50,50,TF1::EAddToList::kNo);
		singlePE->SetParameters(a,sampleShift,beta);
		for (int j = 0 ; j < nSamples ; j++)
		{
			values[j] = static_cast<Short_t>(singlePE->Eval(j+minSample)+randomGenerator->Gaus(0,0.52));
			// values[j] = static_cast<Short_t>(singlePE.Eval(j+MinimumTWE) + 50);
			// std::cout << j << " " << values[j] << std::endl;
		}
		delete singlePE;

		Bool_t detected = false;
		for (int i = 0; i < nSamples-smoothing+1; ++i)
		{
			Double_t averageAmp = 0;
			for (int j = 0; j < smoothing; ++j)
			{
				averageAmp += values[i+j];
			}
			averageAmp /= smoothing;
			if (averageAmp > thresFADC)
			{
				detected = true;
				break;
			}
		}
		
		if (detected)
			hist2->Fill(nPE);
		hist2->SetLineColor(kRed);

	}
	cout << endl;

	Bool_t isFound10 = false;
	Bool_t isFound50 = false;
	Bool_t isFound90 = false;
	Double_t valueAt50P = -1;
	Double_t valueAt10P = -1;
	Double_t valueAt90P = -1;

	for (int i = 0; i < nBins; ++i)
	{
		if (hist->GetBinContent(i) != 0)
		{
			ratio->SetBinContent(i,(double)hist2->GetBinContent(i)/hist->GetBinContent(i));
			if (!isFound10 && (double)hist2->GetBinContent(i)/hist->GetBinContent(i) > 0.1)
			{
				isFound10 = true;
				valueAt10P = hist->GetBinCenter(i);
			}			
			if (!isFound50 && (double)hist2->GetBinContent(i)/hist->GetBinContent(i) > 0.5)
			{
				isFound50 = true;
				valueAt50P = hist->GetBinCenter(i);
			}			
			if (!isFound90 && (double)hist2->GetBinContent(i)/hist->GetBinContent(i) > 0.9)
			{
				isFound90 = true;
				valueAt90P = hist->GetBinCenter(i);
			}			
		}
		// cout << hist2->GetBinContent(i) << " " << hist->GetBinContent(i) << endl;
	}

	// cout << valueAt10P << " " << valueAt50P << " " << valueAt90P << endl;
	width = valueAt90P-valueAt10P;

	// How to obtain x values at 10%, 50% and 90% with Erf(x), but it works just for small thresholds (1.5 p.e.)
	// TF1* erf = new TF1("erf","0.5+0.5*TMath::Erf([0]*x-[1])",0,20);
	// erf->SetParameters(0.1,recentThres);
	// ratio->Fit("erf");
	// cout << erf->GetX(0.1) << " " << erf->GetX(0.5) << " " << erf->GetX(0.9) << endl;

	hist->DrawClone();
	// hist2->DrawClone("same");

	myCan->Modified();
    myCan->Update();

	ratio->DrawClone("same");

	delete hist;
	delete hist2;
	delete ratio;
	delete randomGenerator;

	return valueAt50P;
}

// to model trigger thresholds. Main idea from Bair. Try to find a coefficient between set threshold (1.5 and 4) 
// and observed one caused by smoothing triggering
// Different smoothing windows 1,2,4,8
// Study dependence of the constant on the amplitude 
int modelTriggerThreshold(void)
{
	clock_t begin = clock();
	Int_t nPulses = 1000000;
	Double_t lowestThres = 0.5;
	Int_t nThres = 20;
	Int_t smoothing = 4;

	// myCan->cd();
    // myCan->Modified();
    // myCan->Update();

	TGraph* g_realThres = new TGraph(nThres);
	TGraph* g_width = new TGraph(nThres);
	TGraph* g_widthVsRealThres = new TGraph(nThres);
 
	for (int i = 0; i < nThres; ++i)
	{
		Double_t realWidth = -1;
		Double_t recentThres = lowestThres + i*0.5;
		Double_t realThres = GetRealThres(recentThres,nPulses,smoothing,realWidth);
		cout << recentThres << " " << realThres << " " << realWidth << endl;
		g_realThres->SetPoint(i,recentThres,realThres);
		g_width->SetPoint(i,recentThres,realWidth);
		g_widthVsRealThres->SetPoint(i,recentThres,realWidth/realThres);
	}

	TCanvas* resultsThres = new TCanvas("resultsThres","Results thresholds",800,600);
	g_realThres->SetMarkerSize(2);
	g_realThres->SetMarkerColor(kBlue);
	g_realThres->SetMarkerStyle(2);
	g_realThres->GetXaxis()->SetTitle("Set Threshold [p.e.]");
	g_realThres->GetYaxis()->SetTitle("Real Threshold [p.e.]");
	string resultsThresTitle = "Set vs. Real Threshold, smoothing = " + std::to_string(smoothing);
	g_realThres->SetTitle(resultsThresTitle.c_str());
	g_realThres->Draw("AP");

	TCanvas* resultsWidth = new TCanvas("resultsWidth","Results width",800,600);
	g_width->SetMarkerSize(2);
	g_width->SetMarkerColor(kBlue);
	g_width->SetMarkerStyle(2);
	g_width->GetXaxis()->SetTitle("Set Threshold [p.e.]");
	g_width->GetYaxis()->SetTitle("Width [p.e.]");
	string resultsWidthTitle = "Set Threshold vs. Width_{10%}^{90%}, smoothing = " + std::to_string(smoothing);
	g_width->SetTitle(resultsWidthTitle.c_str());
	g_width->Draw("AP");

	TCanvas* resultsWidthVsRealThres = new TCanvas("resultsWidthVsRealThres","Results width/Real Thres",800,600);
	g_widthVsRealThres->SetMarkerSize(2);
	g_widthVsRealThres->SetMarkerColor(kBlue);
	g_widthVsRealThres->SetMarkerStyle(2);
	g_widthVsRealThres->GetXaxis()->SetTitle("Set Threshold [p.e.]");
	g_widthVsRealThres->GetYaxis()->SetTitle("Width/Real Threshold [#]");
	string resultsWidthVsRealThresTitle = "Set Threshold vs. Width/Real Threshold, smoothing = " + std::to_string(smoothing);
	g_widthVsRealThres->SetTitle(resultsWidthVsRealThresTitle.c_str());
	g_widthVsRealThres->Draw("AP");

	clock_t end = clock();
    cout << endl << "Elapsed time : " << double(end - begin)/CLOCKS_PER_SEC << endl;

	return 1;
}
