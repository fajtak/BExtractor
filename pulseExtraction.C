#include <vector>

#include "BRawMasterHeader.h"
#include "BRawMasterData.h"

#include "BExtractedImpulseTel.h"

#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TMath.h"
#include "TLine.h"
#include "TString.h"

#include "TMVA/Reader.h"

struct DP
{
  Float_t nSamples;     // Number of samples in the whole waveform
  Float_t nSamplesExt;  // Number of samples in the extracted pulse
  Float_t QAratio;    // pulse charge/Amplitude ratio
  Float_t A_dp;       // Amplitude in the position of DP
  Float_t A_1;      // Amplitude of the left peak
  Float_t A_2;      // Amplitude of the right peak
  Float_t L_1;      // distance between DP and left peak
  Float_t L_2;      // distance between DP and right peak
  Float_t Q_1;      // charge of the left peak
  Float_t Q_2;      // charge of the right peak
  Float_t Prel_DP;    // relative position of the DP in the waveform
  Float_t Prel_DPExt;   // relative position of the DP in the extracted pulse
  Float_t P_1;      // position of the first peak
  Float_t P_2;      // position of the second peak
  Float_t P_DP;     // position of the DP
};

// Double pulse analysis variables
TMVA::Reader fReader;
DP fFoundDP;
std::vector<DP> v_foundDPs;
Double_t fBDTCut;

// TString fileName = "/Data/BaikalData/exp16_barsv051/cluster0/0070/f0070.raw.events.212.sorted.root";
Int_t fSizeAvgWindowB = 8; // number of samples that are used to calculate pedestal value at the beginning of the waveform
Int_t fSizeAvgWindowE = 5; // number of samples that are used to calculate pedestal value at the end of the waveform
Int_t fNSigmas = 4; // number of sigmas of pedestal that has to be reached to start pulse processing

// global histograms
TCanvas* myCan = new TCanvas("myCan","Results",1024,720);
TH1F* h_initPedRMS = new TH1F("h_initPedRMS","Initial pedestal RMS; Ped_{init} [FADC channels]; NoE [#]",1000,0,100);
TH1F* h_begEndDiff = new TH1F("h_begEndDiff","First - Last sample; First - Last [FADC channels]; NoE [#]",1000,0,100);

// root file where the visualizations of waveforms can be saved
TFile* outputFile;

void InitializeBDT()
{
    // Create a set of variables and declare them to the reader
    // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
    // fBDTCut = -0.03;
    fBDTCut = 0.4;
    fReader.AddVariable( "MVA_nbins_Extracted", &fFoundDP.nSamplesExt);
    fReader.AddVariable( "MVA_relative_position_Extracted", &fFoundDP.Prel_DPExt);
    fReader.AddVariable( "MVA_position_of_doublepeak_Extracted", &fFoundDP.P_DP);
    fReader.AddVariable( "MVA_position_of_amp1doublepeak_Extracted", &fFoundDP.P_1);
    fReader.AddVariable( "MVA_position_of_amp2doublepeak_Extracted", &fFoundDP.P_2);
    fReader.AddVariable( "MVA_amplitude1_of_doublepeak", &fFoundDP.A_1);
    fReader.AddVariable( "MVA_amplitude2_of_doublepeak", &fFoundDP.A_2);
    fReader.AddVariable( "MVA_amplitude_of_positionDP", &fFoundDP.A_dp);
    fReader.AddVariable( "MVA_Q1_Extracted", &fFoundDP.Q_1);
    fReader.AddVariable( "MVA_Q2_Extracted", &fFoundDP.Q_2);
    fReader.AddVariable( "MVA_dist_Amp1_posDP", &fFoundDP.L_1);
    fReader.AddVariable( "MVA_dist_Amp2_posDP", &fFoundDP.L_2);
    TString weights_path = std::string(std::getenv("BARSSYS"))+"/config/DP12_BDT_weights2.xml";
    fReader.BookMVA( "BDT", weights_path );
}


// Drawing of global histograms that are defined above
int DrawResults()
{
	TCanvas* c_initPed = new TCanvas("c_initPed","InitialPedestal");
	h_initPedRMS->Draw();

	TCanvas* c_begEndDiff = new TCanvas("c_begEndDiff","BegEndDiff");
	h_begEndDiff->Draw();

	return 0;
}

int DrawWaveform(Short_t* data, Int_t nbins, float average, int offset, int t1, int t1new, int t2, int t2new, float amp, double threshold)
{
    myCan->cd();
    TGraph* waveform = new TGraph(nbins);
    for (int i = 0; i < nbins; ++i)
    {
        waveform->SetPoint(i,i+offset,data[i]-average);
    }
	waveform->SetTitle("");
	waveform->GetXaxis()->SetTitle("Time [FADC channels]");
	waveform->GetYaxis()->SetTitle("Amplitude [FADC channels]");
	waveform->Draw("APL*");
	outputFile->cd();
	waveform->Write();

    TLine* pedestalLine = new TLine(0,0,nbins,0);
    pedestalLine->SetLineColor(kRed);
    pedestalLine->SetLineWidth(3);
    pedestalLine->SetLineStyle(6);
    pedestalLine->Draw();

    TLine* t1Line = new TLine(t1,0,t1,amp);
    t1Line->SetLineColor(kGreen);
    t1Line->SetLineWidth(3);
    t1Line->SetLineStyle(6);
    t1Line->Draw();

    TLine* newT1Line = new TLine(t1new,0,t1new,amp);
    newT1Line->SetLineColor(kOrange);
    newT1Line->SetLineWidth(3);
    newT1Line->SetLineStyle(6);
    newT1Line->Draw();

    TLine* t2Line = new TLine(t2,0,t2,amp);
    t2Line->SetLineColor(kGreen);
    t2Line->SetLineWidth(3);
    t2Line->SetLineStyle(6);
    t2Line->Draw();

    TLine* newT2Line = new TLine(t2new,0,t2new,amp);
    newT2Line->SetLineColor(kOrange);
    newT2Line->SetLineWidth(3);
    newT2Line->SetLineStyle(6);
    newT2Line->Draw();

    TLine* thresLine = new TLine(0,threshold-average,nbins,threshold-average);
    thresLine->SetLineColor(kYellow);
    thresLine->SetLineWidth(3);
    thresLine->SetLineStyle(6);
    thresLine->Draw();

    TLine* halfAmplitude = new TLine(0,amp/2,nbins,amp/2);
    halfAmplitude->SetLineColor(kBlue);
    halfAmplitude->SetLineWidth(3);
    halfAmplitude->SetLineStyle(6);
    halfAmplitude->Draw();

    // Label();

    myCan->Modified();
    myCan->Update();
    int dummy;
    // cin >> dummy;
    delete waveform;
    return dummy;
}

int DrawWaveform(Short_t* data, Int_t nbins, float average, int offset, int t1, int t1new, int t2, int t2new, float amp, double threshold, TF1* fitFunc)
{
    myCan->cd();
    TGraph* waveform = new TGraph(nbins);
    for (int i = 0; i < nbins; ++i)
    {
        waveform->SetPoint(i,i+offset,data[i]-average);
    }
	waveform->SetTitle("");
	waveform->GetXaxis()->SetTitle("Time [FADC channels]");
	waveform->GetYaxis()->SetTitle("Amplitude [FADC channels]");
    waveform->SetLineWidth(3);
    waveform->Draw("APL*");
    outputFile->cd();
    waveform->Write();

    TLine* pedestalLine = new TLine(0,0,nbins,0);
    pedestalLine->SetLineColor(kRed);
    pedestalLine->SetLineWidth(3);
    pedestalLine->SetLineStyle(6);
    pedestalLine->Draw();

    TLine* t1Line = new TLine(t1,0,t1,amp);
    t1Line->SetLineColor(kGreen);
    t1Line->SetLineWidth(3);
    t1Line->SetLineStyle(6);
    // t1Line->Draw();

    TLine* newT1Line = new TLine(t1new,0,t1new,amp);
    newT1Line->SetLineColor(kOrange);
    newT1Line->SetLineWidth(3);
    newT1Line->SetLineStyle(6);
    // newT1Line->Draw();

    TLine* t2Line = new TLine(t2,0,t2,amp);
    t2Line->SetLineColor(kGreen);
    t2Line->SetLineWidth(3);
    t2Line->SetLineStyle(6);
    // t2Line->Draw();

    TLine* newT2Line = new TLine(t2new,0,t2new,amp);
    newT2Line->SetLineColor(kOrange);
    newT2Line->SetLineWidth(3);
    newT2Line->SetLineStyle(6);
    // newT2Line->Draw();

    TLine* thresLine = new TLine(0,threshold-average,nbins,threshold-average);
    thresLine->SetLineColor(kOrange);
    thresLine->SetLineWidth(3);
    thresLine->SetLineStyle(6);
    thresLine->Draw();

    TLine* halfAmplitude = new TLine(0,amp/2,nbins,amp/2);
    halfAmplitude->SetLineColor(kBlue);
    halfAmplitude->SetLineWidth(3);
    halfAmplitude->SetLineStyle(6);
    halfAmplitude->Draw();

    TLine* extrTime = new TLine(fitFunc->GetParameter(2)-fitFunc->GetParameter(3)*0.985199,0,fitFunc->GetParameter(2)-fitFunc->GetParameter(3)*0.985199,amp);
    extrTime->SetLineColor(kPink);
    extrTime->SetLineWidth(3);
    extrTime->SetLineStyle(6);
    extrTime->Draw();

    fitFunc->Draw("same");

    // Label();

    myCan->Modified();
    myCan->Update();
    int dummy;
    // cin >> dummy;
    delete waveform;
    return dummy;
}

// It calculates average pedestal and pedestalRMS from given bins [firstBin,lastBin)
// Can be used for pedestal estimation at the beginning as well as end of the pulse
void InitialPed(BRawFADCSample* ch, Double_t& ped, Double_t& pedrms, Int_t firstBin, Int_t lastBin)
{
  Short_t* data = ch->GetData();

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

  // pedrms = 0;
  // for(int i = 0; i < nbins; i++) {
  //  pedrms += (data[i] - ped) * (data[i] - ped);
  // }
  // pedrms = TMath::Sqrt(pedrms/nbins);
}

void FindNeighboringAmplitudes(Short_t* data, int DPposition, int t1, int t2, double average)
{
    for (int k = DPposition + 1; k <= t2; ++k)
    {
        if (data[k] >= data[k-1])
        {
            fFoundDP.A_2 = data[k]-average;
            fFoundDP.P_2 = k - t1;
            fFoundDP.L_2 = k - DPposition;
        }else
        {
            break;
        }
    }
    for (int k = DPposition; k >= t1 ; --k)
    {
        if (data[k-1] >= data[k])
        {
            fFoundDP.A_1 = data[k-1]-average;
            fFoundDP.P_1 = k - 1 - t1;
            fFoundDP.L_1 = DPposition - k + 1;
        }else
        {
            break;
        }
    }
}

void QProcessing(Short_t* data, int posDP, int t1, int t2, double average)
{
    double Q = 0;
    double A = 0;
    fFoundDP.Q_1 = 0;
    fFoundDP.Q_2 = 0;
    for (int i = t1; i <= t2; ++i)
    {
        Q += data[i] - average;
        if (data[i]-average > A)
            A = data[i] - average;
        if (i <= posDP)
            fFoundDP.Q_1 += data[i] - average;
        if (i >= posDP)
            fFoundDP.Q_2 += data[i] - average;
    }
    fFoundDP.QAratio = Q/A;
}

// search for the local minimum and test if it is DP with BDTs
Int_t FindDP(Short_t* data, Int_t t1, Int_t t2, Float_t pedestal, Int_t nbins)
{
    Int_t nFoundDPs = 0;
    Int_t previousDerivative = 1;
    Int_t nZeros = 0;
    for (int i = t1; i < t2; ++i)
    {
        Int_t newDerivative = data[i+1] - data[i];
        if (newDerivative >= 0 && previousDerivative < 0)
        {
            if (newDerivative == 0)
            {
                nZeros++;
                continue;
            }
            // cout << "DP found: " << i << endl;
            fFoundDP.A_2 = data[i]-pedestal;
            fFoundDP.A_1 = data[i]-pedestal;
            fFoundDP.nSamples = nbins;
            fFoundDP.nSamplesExt = t2- t1;
            fFoundDP.A_dp = data[i]-pedestal;
            fFoundDP.P_DP = (i-nZeros/2) - t1;
            fFoundDP.Prel_DPExt = (double)(i-nZeros/2)/fFoundDP.nSamplesExt*100;
            fFoundDP.Prel_DP = (double)(i-nZeros/2)/nbins*100;
            FindNeighboringAmplitudes(data,i-nZeros/2,t1,t2,pedestal);
            QProcessing(data,i-nZeros/2,t1,t2,pedestal);
            // cout << "Amp1: " << amp1 << " Amp2: " << amp2 << " AmpDP: " << data[i]-pedestal << " Dist1: " << dist1 << " Dist2: " << dist2 << endl;
            Double_t BDTvalue = fReader.EvaluateMVA("BDT");
            if (BDTvalue > fBDTCut)
            {
              v_foundDPs.push_back(fFoundDP);
                nFoundDPs++;
            }
        }
        else
          nZeros = 0;
        previousDerivative = newDerivative;
    }
    return nFoundDPs;
}

Double_t ChargeCorrBefore(Short_t* data, double pedestal, int &newT1)
{
  Short_t previousValue = data[newT1];
  double qCorr = 0;
  for (int i = newT1-1; i >= 0; i--)
  {
      if (data[i]-pedestal > 0 && previousValue > data[i])
      {
          qCorr += data[i]-pedestal;
          previousValue = data[i];
      }else
      {
        newT1 = i+1;
          break;
      }
  }
  return qCorr;
}

Double_t ChargeCorrAfter(Short_t* data, double pedestal, int &newT2, double posThres, int nbins)
{
  double qCorr = 0;
  Double_t lowerPosThres = (posThres-pedestal)/fNSigmas + pedestal;
  for (int j = newT2; j < nbins; ++j)
  {
      if (data[j] > lowerPosThres && data[j] < posThres)
      {
          qCorr += data[j]-pedestal;
      }else
      {
          newT2 = j - 1;
          break;
      }
      if (j == nbins - 1)
      {
          newT2 = j;
      }
  }
  return qCorr;
}

// function to fit a single pulse. If the pulse is saturated, other function has to be used.
TF1* FitPulse(Short_t* data, int t1, int t2, double average, Float_t amplitude,Int_t ampPos, Double_t& beta)
{
    // cout << "amplitude: " << amplitude << " position: " << (t2+t1)/2 << endl;
    TGraph waveform(t2-t1+1);
    for (int i = 0; i < t2-t1+1; ++i)
    {
        waveform.SetPoint(i,t1+i,data[t1+i]-average);
    }

    TF1* fitFunc;

    if (amplitude < 1500)
    {
      fitFunc = new TF1("fitFunc","[0]+[1]*exp(-((x-[2])/[3] + exp(-(x-[2])/[3])))",t1,t2);
      fitFunc->SetParameters(0,amplitude/0.37,ampPos,2.0);
      fitFunc->FixParameter(0,0);
    }else
    {
      fitFunc = new TF1("fitFunc","[0]*(1-[1]*exp([2]*(x-[3])))*(x<[3]+0.5)+[0]*exp([2]*(x-[3]))*(x>=[3]+0.5)",t1+2,t2-10);
      fitFunc->SetParameters(amplitude,0.05,-0.25,ampPos);
      // fitFunc->FixParameter(1,);
      cout << amplitude << " " << "10" << " " << "-0.2" << " " << "20" << endl;
    }


    waveform.Fit(fitFunc,"");

    if (amplitude >= 1500)
      // cout << "time: " << TMath::Log((1-(amplitude/(2*fitFunc->GetParameter(0)))/fitFunc->GetParameter(1)))/fitFunc->GetParameter(2)+fitFunc->GetParameter(3) << " charge: " << fitFunc->Integral(t1,t2) << endl;
      cout << "time: " << TMath::Log((1-(amplitude/(2*fitFunc->GetParameter(0))))/fitFunc->GetParameter(1))/fitFunc->GetParameter(2) + fitFunc->GetParameter(3) << " charge: " << fitFunc->Integral(t1,t2) << endl;
    // cout << fitFunc->GetParameter(0) << " " << fitFunc->GetParameter(1) << " " << fitFunc->GetParameter(2) << " " << fitFunc->GetParameter(3) << endl;
  // cout << "chi2/NDF" << fitFunc->GetChisquare()/fitFunc->GetNDF() << endl;
    // cout << "Integrated charge: " << fitFunc->GetParameter(1)*fitFunc->GetParameter(3) << endl;

    beta = fitFunc->GetParameter(3);

    // globalWaveform = data;
    // globalT1=t1;
    // globalT2= t2;
    // globalAverage = average;

    // TMinuit myMinuit(3);
    // myMinuit.SetFCN(chi2);

    // int ierflg = 0;
    // // myMinuit.mnparm(0,"Amplitude",amplitude/0.37,1.0,0,4000,ierflg);
    // // myMinuit.mnparm(1,"Position",(t2+t1)/2,0.1,0,1024,ierflg);
    // // myMinuit.mnparm(2,"Sigma",2.0,0.01,0,5,ierflg);

    // double arglist[10];
    // myMinuit.mnexcm("MIGRAD",arglist,1,ierflg);

    return fitFunc;
}

// fucntion for the DP fitting. The function is made based on the number of found DPs (global)
// initial parameters are set based on the studies of Gumbel
TF1* FitDP(Short_t* data, double initPed, int t1, int t2)
{
    // cout << "amplitude: " << amplitude << " position: " << (t2+t1)/2 << endl;
    TGraph waveform(t2-t1+1);
    for (int i = 0; i < t2-t1+1; ++i)
    {
        waveform.SetPoint(i,t1+i,data[t1+i]-initPed);
    }

    TString funcEquation = "[0]+[1]*exp(-((x-[2])/[3] + exp(-(x-[2])/[3])))+[4]*exp(-((x-[5])/[6] + exp(-(x-[5])/[6])))";
    if (v_foundDPs.size() > 1)
    {
      for (int i = 0; i < v_foundDPs.size()-1; ++i)
      {
        funcEquation += Form("+[%d]*exp(-((x-[%d])/[%d] + exp(-(x-[%d])/[%d])))",7+i*3,8+i*3,9+i*3,8+i*3,9+i*3);
      }
    }

    TF1* fitFunc = new TF1("fitFunc",funcEquation.Data(),t1,t2);
    fitFunc->SetParameters(0,v_foundDPs[0].A_1/0.37,v_foundDPs[0].P_1+t1,2.0,v_foundDPs[0].A_2/0.37,v_foundDPs[0].P_2+t1,2.0);
    fitFunc->FixParameter(0,0);
    fitFunc->SetParLimits(1,v_foundDPs[0].A_1/0.37*0.75,v_foundDPs[0].A_1/0.37*1.5);
    fitFunc->SetParLimits(2,v_foundDPs[0].P_1-1.5+t1,v_foundDPs[0].P_1+1.5+t1);
    fitFunc->SetParLimits(4,v_foundDPs[0].A_2/0.37*0.2,v_foundDPs[0].A_2/0.37*2.0);
    fitFunc->SetParLimits(5,v_foundDPs[0].P_2-1.5+t1,v_foundDPs[0].P_2+1.5+t1);
    fitFunc->SetParLimits(3,1.5,3);
    fitFunc->SetParLimits(6,1.5,3);

    if (v_foundDPs.size() > 1)
    {
      for (int i = 0; i < v_foundDPs.size()-1; ++i)
      {
        fitFunc->SetParameter(7+i*3,v_foundDPs[1+i].A_2/0.37);
        fitFunc->SetParameter(8+i*3,v_foundDPs[1+i].P_2+t1);
        fitFunc->SetParLimits(8+i*3,v_foundDPs[1+i].P_2+t1-1.5,v_foundDPs[1+i].P_2+t1+1.5);
        fitFunc->SetParameter(9+i*3,2);
        fitFunc->SetParLimits(9+i*3,1.5,3);
      }
    }

    waveform.Fit(fitFunc,"W");
    cout << v_foundDPs[0].A_1/0.37 << " " << v_foundDPs[0].P_1+t1 << " " << 2.0 << " " << v_foundDPs[0].A_2/0.37 << " " << v_foundDPs[0].P_2+t1 << " " << 2.0 << endl;
    // cout << fitFunc->GetParameter(0) << " " << fitFunc->GetParameter(1) << " " << fitFunc->GetParameter(2) << " " << fitFunc->GetParameter(3) << endl;
	// cout << "chi2/NDF" << fitFunc->GetChisquare()/fitFunc->GetNDF() << endl;
    // cout << "Integrated charge: " << fitFunc->GetParameter(1)*fitFunc->GetParameter(3) << endl;

  waveform.SetTitle("");
  waveform.GetXaxis()->SetTitle("Time [FADC channels]");
  waveform.GetYaxis()->SetTitle("Amplitude [FADC channels]");
  waveform.Draw("APL*");

  TF1* functions[v_foundDPs.size()+1];

  for (int i = 0; i < v_foundDPs.size()+1; ++i)
  {
    functions[i] = new TF1(Form("fitfunc%d",i),"[0]+[1]*exp(-((x-[2])/[3] + exp(-(x-[2])/[3])))",t1,t2);
    functions[i]->SetParameters(0,fitFunc->GetParameter(1+3*i),fitFunc->GetParameter(2+3*i),fitFunc->GetParameter(3+3*i));
    functions[i]->SetLineColor(kBlue);
    functions[i]->SetLineStyle(3);
    functions[i]->Draw("same");
  }
  outputFile->cd();

  myCan->Modified();
    myCan->Update();
    int dummy;
    // cin >> dummy;
    // delete waveform;
    // return dummy;
    for (int i = 0; i < v_foundDPs.size()+1; ++i)
    {
      delete functions[i];
      /* code */
    }

    // beta = fitFunc->GetParameter(3);

    // globalWaveform = data;
    // globalT1=t1;
    // globalT2= t2;
    // globalAverage = initPed;

    // TMinuit myMinuit(3);
    // myMinuit.SetFCN(chi2);

    int ierflg = 0;
    // myMinuit.mnparm(0,"Amplitude",amplitude/0.37,1.0,0,4000,ierflg);
    // myMinuit.mnparm(1,"Position",(t2+t1)/2,0.1,0,1024,ierflg);
    // myMinuit.mnparm(2,"Sigma",2.0,0.01,0,5,ierflg);

    double arglist[10];
    // myMinuit.mnexcm("MIGRAD",arglist,1,ierflg);

    return fitFunc;
}

// the same like FitDP but here all the betas are fixed as one value
TF1* FitDPEqual(Short_t* data, double initPed, int t1, int t2)
{
    // cout << "amplitude: " << amplitude << " position: " << (t2+t1)/2 << endl;
    TGraph waveform(t2-t1+1);
    for (int i = 0; i < t2-t1+1; ++i)
    {
        waveform.SetPoint(i,t1+i,data[t1+i]-initPed);
    }

    TString funcEquation = "[1]*exp(-((x-[2])/[0] + exp(-(x-[2])/[0])))+[3]*exp(-((x-[4])/[0] + exp(-(x-[4])/[0])))";
    if (v_foundDPs.size() > 1)
    {
      for (int i = 0; i < v_foundDPs.size()-1; ++i)
      {
        funcEquation += Form("+[%d]*exp(-((x-[%d])/[%d] + exp(-(x-[%d])/[%d])))",5+i*2,6+i*2,0,6+i*2,0);
      }
    }

    TF1* fitFunc = new TF1("fitFunc",funcEquation.Data(),t1,t2);
    fitFunc->SetParameters(2.0,v_foundDPs[0].A_1/0.37,v_foundDPs[0].P_1+t1,v_foundDPs[0].A_2/0.37,v_foundDPs[0].P_2+t1);
    fitFunc->SetParLimits(0,1.5,3);
    fitFunc->SetParLimits(1,v_foundDPs[0].A_1/0.37*0.75,v_foundDPs[0].A_1/0.37*1.5);
    fitFunc->SetParLimits(2,v_foundDPs[0].P_1-1.5+t1,v_foundDPs[0].P_1+1.5+t1);
    fitFunc->SetParLimits(3,v_foundDPs[0].A_2/0.37*0.2,v_foundDPs[0].A_2/0.37*2.0);
    fitFunc->SetParLimits(4,v_foundDPs[0].P_2-1.5+t1,v_foundDPs[0].P_2+1.5+t1);
    // fitFunc->SetParLimits(6,1.5,3);

    if (v_foundDPs.size() > 1)
    {
      for (int i = 0; i < v_foundDPs.size()-1; ++i)
      {
        fitFunc->SetParameter(5+i*2,v_foundDPs[1+i].A_2/0.37);
        fitFunc->SetParameter(6+i*2,v_foundDPs[1+i].P_2+t1);
        fitFunc->SetParLimits(6+i*2,v_foundDPs[1+i].P_2+t1-1.5,v_foundDPs[1+i].P_2+t1+1.5);
        // fitFunc->SetParameter(9+i*3,2);
        // fitFunc->SetParLimits(9+i*3,1.5,3);
      }
    }

    waveform.Fit(fitFunc,"W");
    cout << v_foundDPs[0].A_1/0.37 << " " << v_foundDPs[0].P_1+t1 << " " << 2.0 << " " << v_foundDPs[0].A_2/0.37 << " " << v_foundDPs[0].P_2+t1 << " " << 2.0 << endl;
    // cout << fitFunc->GetParameter(0) << " " << fitFunc->GetParameter(1) << " " << fitFunc->GetParameter(2) << " " << fitFunc->GetParameter(3) << endl;
  // cout << "chi2/NDF" << fitFunc->GetChisquare()/fitFunc->GetNDF() << endl;
    // cout << "Integrated charge: " << fitFunc->GetParameter(1)*fitFunc->GetParameter(3) << endl;

  waveform.SetTitle("");
  waveform.GetXaxis()->SetTitle("Time [FADC channels]");
  waveform.GetYaxis()->SetTitle("Amplitude [FADC channels]");
  waveform.Draw("APL*");

  TF1* functions[v_foundDPs.size()+1];

  for (int i = 0; i < v_foundDPs.size()+1; ++i)
  {
    functions[i] = new TF1(Form("fitfunc%d",i),"[1]*exp(-((x-[2])/[0] + exp(-(x-[2])/[0])))",t1,t2);
    functions[i]->SetParameters(fitFunc->GetParameter(0),fitFunc->GetParameter(1+2*i),fitFunc->GetParameter(2+2*i));
    functions[i]->SetLineColor(kBlue);
    functions[i]->SetLineStyle(3);
    functions[i]->Draw("same");
  }
  outputFile->cd();

  myCan->Modified();
    myCan->Update();
    int dummy;
    // cin >> dummy;
    // delete waveform;
    // return dummy;
    for (int i = 0; i < v_foundDPs.size()+1; ++i)
    {
      delete functions[i];
      /* code */
    }

    // beta = fitFunc->GetParameter(3);

    // globalWaveform = data;
    // globalT1=t1;
    // globalT2= t2;
    // globalAverage = initPed;

    // TMinuit myMinuit(3);
    // myMinuit.SetFCN(chi2);

    int ierflg = 0;
    // myMinuit.mnparm(0,"Amplitude",amplitude/0.37,1.0,0,4000,ierflg);
    // myMinuit.mnparm(1,"Position",(t2+t1)/2,0.1,0,1024,ierflg);
    // myMinuit.mnparm(2,"Sigma",2.0,0.01,0,5,ierflg);

    double arglist[10];
    // myMinuit.mnexcm("MIGRAD",arglist,1,ierflg);

    return fitFunc;
}

Float_t X(BRawFADCSample* ch, BExtractedImpulse* imp, Float_t y)
{
  // get FADC sample
  Short_t* data = ch->GetData();

  Int_t t1 = imp->GetT1();
  Int_t t2 = imp->GetT2();

  // cout << " " << t1 << " " << t2 << " " << y << endl;

  Int_t x = -1;

  for(int i = t1; i <= t2; i++) {
    if(data[i] >= y) {
      x = i;
      break;
    }
  }

  if(x >= 0) {
    if(x != 0 && data[x] != data[x-1])
      return x - 1 + (y - data[x-1]) / (data[x] - data[x-1]);
    else
      return x;
  }
  else {
    cerr << "Time extractor X() Error (x < 0)" << endl;
    return -1;
  }
}

Float_t dT(BRawFADCSample* ch, BExtractedImpulse* imp, Float_t y, double pedestal)
{
  // get FADC sample
  Int_t nbins = ch->GetNbins();
  Short_t* data = ch->GetData();

  Int_t ta = imp->GetTa();
  Float_t a  = imp->GetA();

  Float_t dt = 0;
  if(a + pedestal < y) {
    cout << "Time extractor dT() Error (y > pedestal + amplitude)" << endl;
    exit(1);
  }

  for(int i = ta + 1; i < nbins; i++) {
    if(data[i] < y) {
      dt += (data[i-1] - y) / (data[i-1] - data[i]);
      break;
    }
    else {
      dt += 1;
    }
  }
  for(int i = ta - 1; i >= 0; i--) {
    if(data[i] < y) {
      dt += (data[i+1] - y) / (data[i+1] - data[i]);
      break;
    }
    else {
      dt += 1;
    }
  }

  return dt;
}

void SetTimeParameters(BRawFADCSample* ch, BExtractedImpulse* imp, double pedestal, double posThres)
{
  Float_t a  = imp->GetA();
  // cout << a << " " << pedestal << endl;

  Float_t ymark = a / 2 + pedestal;
  Float_t t = X(ch, imp, ymark);
  Float_t fwhm = dT(ch, imp, ymark,pedestal);

  // front analysis
    Float_t y1 = 0.1 * a + pedestal;
    Float_t y2 = 0.9 * a + pedestal;

    Float_t x1 = X(ch, imp, y1);
    Float_t x2 = X(ch, imp, y2);

    Float_t tfront = x2 - x1;
    imp->SetTimeParameters(t + ch->GetOffset(), fwhm, tfront);
    // TimeOverThreshold(ch,imp,posThres);
}

Int_t ProcessPositivePulse(BRawFADCSample* ch, int &t1, Double_t posThres, Double_t initPed, Double_t initPedRMS)
{
    v_foundDPs.clear();
    // cout << "Positive pulse" << endl;
    Int_t nbins = ch->GetNbins(); // nSamples in  the waveform
    Short_t* data = ch->GetData(); // array of waveform samples
    Int_t offset = ch->GetOffset(); // time offset in the 1024 window

    Float_t q = 0, a = 0; // q = charge of the pulse, a = amplitude of the pulse
    // t1,2 = time of the first and last point above threshold, newT1,2 = the same with new extended method
    Int_t t2 = 0, ta = 0, newT2 = 0, newT1 = t1;

    // if (initPedRMS > 2)
    // 	DrawWaveform(data,nbins,initPed,offset,0,0,0,0,0,0);

    // add also a charge of the pulse before the 4sigma of the pedestal was exceeded
    q += ChargeCorrBefore(data,initPed,newT1);

    // integrate charge of the pulse, find amplitude and time of the amplitude
    for (int i = t1; i < nbins; ++i)
    {
        if (data[i] >= posThres)
        {
            double dataNorm = data[i]-initPed;
            q += dataNorm;
            if (dataNorm > a)
            {
                a = dataNorm;
                ta = i;
            }
            if (i == nbins - 1)
            {
                t2 = i;
                newT2 = i;
            }
        }else
        {
          newT2 = i;
            q += ChargeCorrAfter(data,initPed,newT2,posThres,nbins);
            t2 = i - 1;
            break;
        }
    }

    // get rid off the pulses with just one sample
    if (t1 == t2)
        return 0;

    // create and fill output class
    BExtractedImpulse* imp = new BExtractedImpulse();
    imp->SetAmpParameters(q, a);
    imp->SetTimeParameters(newT1, newT2, ta);
    SetTimeParameters(ch,imp,initPed,posThres);
    imp->SetTimeParameters(newT1+offset, newT2+offset, ta+offset);

    // search for DP
    Int_t nDPfound = FindDP(data,newT1,newT2,initPed,nbins);
    // if at least one DP found, fit the waveform with DP function otherwise fit with single
    if (nDPfound > 0)
    {
      // if (nDPfound > 10)
        cout << nDPfound << endl;
      // imp->SetDP(DPfound);//SeparatePulses();
    // imp->SetGumpbelParameters(0,0,0,0,0,0);
    // }else{
    	// TF1* fitFunc = FitDPEqual(data,initPed,newT1,newT2);
        TF1* fitFunc = FitDP(data,initPed,newT1,newT2);
    	// DrawWaveform(data,nbins,initPed,0,t1,newT1,t2,newT2,a,posThres,fitFunc);
    	delete fitFunc;
    }else
    {
      // if (nDPfound == 0 && a > 1500)
      {
        double beta = 0;
        TF1* fitFunc = FitPulse(data,newT1,newT2,initPed,a,ta,beta);
        // cout << offset << endl;
        // cout << "Real time: " << imp->GetT() << " real charge: " << imp->GetQ() << endl;
        // DrawWaveform(data,nbins,initPed,0,t1,newT1,t2,newT2,a,posThres,fitFunc);
      }
    }

    // if (posThres-initPed > a/2)
    // if (DPfound)
        // int option = DrawWaveform(data,nbins,initPed,0,t1,newT1,t2,newT2,a,posThres);
    t1 = t2;
    return 1;
}

// Method for pulse extraction. It transform raw BRawFADCSamples to BExtractedImpulseTel
Int_t AnalyzeFADCSample(BRawFADCSample* ch)
{
    //Set parameters which are waveform shape independent and the same for all pulses in the waveform
    // Int_t nchgeom = fGeomTel->GetNchGeom(ch->GetNum(), fRawMasterHeader->GetSdc());
    Int_t nfilter = ch->GetNfilter();
    Int_t nbins = ch->GetNbins(); // nSamples in the waveform
    Short_t* data = ch->GetData(); //array of sampled values of waveform
    Short_t offset = ch->GetOffset(); // offset position in the 1024 FADC time window

    // if (offset != 0)
    	// return 0;

    //Calculate initial values for pedestal and pedestalRMS
    Double_t initPed;
    Double_t initPedRMS;
    InitialPed(ch,initPed,initPedRMS,0,fSizeAvgWindowB);
    h_initPedRMS->Fill(initPedRMS);
    h_begEndDiff->Fill(data[0]-data[nbins-1]);
    // if (data[0]-data[nbins-1] > 2 && offset == 0)
    // if (initPedRMS > 2)
    	// DrawWaveform(data,nbins,initPed,offset,0,0,0,0,0,0);

    // cout << initPed << " " << initPedRMS << endl;

    // For a very long waveforms, the pedestal value can float a lot. The current value of the pedestal
    // has to be calculated in the waveform processing. The samples are used for the calculation only
    // when they are not in the region of positive or negative pulse
    Double_t currPed = initPed;
    Double_t currPedRMS = initPedRMS;
    Int_t nPedBins = 1;

    for (int i = 0; i < nbins; ++i)
    {
        Double_t negThres = (currPedRMS > 1) ? currPed - fNSigmas * currPedRMS : currPed - fNSigmas * 1;
        Double_t posThres = (currPedRMS > 1) ? currPed + fNSigmas * currPedRMS : currPed + fNSigmas * 1;

        // if the current sample is above the threshold defined for positive pulses, process it
        // otherwise use the value to update the current pedestal value
        if (data[i] > posThres) //
        {
            if (ProcessPositivePulse(ch,i,posThres,currPed,initPedRMS))
            {
              // fImpulses->At(fImpulses->GetNimpulse()-1)->SetPedParameters(initPed,initPedRMS,currPed,currPedRMS);
              // fImpulses->At(fImpulses->GetNimpulse()-1)->SetNfilter(nfilter);
            }
        }else if (data[i] < negThres)
        {
            // ProcessNegativePulse(ch,i,negThres,currPed);
        }else // Pedestal region
        {
            currPed = (currPed * nPedBins + data[i]) / (nPedBins + 1);
            currPedRMS = (currPedRMS * currPedRMS * nPedBins + (data[i] - currPed) * (data[i] - currPed)) / (nPedBins + 1);
            currPedRMS = TMath::Sqrt(currPedRMS);
            nPedBins++;
        }
    }

    return 1;
}

// Main function responsible for the reading of the files with raw waveforms, waveform processing and results drawing
int pulseExtraction(TString fileName)
{
    //open the right file and set the proper variables in the TTree
	TFile *file = new TFile(fileName.Data(), "READ");
	TTree *tree = (TTree*)file->Get("Events");

    TString outputPath = fileName(0,fileName.Last('/')+1);

    outputFile = new TFile(Form("%swaveforms.root",outputPath.Data()), "RECREATE");

	BRawMasterData *nsdmaster = 0;
	tree->SetBranchAddress("BRawMasterData.", &nsdmaster);
	BRawMasterHeader *nsdheader = 0;
	tree->SetBranchAddress("BRawMasterHeader.", &nsdheader);

    // initialize the BDTs
	InitializeBDT();

	// loop over all the events
    for (int i = 0; i < tree->GetEntries(); ++i)
    {
        // if (i%1000 == 0)
            // cout << "Event: " << i << " processed!" << endl;

        tree->GetEntry(i);
        // cout << "Event: " << i << " with: " << nsdmaster->GetNumSamples() << " raw waveforms" << endl;
        // loop over all the waveforms in given events
        for (int j = 0; j < nsdmaster->GetNumSamples(); ++j)
        {
            // cout << "Waveform: " << j << " length: " << nsdmaster->GetFADCSample(j)->GetNbins() << " offset: " << nsdmaster->GetFADCSample(j)->GetOffset() << endl;
            // badPulses += ProcessHit(nsdmaster->GetFADCSample(j),i,j,processedPulses);
            // cout << "Event: " << i << " waveform: " << j << endl;
            AnalyzeFADCSample(nsdmaster->GetFADCSample(j));
        }
    }

    DrawResults();

	return 0;
}