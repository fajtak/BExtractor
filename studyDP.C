#include "BExtractedImpulseTel.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TF1.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMaterial.h"
#include "THStack.h"
#include "TGeoMedium.h"
#include "TSystem.h"
#include "TThread.h"
#include "TProfile.h"

#include <iostream>

const int c_nStrings = 8;
const int c_nOMsPerString = 36;
const double c_xPos[c_nStrings]{5.58,51.72,57.78,25.21,-29.56,-53.45,-41.94,0};
const double c_yPos[c_nStrings]{61.80,37.00,-14.9,-52.8,-52.9,-8.3,41.7,0};
int dummy;

TCanvas* myCan = new TCanvas("myCan","Results",800,600);

Bool_t VisEvent(BExtractedImpulseTel* extrImpTel)
{
    gSystem->Load("libGeom");
    gSystem->Load("/home/fajtak/work/BARS/trunk/libmars");
    gStyle->SetPalette(1, 0);

    TGeoManager *geom=new TGeoManager("geom","My first 3D geometry");

    TGeoMaterial *vacuum=new TGeoMaterial("vacuum",0,0,0);
    TGeoMaterial *Fe=new TGeoMaterial("Fe",55.845,26,7.87);

    TGeoMedium *Air=new TGeoMedium("Vacuum",0,vacuum);
    TGeoMedium *Iron=new TGeoMedium("Iron",1,Fe);

    TGeoVolume *top=geom->MakeBox("top",Air,1000,1000,1000);
    geom->SetTopVolume(top);
    geom->SetTopVisible(0);

    // draw OMs

    for(int n = 0; n < c_nOMsPerString; n++)
    {
        for (int j = 0; j < c_nStrings; ++j)
        {
            TGeoVolume* sphere = geom->MakeSphere(Form("Sphere%d", n+j*c_nOMsPerString),Iron,0,1);
            sphere->SetFillColor(0);
            sphere->SetLineColor(0);
            top->AddNodeOverlap(sphere,1,new TGeoTranslation(c_xPos[j],c_yPos[j],n*15));    
        }
    }

    // find minimal and maximal times to define color scale
    Float_t timmin =  0;
    Float_t timmax =  1024;

    geom->CloseGeometry();
    top->Draw("ogl");

    // draw hit channels
    for(int i = 0; i < extrImpTel->GetNimpulse(); i++)
    {
        // if (extrImpTel->GetQ(i) < 150 || extrImpTel->GetT(i)<350 || extrImpTel->GetT(i)>750 )
        // if (extrImpTel->GetT(i)<350 || extrImpTel->GetT(i)>750 )
            // continue;
        int OM = extrImpTel->GetNch(i);
        int zID = OM%36;
        int stringID = OM/36;
        // TGeoVolume* sphere = geom->MakeSphere(Form("SphereEvent%d", i),Iron,0,(TMath::Log(extrImpTel->GetQ(i)/150)+1)*5);
        int pointSize = extrImpTel->At(i)->GetIsDP() ? 10 : 5;
        TGeoVolume* sphere = geom->MakeSphere(Form("SphereEvent%d", i),Iron,0,pointSize);
        sphere->SetLineColor(gStyle->GetColorPalette((gStyle->GetNumberOfColors() - 1) * (1 - (extrImpTel->GetT(i) - timmin)/(timmax-timmin))));
        // sphere->SetFillColor(gStyle->GetColorPalette((gStyle->GetNumberOfColors() - 1) * (1 - (extrImpTel->GetT(i) - timmin)/(timmax-timmin))));
        top->AddNodeOverlap(sphere,1,new TGeoTranslation(c_xPos[stringID],c_yPos[stringID],zID*15));
        //cout << imp->GetChannelID() << " " << amp << " " << tim << endl;
    }

    top->SetVisibility(0);
    geom->CloseGeometry();

    top->DrawClone("ogl");
    // int dummy;
    // cin >> dummy;
    return kTRUE;
}

Bool_t VisTime(BExtractedImpulseTel* extrImpTel)
{
    int nPulses = extrImpTel->GetNimpulse();
    
    Float_t* x = new Float_t[nPulses];
    Float_t* y = new Float_t[nPulses];

    for (int i = 0; i < nPulses; ++i)
    {
        x[i] = extrImpTel->GetT(i);
        y[i] = extrImpTel->GetNch(i);
    }

    TGraph* timePlot = new TGraph(nPulses,y,x);
    timePlot->SetMarkerColor(kRed);
    timePlot->SetMarkerStyle(3);
    myCan->cd();
    timePlot->Draw("AXP");
    myCan->Update();

    return kTRUE;
}

int nDPsPerEvent(BExtractedImpulseTel* extrImpTel)
{
    int nDPs = 0;
    for (int i = 0; i < extrImpTel->GetNimpulse(); ++i)
    {
        if (extrImpTel->At(i)->GetIsDP())
            nDPs++;
    }
    return nDPs;
}

int nDPsPerEventInMiddle(BExtractedImpulseTel* extrImpTel)
{
    int nDPs = 0;
    for (int i = 0; i < extrImpTel->GetNimpulse(); ++i)
    {
        if (extrImpTel->At(i)->GetIsDP() && extrImpTel->GetT(i) > 400 && extrImpTel->GetT(i) < 700)
            nDPs++;
    }
    return nDPs;
}

int QPerEvent(BExtractedImpulseTel* extrImpTel)
{
    double Q = 0;
    for (int i = 0; i < extrImpTel->GetNimpulse(); ++i)
    {
        Q += extrImpTel->GetQ(i);
    }
    return Q;
}

void *handle(void *ptr) {
    cout << "Enter a number to continue: ";
    cin >> dummy;
    return 0;
}

TH1F* h_timeImpulses = new TH1F("h_timeImpulses","Time distribution of all impulses;Time [FADC channels];NoE [#]",1024,0,1024);
TH1F* h_timeDPs = new TH1F("h_timeDPs","Time distribution of DPs;Time [FADC channels];NoE [#]",1024,0,1024);
TH1F* h_qImpulses = new TH1F("h_qImpulses","Charge distribution of all impulses;Charge [FADC channels];NoE [#]",6000,0,60000);
TH1F* h_qDPs = new TH1F("h_qDPs","Charge distribution of DPs;Charge [FADC channels];NoE [#]",6000,0,60000);


bool TimeAnalysis(BExtractedImpulseTel* extrImpTel)
{
    for (int i = 0; i < extrImpTel->GetNimpulse(); ++i)
    {
        h_qImpulses->Fill(extrImpTel->GetQ(i));
        h_timeImpulses->Fill(extrImpTel->GetT(i));
        if (extrImpTel->At(i)->GetIsDP())
        {
            h_qDPs->Fill(extrImpTel->GetQ(i));
            h_timeDPs->Fill(extrImpTel->GetT(i));
        }
    }
    return kTRUE;
}

int studyDP()
{
    gStyle->SetOptStat(111111);
    gStyle->SetOptFit(1);

    // definition of the important pointers ---------------------------------------------------------------------------------------------
    // TFile *file = new TFile("/Data/BaikalData/2016/cluster-0/0102/f0102.joint.events.root", "READ");
    // TFile *file = new TFile("/Data/BaikalData/2018/cluster-2/0005/h0005.joint.events.root", "READ");
    TFile *file = new TFile("/Data/BaikalData/2018/cluster-2/0005_DP12/h0005.joint.events.dp12.root", "READ");
    // TFile *file = new TFile("/Data/BaikalData/2018/cluster-2/0005_DP7/h0005.joint.events.dp7.root", "READ");
    // TFile *file = new TFile("/Data/BaikalData/2018/cluster-2/0005_DP10/h0005.joint.events.dp10.root", "READ");
    // TFile *file = new TFile("/Data/BaikalData/2018/cluster-2/0006/h0006.joint.events.root", "READ");
    TTree *tree = (TTree*)file->Get("Events");

    BExtractedImpulseTel* extrImpTel = NULL;
    tree->SetBranchAddress("BJointImpulseTel.", &extrImpTel);

    TH1F* h_nImpulses = new TH1F("nImpulses","All events;Number of impulses [#];NoE[#]",500,0,500);
    TH1F* h_nImpulsesDP = new TH1F("nImpulsesDP","Events with nDP > 30;Number of impulses [#];NoE[#]",500,0,500);
    TH2F* h_nImpVsnDP = new TH2F("nImpBsnDP","Number of Impulses vs. number of DPs;Number of Impulses [#];Number of DPs [#]",500,0,500,500,0,500);
    TH1F* h_nDPs = new TH1F("nDPs","Number of DPs;Number of DPs [#];NoE[#]",500,0,500);
    TH1F* h_QAll = new TH1F("h_QAll","Overall charge per event",1000,0,1000000);
    TH1F* h_QSP = new TH1F("h_QSP","Overall charge per event with only SPs",1000,0,1000000);
    h_nDPs->SetLineColor(kGreen);
    h_nDPs->SetLineWidth(2);
    h_nImpulsesDP->SetLineColor(kGreen);
    h_nImpulsesDP->SetLineWidth(2);
    h_nImpulses->SetLineWidth(2);

    cout << "Number of entries: " << tree->GetEntries() << endl;

    int nPulsesTotal = 0;
    int nDPsTotal = 0;

    for(Int_t i = 0; i < tree->GetEntries(); i++)
    {
        if (i%100000 == 0)
            cout << i << endl;
        tree->GetEntry(i);
        TimeAnalysis(extrImpTel);
        int nPulses = extrImpTel->GetNimpulse();
        nPulsesTotal += nPulses;
        h_nImpulses->Fill(nPulses);
        int nDPs = nDPsPerEvent(extrImpTel);
        // int nDPs = nDPsPerEventInMiddle(extrImpTel);
        nDPsTotal += nDPs;
        h_nDPs->Fill(nDPs);
        h_nImpVsnDP->Fill(nPulses,nDPs);
        if (nDPs >= 30)
            h_nImpulsesDP->Fill(nPulses);
        double eventCharge = QPerEvent(extrImpTel);
        h_QAll->Fill(eventCharge);
        if (nDPs == 0)
            h_QSP->Fill(eventCharge);
        if (nPulses > 120 && nDPs > 30)
        // if (nPulses > 33 && nPulses < 53 && nDPs == 0 && eventCharge > 7500 && eventCharge < 15000)
        {
            cout << "Event No. : " << i << endl;  
            TThread *t = new TThread("t", handle, (void*) 1);
            VisEvent(extrImpTel);
            // VisTime(extrImpTel);
            // extrImpTel->Print();
            t->Run();
            t->Join();
            delete t;
        }
    }

    cout << "DPs/All pulses (%): " << nDPsTotal << "/" << nPulsesTotal << " (" << (double)nDPsTotal/nPulsesTotal*100 << "%)" << endl;
    cout << "End" << endl;

    TCanvas* c_Results1 = new TCanvas("c_Results1","Results",800,600);
    c_Results1->cd();
    THStack *hs1 = new THStack("hs1","Number of impulses & DPs;Number of pulses [#];NoE [#]");
    hs1->Add(h_nImpulses);
    hs1->Add(h_nDPs);
    hs1->Draw("nostack");
    gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

    TCanvas* c_Results2 = new TCanvas("c_Results2","Results",800,600);
    c_Results2->cd();
    THStack *hs2 = new THStack("hs2","Number of impulses per Event;Number of pulses [#];NoE [#]");
    hs2->Add(h_nImpulses);
    hs2->Add(h_nImpulsesDP);
    hs2->Draw("nostack");
    gPad->BuildLegend(0.65,0.55,0.9,0.9,"");

    TCanvas* c_Results3 = new TCanvas("c_Results3","Results",800,600);
    c_Results3->cd();
    THStack *hs3 = new THStack("hs3","Pulse time distribution;Time [FADC channels];NoE [#]");
    h_timeDPs->SetLineColor(kGreen);
    hs3->Add(h_timeImpulses);
    hs3->Add(h_timeDPs);
    hs3->Draw("nostack");
    gPad->BuildLegend(0.65,0.55,0.9,0.9,"");

    TCanvas* c_Results4 = new TCanvas("c_Results4","Results",800,600);
    c_Results4->cd();
    THStack *hs4 = new THStack("hs4","Pulse Charge distribution;Charge [FADC channels];NoE [#]");
    h_timeDPs->SetLineColor(kGreen);
    hs4->Add(h_qImpulses);
    hs4->Add(h_qDPs);
    hs4->Draw("nostack");
    gPad->BuildLegend(0.65,0.55,0.9,0.9,"");

    TCanvas* c_Results5 = new TCanvas("c_Results5","Results",800,600);
    c_Results5->cd();
    TH1F *h3 = (TH1F*)h_timeDPs->Clone("h3");
    TH1F *h4 = (TH1F*)h_timeImpulses->Clone("h4");
    h3->Rebin(4);
    h4->Rebin(4);
    h3->SetLineColor(kGreen);
    h3->Sumw2();
    // h3->SetStats(0);      // No statistics on lower plot
    h3->Divide(h3,h4,100,1);
    h3->SetTitle("Ratio DPs vs All Pulses in time;Time [FADC channels]; DPs/All pulses [%]");
    h3->SetMarkerStyle(21);
    h3->Draw("pl");       // Draw the ratio plot

    TCanvas* c_Results6 = new TCanvas("c_Results6","Results",800,600);
    c_Results6->cd();
    h_nImpVsnDP->Draw("colz");


    TProfile *newProfile = h_nImpVsnDP->ProfileX();
    TH1F* DPratio = new TH1F("DPratio","Relative number of DPs",newProfile->GetNbinsX(),0,newProfile->GetNbinsX());

    for (int i = 0; i < newProfile->GetNbinsX(); ++i)
    {
        if (newProfile->GetBinContent(i) != 0)
        {
            DPratio->Fill(i,newProfile->GetBinContent(i)/i*100);
        }
        // cout << i << " " << newProfile->GetBinCenter(i) << " " << newProfile->GetBinContent(i) << " " << newProfile->GetBinContent(i)/i*100 << endl;
    }

    TCanvas* c_Results7 = new TCanvas("c_Results7","Results",800,600);
    c_Results7->cd();
    newProfile->Draw();

    TCanvas* c_Results8 = new TCanvas("c_Results8","Results",800,600);
    c_Results8->cd();
    DPratio->Draw();

    TCanvas* c_Results9 = new TCanvas("c_Results9","Results",800,600);
    c_Results9->cd();
    h_QAll->Draw();
    h_QSP->Draw("same");

    return 1;
}