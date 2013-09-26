#include "TFile.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"

#include <iostream>

using namespace std;

void plotSyst(){

   //TFile* fcent = new TFile("results/plotsZmumu_nosmear.root");
   //TFile* fup = new TFile("results/plotsZmumu_mcup_nosmear.root");
   //TFile* fdown = new TFile("results/plotsZmumu_mcdown_nosmear.root");
   TFile* fcent = new TFile("plots20130926/plotsZmumu.root");
   TFile* fup = new TFile("plots20130926/plotsZmumu_up.root");
   TFile* fdown = new TFile("plots20130926/plotsZmumu_down.root");

   TFile* fsmeared = new TFile("results/plotsZmumu.root");

   // get canvases
   TCanvas* cmet = (TCanvas*)fcent->Get("sig_100");
   TCanvas* cmet_up = (TCanvas*)fup->Get("sig_100");
   TCanvas* cmet_down = (TCanvas*)fdown->Get("sig_100");
   TCanvas* cmet_smeared = (TCanvas*)fsmeared->Get("sig_100");
   
   // get histograms
   TH1D* hmet_mc = (TH1D*)cmet->FindObject("sig_100_MC");
   TH1D* hmet_data = (TH1D*)cmet->FindObject("sig_100_Data");
   TH1D* hmet_mc_up = (TH1D*)cmet_up->FindObject("sig_100_MC");
   TH1D* hmet_mc_down = (TH1D*)cmet_down->FindObject("sig_100_MC");
   TH1D* hmet_mc_smeared = (TH1D*)cmet_smeared->FindObject("sig_100_MC");

   TH1D* hmetratio = (TH1D*)cmet->FindObject("hratio");
   TH1D* hmetratio_up = (TH1D*)cmet_up->FindObject("hratio");
   TH1D* hmetratio_down = (TH1D*)cmet_down->FindObject("hratio");
   TH1D* hmetratio_smeared = (TH1D*)cmet_smeared->FindObject("hratio");

   // compute error bands
   TGraphAsymmErrors* gmet = new TGraphAsymmErrors();
   for(int i=1; i <= hmet_mc->GetNbinsX(); i++){

      double maxvar = max(max(hmet_mc_up->GetBinContent(i),hmet_mc_down->GetBinContent(i)),
            hmet_mc->GetBinContent(i));
      double minvar = min(min(hmet_mc_up->GetBinContent(i),hmet_mc_down->GetBinContent(i)),
            hmet_mc->GetBinContent(i));
      double ehigh = maxvar - hmet_mc->GetBinContent(i);
      double elow = hmet_mc->GetBinContent(i) - minvar;

      gmet->SetPoint(i-1, hmet_mc_smeared->GetBinCenter(i), hmet_mc_smeared->GetBinContent(i));
      gmet->SetPointError(i-1, 0, 0, elow, ehigh);
      
   }

   TGraphAsymmErrors* gmetratio = new TGraphAsymmErrors();
   for(int i=1; i <= hmetratio->GetNbinsX(); i++){

      double maxvar = max(max(hmetratio_up->GetBinContent(i),hmetratio_down->GetBinContent(i)),
            hmetratio->GetBinContent(i));
      double minvar = min(min(hmetratio_up->GetBinContent(i),hmetratio_down->GetBinContent(i)),
            hmetratio->GetBinContent(i));
      double ehigh = maxvar - hmetratio->GetBinContent(i);
      double elow = hmetratio->GetBinContent(i) - minvar;

      gmetratio->SetPoint(i-1, hmetratio_smeared->GetBinCenter(i), hmetratio_smeared->GetBinContent(i));
      gmetratio->SetPointError(i-1, 0, 0, elow, ehigh);

      //cout << hmetratio->GetBinCenter(i) << ": " << hmetratio->GetBinContent(i)
      //   << " + " << maxvar << " - " << minvar << endl;
      
   }

   // add bands to plot, print to file
   TPad* p1 = (TPad*)cmet_smeared->FindObject("pad1");
   p1->cd();
   gmet->SetFillColor(17);
   gmet->SetFillStyle(3001);
   gmet->Draw("4");
   hmet_data->Draw("same");

   TPad* p2 = (TPad*)cmet_smeared->FindObject("pad2");
   p2->cd();
   gmetratio->SetFillColor(17);
   gmetratio->SetFillStyle(3001);
   gmetratio->Draw("3");
   hmetratio_smeared->Draw("same");

   cmet_smeared->cd();
   cmet_smeared->Draw();
   
   //cmet_down->SetName("down");
   //cmet_down->Draw();
   //cmet_up->SetName("up");
   //cmet_up->Draw();
   

   //TCanvas * c = new TCanvas("c","c",800,800);
   //c->cd();
   //hmetratio_up->Draw();
   //hmetratio_down->SetMarkerColor(2);
   //hmetratio_down->SetLineColor(2);
   //hmetratio_down->Draw("same");
   //c->Draw();

   return;
}
