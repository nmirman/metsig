#include "TFile.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"
#include "TKey.h"

#include <iostream>

using namespace std;

void plotSyst(){

   TFile* fsmeared = new TFile("plots20130930/plotsWenu.root");

   TFile* fcent = new TFile("plots20130930/plotsWenu_nosmear.root");
   TFile* fup = new TFile("plots20130930/plotsWenu_nosmear_up.root");
   TFile* fdown = new TFile("plots20130930/plotsWenu_nosmear_down.root");

   TFile* fout = new TFile("plots20130930/plotsWenuJEC_nosmear.root","RECREATE");
   fout->cd();

   // loop through objects in file
   TIter nextkey(fcent->GetListOfKeys());
   TKey *key;
   while( (key = (TKey*)nextkey()) ){
      string name = key->GetName();
      string classname = key->GetClassName();

      cout << "Opening key, class " << name << ", " << classname << endl;
      if( classname.compare("TCanvas") != 0 ) continue;
      if( name.compare("pchi2_old") == 0 ) continue;

      // get canvases
      TCanvas* cmet = (TCanvas*)fcent->Get( name.c_str() );
      TCanvas* cmet_up = (TCanvas*)fup->Get( name.c_str() );
      TCanvas* cmet_down = (TCanvas*)fdown->Get( name.c_str() );
      TCanvas* cmet_smeared = (TCanvas*)fsmeared->Get( name.c_str() );

      bool containsratio = cmet->GetListOfPrimitives()->Contains("pad2");
      if( !containsratio ) continue;
      string objclass = cmet->FindObject( (name+"_MC").c_str() )->ClassName();
      if( objclass.compare("TH1D") != 0 ) continue;

      // get histograms
      TH1D* hmet_mc = (TH1D*)cmet->FindObject( (name+"_MC").c_str() );
      TH1D* hmet_data = (TH1D*)cmet->FindObject( (name+"_Data").c_str() );
      TH1D* hmet_mc_up = (TH1D*)cmet_up->FindObject( (name+"_MC").c_str() );
      TH1D* hmet_mc_down = (TH1D*)cmet_down->FindObject( (name+"_MC").c_str() );
      TH1D* hmet_mc_smeared = (TH1D*)cmet_smeared->FindObject( (name+"_MC").c_str() );

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

         gmetratio->SetPoint(i-1,
               hmetratio_smeared->GetBinCenter(i), hmetratio_smeared->GetBinContent(i));
         gmetratio->SetPointError(i-1, 0, 0, elow, ehigh);

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
      cmet_smeared->Write();

   }

   return;
}
