#include "METSigFit.h"

#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include <iostream>
#include <map>

int main(){

   // declarations
   Fitter fitter;
   std::vector<event> eventvec_MC;
   std::vector<event> eventvec_data;
   eventvec_MC.reserve( 100000 );
   eventvec_data.reserve( 1000000 );

   // fill eventvecs
   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Zmumu/Zmumu_MC_DYJettoLL_TuneZ2_M-50_7TeV_madgraph_tauola_20121104.root",
         eventvec_MC, true);
   fitter.MatchMCjets( eventvec_MC );

   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Zmumu/Zmumu_data_DoubleMu_Run2011A_08Nov2011_v1_20121104.root",
         eventvec_data, false);

   std::cout << "SIZE MC = " << eventvec_MC.size() << std::endl;
   std::cout << "SIZE Data = " << eventvec_data.size() << std::endl;
   // minimize
   std::cout << " ########### MC ########### " << std::endl;
   fitter.RunMinimizer( eventvec_MC );
   std::cout << " ########### Data ########### " << std::endl;
   fitter.RunMinimizer( eventvec_data );
   //fitter.FindSignificance( fitter.parameter, eventvec_data );


   // pileup reweighting
   int verts_MC [100] = {0};
   int verts_Data [100] = {0};
   double weights_MC [100];
   double weights_Data [100];
   std::fill_n(weights_MC,100,1.0);
   std::fill_n(weights_Data,100,1.0);
   for( std::vector<event>::iterator ev = eventvec_MC.begin(); ev < eventvec_MC.end(); ev++){
      verts_MC[ ev->nvertices ] += 1;
   }
   for( std::vector<event>::iterator ev = eventvec_data.begin(); ev < eventvec_data.end(); ev++){
      verts_Data[ ev->nvertices ] += 1;
   }
   for(int i=0; i < 100; i++){
      if( verts_MC[i] != 0){
         weights_MC[i] = (1.0*eventvec_MC.size()/eventvec_data.size())
            * (1.0*verts_Data[i]/verts_MC[i]);
      }
   }

   // plots
   std::map<std::string, TH1*> histsMC_;
   std::map<std::string, TH1*> histsMCnoPR_;
   std::map<std::string, TH1*> histsData_;

   histsMC_["muon_pt"] = new TH1D("muon_pt_MC", "Muon p_{T}", 100, 0, 200);
   histsMC_["muon_invmass"] = new TH1D("muon_invmass_MC", "M_{#mu#mu}", 100, 0, 200);
   histsMC_["jet_pt" ] = new TH1D("jet_pt_MC", "Jet p_{T}", 100, 0, 200);
   histsMC_["jet1_pt"] = new TH1D("jet1_pt_MC", "Jet 1 p_{T}", 100, 0, 200);
   histsMC_["jet2_pt"] = new TH1D("jet2_pt_MC", "Jet 2 p_{T}", 100, 0, 200);
   histsMC_["jet3_pt"] = new TH1D("jet3_pt_MC", "Jet 3 p_{T}", 100, 0, 200);
   histsMC_["njets"  ] = new TH1D("njets_MC", "N jets", 100, 0, 100);
   histsMC_["pjet_size"  ] = new TH1D("pjet_size_MC", "N jets", 100, 0, 500);
   histsMC_["pjet_scalpt"] = new TH1D("pjet_scalpt_MC", "Pseudojet Scalar p_{T}", 100, 0, 500);
   histsMC_["pjet_vectpt"] = new TH1D("pjet_vectpt_MC", "Pseudojet Scalar p_{T}", 100, 0, 200);
   histsMC_["qt"] = new TH1D("qt_MC", "q_{T}", 100, 0, 200);
   histsMC_["ut_par"] = new TH1D("ut_par_MC", "(u_{T})_{#parallel}", 100, 0, 100);
   histsMC_["nvert"] = new TH1D("nvert_MC", "N Vertices", 100, 0, 100);
   histsMC_["cov_xx"] = new TH1D("cov_xx_MC", "Cov_{xx}", 100, 0, 300);
   histsMC_["cov_xy"] = new TH1D("cov_xy_MC", "Cov_{xy}", 100, -100, 100);
   histsMC_["cov_yy"] = new TH1D("cov_yy_MC", "Cov_{yy}", 100, 0, 300);
   histsMC_["sig"] = new TH1D("sig_MC", "Significance", 100, 0, 50);
   histsMC_["det"] = new TH1D("det_MC", "Determinant", 100, 0, 100000);
   histsMC_["pchi2"] = new TH1D("pchi2_MC", "P(#chi^{2})", 100, 0, 1);

   histsMCnoPR_["muon_pt"] = new TH1D("muon_pt_MCnoPR", "Muon p_{T}", 100, 0, 200);
   histsMCnoPR_["muon_invmass"] = new TH1D("muon_invmass_MCnoPR", "M_{#mu#mu}", 100, 0, 200);
   histsMCnoPR_["jet_pt" ] = new TH1D("jet_pt_MCnoPR", "Jet p_{T}", 100, 0, 200);
   histsMCnoPR_["jet1_pt"] = new TH1D("jet1_pt_MCnoPR", "Jet 1 p_{T}", 100, 0, 200);
   histsMCnoPR_["jet2_pt"] = new TH1D("jet2_pt_MCnoPR", "Jet 2 p_{T}", 100, 0, 200);
   histsMCnoPR_["jet3_pt"] = new TH1D("jet3_pt_MCnoPR", "Jet 3 p_{T}", 100, 0, 200);
   histsMCnoPR_["njets"  ] = new TH1D("njets_MCnoPR", "N jets", 100, 0, 100);
   histsMCnoPR_["pjet_size"  ] = new TH1D("pjet_size_MCnoPR", "N jets", 100, 0, 500);
   histsMCnoPR_["pjet_scalpt"] = new TH1D("pjet_scalpt_MCnoPR", "Pseudojet Scalar p_{T}", 100, 0, 500);
   histsMCnoPR_["pjet_vectpt"] = new TH1D("pjet_vectpt_MCnoPR", "Pseudojet Scalar p_{T}", 100, 0, 200);
   histsMCnoPR_["qt"] = new TH1D("qt_MCnoPR", "q_{T}", 100, 0, 200);
   histsMCnoPR_["ut_par"] = new TH1D("ut_par_MCnoPR", "(u_{T})_{#parallel}", 100, 0, 100);
   histsMCnoPR_["nvert"] = new TH1D("nvert_MCnoPR", "N Vertices", 100, 0, 100);
   histsMCnoPR_["cov_xx"] = new TH1D("cov_xx_MCnoPR", "Cov_{xx}", 100, 0, 300);
   histsMCnoPR_["cov_xy"] = new TH1D("cov_xy_MCnoPR", "Cov_{xy}", 100, -100, 100);
   histsMCnoPR_["cov_yy"] = new TH1D("cov_yy_MCnoPR", "Cov_{yy}", 100, 0, 300);
   histsMCnoPR_["sig"] = new TH1D("sig_MCnoPR", "Significance", 100, 0, 50);
   histsMCnoPR_["det"] = new TH1D("det_MCnoPR", "Determinant", 100, 0, 100000);
   histsMCnoPR_["pchi2"] = new TH1D("pchi2_MCnoPR", "P(#chi^{2})", 100, 0, 1);

   histsData_["muon_pt"] = new TH1D("muon_pt_Data", "Muon p_{T}", 100, 0, 200);
   histsData_["muon_invmass"] = new TH1D("muon_invmass_Data", "M_{#mu#mu}", 100, 0, 200);
   histsData_["jet_pt" ] = new TH1D("jet_pt_Data", "Jet p_{T}", 100, 0, 200);
   histsData_["jet1_pt"] = new TH1D("jet1_pt_Data", "Jet 1 p_{T}", 100, 0, 200);
   histsData_["jet2_pt"] = new TH1D("jet2_pt_Data", "Jet 2 p_{T}", 100, 0, 200);
   histsData_["jet3_pt"] = new TH1D("jet3_pt_Data", "Jet 3 p_{T}", 100, 0, 200);
   histsData_["njets"  ] = new TH1D("njets_Data", "N jets", 100, 0, 100);
   histsData_["pjet_size"  ] = new TH1D("pjet_size_Data", "N jets", 100, 0, 500);
   histsData_["pjet_scalpt"] = new TH1D("pjet_scalpt_Data", "Pseudojet Scalar p_{T}", 100, 0, 500);
   histsData_["pjet_vectpt"] = new TH1D("pjet_vectpt_Data", "Pseudojet Scalar p_{T}", 100, 0, 200);
   histsData_["qt"] = new TH1D("qt_Data", "q_{T}", 100, 0, 200);
   histsData_["ut_par"] = new TH1D("ut_par_Data", "(u_{T})_{#parallel}", 100, 0, 100);
   histsData_["nvert"] = new TH1D("nvert_Data", "N Vertices", 100, 0, 100);
   histsData_["cov_xx"] = new TH1D("cov_xx_Data", "Cov_{xx}", 100, 0, 500);
   histsData_["cov_xy"] = new TH1D("cov_xy_Data", "Cov_{xy}", 100, -150, 150);
   histsData_["cov_yy"] = new TH1D("cov_yy_Data", "Cov_{yy}", 100, 0, 500);
   histsData_["sig"] = new TH1D("sig_Data", "Significance", 100, 0, 50);
   histsData_["det"] = new TH1D("det_Data", "Determinant", 100, 0, 100000);
   histsData_["pchi2"] = new TH1D("pchi2_Data", "P(#chi^{2})", 100, 0, 1);

   // fill hists
   for( int i=0; i < 3; i++ ){

      std::map<std::string, TH1*> hists_;
      std::vector<event>::iterator iter_begin;
      std::vector<event>::iterator iter_end;
      double *weights;

      if( i==0 ){
         hists_ = histsMC_;
         iter_begin = eventvec_MC.begin();
         iter_end = eventvec_MC.end();
         weights = weights_MC;
      }
      if( i==1 ){
         hists_ = histsData_;
         iter_begin = eventvec_data.begin();
         iter_end = eventvec_data.end();
         weights = weights_Data;
      }
      if( i==2 ){
         hists_ = histsMCnoPR_;
         iter_begin = eventvec_MC.begin();
         iter_end = eventvec_MC.end();
         weights = weights_Data;
      }

      for( std::vector<event>::iterator ev = iter_begin; ev < iter_end; ev++ ){
         int nvert = ev->nvertices;

         // muons
         for( int j=0; j < int(ev->muon_pt.size()); j++){
            hists_["muon_pt"]->Fill( ev->muon_pt[j] , weights[nvert]);
         }
         hists_["muon_invmass"]->Fill( ((ev->muon_4vect[0])+(ev->muon_4vect[1])).M(), 
               weights[nvert] );
         // jets
         hists_["njets"]->Fill( ev->jet_ptL123.size() , weights[nvert] );
         for( int j=0; j < int(ev->jet_ptL123.size()); j++){
            hists_["jet_pt"]->Fill( ev->jet_ptL123[j] , weights[nvert] );
         }
         if( ev->jet_ptL123.size() > 0 )
            hists_["jet1_pt"]->Fill( ev->jet_ptL123[0] , weights[nvert] );
         if( ev->jet_ptL123.size() > 1 )
            hists_["jet2_pt"]->Fill( ev->jet_ptL123[1] , weights[nvert] );
         if( ev->jet_ptL123.size() > 2 )
            hists_["jet3_pt"]->Fill( ev->jet_ptL123[2] , weights[nvert] );

         hists_["pjet_scalpt"]->Fill( ev->pjet_scalpt , weights[nvert] );
         hists_["pjet_vectpt"]->Fill( ev->pjet_vectpt , weights[nvert] );
         hists_["pjet_size"]->Fill( ev->pjet_size , weights[nvert] );

         hists_["qt"]->Fill( ev->qt, weights[nvert] );
         hists_["ut_par"]->Fill( ev->ut_par, weights[nvert] );

         hists_["nvert"]->Fill( nvert , weights[nvert] );

         hists_["cov_xx"]->Fill( ev->cov_xx, weights[nvert] );
         hists_["cov_xy"]->Fill( ev->cov_xy, weights[nvert] );
         hists_["cov_yy"]->Fill( ev->cov_yy, weights[nvert] );
         hists_["sig"]->Fill( ev->sig, weights[nvert] );
         hists_["det"]->Fill( ev->det, weights[nvert] );
         hists_["pchi2"]->Fill( TMath::Prob(ev->sig,2), weights[nvert] );
      }
   }

   TFile *file = new TFile("plots.root","RECREATE");
   file->cd();

   TCanvas *cmuon_pt = new TCanvas("cmuon_pt","cmuon_pt",700,700);
   cmuon_pt->cd();
   histsMC_["muon_pt"]->SetLineColor(2);
   histsMC_["muon_pt"]->Scale( 1.0/eventvec_MC.size() );
   histsMC_["muon_pt"]->Draw();
   histsMCnoPR_["muon_pt"]->SetLineColor(4);
   histsMCnoPR_["muon_pt"]->Scale( 1.0/eventvec_MC.size() );
   histsMCnoPR_["muon_pt"]->Draw("same");
   histsData_["muon_pt"]->SetLineColor(1);
   histsData_["muon_pt"]->SetMarkerStyle(20);
   histsData_["muon_pt"]->Sumw2();
   histsData_["muon_pt"]->Scale( 1.0/eventvec_data.size() );
   histsData_["muon_pt"]->Draw("EP same");
   cmuon_pt->Write();

   TCanvas *cmuon_invmass = new TCanvas("cmuon_invmass","cmuon_invmass",700,700);
   cmuon_invmass->cd();
   histsMC_["muon_invmass"]->SetLineColor(2);
   histsMC_["muon_invmass"]->Scale( 1.0/eventvec_MC.size() );
   histsMC_["muon_invmass"]->Draw();
   histsMCnoPR_["muon_invmass"]->SetLineColor(4);
   histsMCnoPR_["muon_invmass"]->Scale( 1.0/eventvec_MC.size() );
   histsMCnoPR_["muon_invmass"]->Draw("same");
   histsData_["muon_invmass"]->SetLineColor(1);
   histsData_["muon_invmass"]->SetMarkerStyle(20);
   histsData_["muon_invmass"]->Sumw2();
   histsData_["muon_invmass"]->Scale( 1.0/eventvec_data.size() );
   histsData_["muon_invmass"]->Draw("EP same");
   cmuon_invmass->Write();

   TCanvas *cjet_pt = new TCanvas("cjet_pt","cjet_pt",700,700);
   cjet_pt->cd();
   histsMC_["jet_pt"]->SetLineColor(2);
   histsMC_["jet_pt"]->Scale( 1.0/eventvec_MC.size() );
   histsMC_["jet_pt"]->Draw();
   histsMCnoPR_["jet_pt"]->SetLineColor(4);
   histsMCnoPR_["jet_pt"]->Scale( 1.0/eventvec_MC.size() );
   histsMCnoPR_["jet_pt"]->Draw("same");
   histsData_["jet_pt"]->SetLineColor(1);
   histsData_["jet_pt"]->SetMarkerStyle(20);
   histsData_["jet_pt"]->Sumw2();
   histsData_["jet_pt"]->Scale( 1.0/eventvec_data.size() );
   histsData_["jet_pt"]->Draw("EP same");
   cjet_pt->Write();

   TCanvas *cjet1_pt = new TCanvas("cjet1_pt","cjet1_pt",700,700);
   cjet1_pt->cd();
   histsMC_["jet1_pt"]->SetLineColor(2);
   histsMC_["jet1_pt"]->Scale( 1.0/eventvec_MC.size() );
   histsMC_["jet1_pt"]->Draw();
   histsMCnoPR_["jet1_pt"]->SetLineColor(4);
   histsMCnoPR_["jet1_pt"]->Scale( 1.0/eventvec_MC.size() );
   histsMCnoPR_["jet1_pt"]->Draw("same");
   histsData_["jet1_pt"]->SetLineColor(1);
   histsData_["jet1_pt"]->SetMarkerStyle(20);
   histsData_["jet1_pt"]->Sumw2();
   histsData_["jet1_pt"]->Scale( 1.0/eventvec_data.size() );
   histsData_["jet1_pt"]->Draw("EP same");
   cjet1_pt->Write();

   TCanvas *cjet2_pt = new TCanvas("cjet2_pt","cjet2_pt",700,700);
   cjet2_pt->cd();
   histsMC_["jet2_pt"]->SetLineColor(2);
   histsMC_["jet2_pt"]->Scale( 1.0/eventvec_MC.size() );
   histsMC_["jet2_pt"]->Draw();
   histsMCnoPR_["jet2_pt"]->SetLineColor(4);
   histsMCnoPR_["jet2_pt"]->Scale( 1.0/eventvec_MC.size() );
   histsMCnoPR_["jet2_pt"]->Draw("same");
   histsData_["jet2_pt"]->SetLineColor(1);
   histsData_["jet2_pt"]->SetMarkerStyle(20);
   histsData_["jet2_pt"]->Sumw2();
   histsData_["jet2_pt"]->Scale( 1.0/eventvec_data.size() );
   histsData_["jet2_pt"]->Draw("EP same");
   cjet2_pt->Write();

   TCanvas *cjet3_pt = new TCanvas("cjet3_pt","cjet3_pt",700,700);
   cjet3_pt->cd();
   histsMC_["jet3_pt"]->SetLineColor(2);
   histsMC_["jet3_pt"]->Scale( 1.0/eventvec_MC.size() );
   histsMC_["jet3_pt"]->Draw();
   histsMCnoPR_["jet3_pt"]->SetLineColor(4);
   histsMCnoPR_["jet3_pt"]->Scale( 1.0/eventvec_MC.size() );
   histsMCnoPR_["jet3_pt"]->Draw("same");
   histsData_["jet3_pt"]->SetLineColor(1);
   histsData_["jet3_pt"]->SetMarkerStyle(20);
   histsData_["jet3_pt"]->Sumw2();
   histsData_["jet3_pt"]->Scale( 1.0/eventvec_data.size() );
   histsData_["jet3_pt"]->Draw("EP same");
   cjet3_pt->Write();

   TCanvas *cnjets = new TCanvas("cnjets","cnjets",700,700);
   cnjets->cd();
   histsMC_["njets"]->SetLineColor(2);
   histsMC_["njets"]->Scale( 1.0/eventvec_MC.size() );
   histsMC_["njets"]->Draw();
   histsMCnoPR_["njets"]->SetLineColor(4);
   histsMCnoPR_["njets"]->Scale( 1.0/eventvec_MC.size() );
   histsMCnoPR_["njets"]->Draw("same");
   histsData_["njets"]->SetLineColor(1);
   histsData_["njets"]->SetMarkerStyle(20);
   histsData_["njets"]->Sumw2();
   histsData_["njets"]->Scale( 1.0/eventvec_data.size() );
   histsData_["njets"]->Draw("EP same");
   cnjets->Write();

   TCanvas *cpjet_scalpt = new TCanvas("cpjet_scalpt","cpjet_scalpt",700,700);
   cpjet_scalpt->cd();
   histsMC_["pjet_scalpt"]->SetLineColor(2);
   histsMC_["pjet_scalpt"]->Scale( 1.0/eventvec_MC.size() );
   histsMC_["pjet_scalpt"]->Draw();
   histsMCnoPR_["pjet_scalpt"]->SetLineColor(4);
   histsMCnoPR_["pjet_scalpt"]->Scale( 1.0/eventvec_MC.size() );
   histsMCnoPR_["pjet_scalpt"]->Draw("same");
   histsData_["pjet_scalpt"]->SetLineColor(1);
   histsData_["pjet_scalpt"]->SetMarkerStyle(20);
   histsData_["pjet_scalpt"]->Sumw2();
   histsData_["pjet_scalpt"]->Scale( 1.0/eventvec_data.size() );
   histsData_["pjet_scalpt"]->Draw("EP same");
   cpjet_scalpt->Write();

   TCanvas *cpjet_vectpt = new TCanvas("cpjet_vectpt","cpjet_vectpt",700,700);
   cpjet_vectpt->cd();
   histsMC_["pjet_vectpt"]->SetLineColor(2);
   histsMC_["pjet_vectpt"]->Scale( 1.0/eventvec_MC.size() );
   histsMC_["pjet_vectpt"]->Draw();
   histsMCnoPR_["pjet_vectpt"]->SetLineColor(4);
   histsMCnoPR_["pjet_vectpt"]->Scale( 1.0/eventvec_MC.size() );
   histsMCnoPR_["pjet_vectpt"]->Draw("same");
   histsData_["pjet_vectpt"]->SetLineColor(1);
   histsData_["pjet_vectpt"]->SetMarkerStyle(20);
   histsData_["pjet_vectpt"]->Sumw2();
   histsData_["pjet_vectpt"]->Scale( 1.0/eventvec_data.size() );
   histsData_["pjet_vectpt"]->Draw("EP same");
   cpjet_vectpt->Write();

   TCanvas *cpjet_size = new TCanvas("cpjet_size","cpjet_size",700,700);
   cpjet_size->cd();
   histsMC_["pjet_size"]->SetLineColor(2);
   histsMC_["pjet_size"]->Scale( 1.0/eventvec_MC.size() );
   histsMC_["pjet_size"]->Draw();
   histsMCnoPR_["pjet_size"]->SetLineColor(4);
   histsMCnoPR_["pjet_size"]->Scale( 1.0/eventvec_MC.size() );
   histsMCnoPR_["pjet_size"]->Draw("same");
   histsData_["pjet_size"]->SetLineColor(1);
   histsData_["pjet_size"]->SetMarkerStyle(20);
   histsData_["pjet_size"]->Sumw2();
   histsData_["pjet_size"]->Scale( 1.0/eventvec_data.size() );
   histsData_["pjet_size"]->Draw("EP same");
   cpjet_size->Write();

   TCanvas *cnvert = new TCanvas("cnvert","cnvert",700,700);
   cnvert->cd();
   histsMC_["nvert"]->SetLineColor(2);
   histsMC_["nvert"]->Scale( 1.0/eventvec_MC.size() );
   histsMC_["nvert"]->Draw();
   histsMCnoPR_["nvert"]->SetLineColor(4);
   histsMCnoPR_["nvert"]->Scale( 1.0/eventvec_MC.size() );
   histsMCnoPR_["nvert"]->Draw("same");
   histsData_["nvert"]->SetLineColor(1);
   histsData_["nvert"]->SetMarkerStyle(20);
   histsData_["nvert"]->Sumw2();
   histsData_["nvert"]->Scale( 1.0/eventvec_data.size() );
   histsData_["nvert"]->Draw("EP same");
   cnvert->Write();

   TCanvas *ccov_xx = new TCanvas("ccov_xx","ccov_xx",700,700);
   ccov_xx->cd();
   histsMC_["cov_xx"]->SetLineColor(2);
   histsMC_["cov_xx"]->Scale( 1.0/eventvec_MC.size() );
   histsMC_["cov_xx"]->Draw();
   histsMCnoPR_["cov_xx"]->SetLineColor(4);
   histsMCnoPR_["cov_xx"]->Scale( 1.0/eventvec_MC.size() );
   histsMCnoPR_["cov_xx"]->Draw("same");
   histsData_["cov_xx"]->SetLineColor(1);
   histsData_["cov_xx"]->SetMarkerStyle(20);
   histsData_["cov_xx"]->Sumw2();
   histsData_["cov_xx"]->Scale( 1.0/eventvec_data.size() );
   histsData_["cov_xx"]->Draw("EP same");
   ccov_xx->Write();

   TCanvas *ccov_xy = new TCanvas("ccov_xy","ccov_xy",700,700);
   ccov_xy->cd();
   histsMC_["cov_xy"]->SetLineColor(2);
   histsMC_["cov_xy"]->Scale( 1.0/eventvec_MC.size() );
   histsMC_["cov_xy"]->Draw();
   histsMCnoPR_["cov_xy"]->SetLineColor(4);
   histsMCnoPR_["cov_xy"]->Scale( 1.0/eventvec_MC.size() );
   histsMCnoPR_["cov_xy"]->Draw("same");
   histsData_["cov_xy"]->SetLineColor(1);
   histsData_["cov_xy"]->SetMarkerStyle(20);
   histsData_["cov_xy"]->Sumw2();
   histsData_["cov_xy"]->Scale( 1.0/eventvec_data.size() );
   histsData_["cov_xy"]->Draw("EP same");
   ccov_xy->Write();

   TCanvas *ccov_yy = new TCanvas("ccov_yy","ccov_yy",700,700);
   ccov_yy->cd();
   histsMC_["cov_yy"]->SetLineColor(2);
   histsMC_["cov_yy"]->Scale( 1.0/eventvec_MC.size() );
   histsMC_["cov_yy"]->Draw();
   histsMCnoPR_["cov_yy"]->SetLineColor(4);
   histsMCnoPR_["cov_yy"]->Scale( 1.0/eventvec_MC.size() );
   histsMCnoPR_["cov_yy"]->Draw("same");
   histsData_["cov_yy"]->SetLineColor(1);
   histsData_["cov_yy"]->SetMarkerStyle(20);
   histsData_["cov_yy"]->Sumw2();
   histsData_["cov_yy"]->Scale( 1.0/eventvec_data.size() );
   histsData_["cov_yy"]->Draw("EP same");
   ccov_yy->Write();

   TCanvas *csig = new TCanvas("csig","csig",700,700);
   csig->cd();
   histsMC_["sig"]->SetLineColor(2);
   histsMC_["sig"]->Scale( 1.0/eventvec_MC.size() );
   histsMC_["sig"]->Draw();
   histsMCnoPR_["sig"]->SetLineColor(4);
   histsMCnoPR_["sig"]->Scale( 1.0/eventvec_MC.size() );
   histsMCnoPR_["sig"]->Draw("same");
   histsData_["sig"]->SetLineColor(1);
   histsData_["sig"]->SetMarkerStyle(20);
   histsData_["sig"]->Sumw2();
   histsData_["sig"]->Scale( 1.0/eventvec_data.size() );
   histsData_["sig"]->Draw("EP same");
   csig->Write();

   TCanvas *cdet = new TCanvas("cdet","cdet",700,700);
   cdet->cd();
   histsMC_["det"]->SetLineColor(2);
   histsMC_["det"]->Scale( 1.0/eventvec_MC.size() );
   histsMC_["det"]->Draw();
   histsMCnoPR_["det"]->SetLineColor(4);
   histsMCnoPR_["det"]->Scale( 1.0/eventvec_MC.size() );
   histsMCnoPR_["det"]->Draw("same");
   histsData_["det"]->SetLineColor(1);
   histsData_["det"]->SetMarkerStyle(20);
   histsData_["det"]->Sumw2();
   histsData_["det"]->Scale( 1.0/eventvec_data.size() );
   histsData_["det"]->Draw("EP same");
   cdet->Write();

   TCanvas *cpchi2 = new TCanvas("cpchi2","cpchi2",700,700);
   cpchi2->cd();
   histsMC_["pchi2"]->SetLineColor(2);
   histsMC_["pchi2"]->Scale( 1.0/eventvec_MC.size() );
   histsMC_["pchi2"]->Draw();
   histsMCnoPR_["pchi2"]->SetLineColor(4);
   histsMCnoPR_["pchi2"]->Scale( 1.0/eventvec_MC.size() );
   histsMCnoPR_["pchi2"]->Draw("same");
   histsData_["pchi2"]->SetLineColor(1);
   histsData_["pchi2"]->SetMarkerStyle(20);
   histsData_["pchi2"]->Sumw2();
   histsData_["pchi2"]->Scale( 1.0/eventvec_data.size() );
   histsData_["pchi2"]->Draw("EP same");
   histsMC_["pchi2"]->SetMaximum( 0.025 );
   histsMC_["pchi2"]->SetMinimum( 0.0 );
   cpchi2->Write();

   TCanvas *cqt = new TCanvas("cqt","cqt",700,700);
   cqt->cd();
   histsMC_["qt"]->SetLineColor(2);
   histsMC_["qt"]->Scale( 1.0/eventvec_MC.size() );
   histsMC_["qt"]->Draw();
   histsMCnoPR_["qt"]->SetLineColor(4);
   histsMCnoPR_["qt"]->Scale( 1.0/eventvec_MC.size() );
   histsMCnoPR_["qt"]->Draw("same");
   histsData_["qt"]->SetLineColor(1);
   histsData_["qt"]->SetMarkerStyle(20);
   histsData_["qt"]->Sumw2();
   histsData_["qt"]->Scale( 1.0/eventvec_data.size() );
   histsData_["qt"]->Draw("EP same");
   cqt->Write();

   TCanvas *cut_par = new TCanvas("cut_par","cut_par",700,700);
   cut_par->cd();
   histsMC_["ut_par"]->SetLineColor(2);
   histsMC_["ut_par"]->Scale( 1.0/eventvec_MC.size() );
   histsMC_["ut_par"]->Draw();
   histsMCnoPR_["ut_par"]->SetLineColor(4);
   histsMCnoPR_["ut_par"]->Scale( 1.0/eventvec_MC.size() );
   histsMCnoPR_["ut_par"]->Draw("same");
   histsData_["ut_par"]->SetLineColor(1);
   histsData_["ut_par"]->SetMarkerStyle(20);
   histsData_["ut_par"]->Sumw2();
   histsData_["ut_par"]->Scale( 1.0/eventvec_data.size() );
   histsData_["ut_par"]->Draw("EP same");
   cut_par->Write();

   return 0;
}
