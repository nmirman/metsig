#include "TRandom3.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TVector2.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "METSigFit.h"

#define JETPT_LOW 10.0
#define JETPT_HIGH 30.0

Fitter::Fitter(){
}

Fitter::~Fitter(){
}

void Fitter::ReadNtuple(const char* filename, bool isMC){
   std::cout << __LINE__ << std::endl;   
   int v_size;
   float met_et;

   int pfj_size;
   float pfj_pt[1000];
   float pfj_phi[1000];
   float pfj_eta[1000];

   int mu_size;
   float mu_pt[100];
   float mu_phi[100];

   float pfj_l1[1000];
   float pfj_l1l2l3[1000];

   int genj_size;
   float genj_pt[1000];
   float genj_phi[1000];
   float genj_eta[1000];

   TFile *file = new TFile(filename);
   TTree *tree = (TTree*)file->Get("events");
   tree->SetBranchAddress("v_size", &v_size);
   tree->SetBranchAddress("met_et",&met_et);

   tree->SetBranchAddress("pfj_size", &pfj_size);
   tree->SetBranchAddress("pfj_pt", pfj_pt);
   tree->SetBranchAddress("pfj_phi", pfj_phi);
   tree->SetBranchAddress("pfj_eta", pfj_eta);

   tree->SetBranchAddress("mu_size", &mu_size);
   tree->SetBranchAddress("mu_pt", mu_pt);
   tree->SetBranchAddress("mu_phi", mu_phi);

   tree->SetBranchAddress("pfj_l1", &pfj_l1);
   tree->SetBranchAddress("pfj_l1l2l3", &pfj_l1l2l3);

   std::cout << __LINE__ << std::endl;   
   if(isMC){
      tree->SetBranchAddress("genj_size", &genj_size);
      tree->SetBranchAddress("genj_pt", genj_pt);
      tree->SetBranchAddress("genj_phi", genj_phi);
      tree->SetBranchAddress("genj_eta", genj_eta);
   }

   std::cout << __LINE__ << std::endl;   
   for( int ev=0; ev<tree->GetEntries(); ev++){
      tree->GetEntry(ev);
      if( mu_size != 2 or mu_pt[0] < 25 or mu_pt[1] < 20) continue;

      double pjet_ht_temp = 0;
      double pjet_px_temp = 0;
      double pjet_py_temp = 0;

      event evtemp;

      // vertices
      evtemp.nvertices = v_size;
      // muons
      for( int i=0; i < mu_size; i++){
         evtemp.muon_pt.push_back( mu_pt[i] );
         evtemp.muon_phi.push_back( mu_phi[i] );
      }
      // jets
      for( int i=0; i < pfj_size; i++){
         
         if( pfj_pt[i] > JETPT_LOW ){
            // clustered jets
            
            evtemp.jet_pt.push_back( pfj_pt[i] );
            evtemp.jet_phi.push_back( pfj_phi[i] );
            evtemp.jet_eta.push_back( pfj_eta[i] );

            evtemp.jet_corrL123.push_back( pfj_l1l2l3[i] );
            evtemp.jet_corrL1.push_back( pfj_l1[i] );

         } else {
            // pseudojet with unclustered energy

            pjet_ht_temp += pfj_pt[i];
            pjet_px_temp += pfj_pt[i]*cos(pfj_phi[i]);
            pjet_py_temp += pfj_pt[i]*sin(pfj_phi[i]);

         }

      } // pfj loop

      // fill pseudojet quantities
      evtemp.pjet_ht = pjet_ht_temp;
      evtemp.pjet_pt = sqrt( pjet_px_temp*pjet_px_temp + pjet_py_temp*pjet_py_temp );
      evtemp.pjet_phi = atan2( pjet_px_temp, pjet_py_temp );

      // genjets
      if(isMC){
         for( int i=0; i < genj_size; i++){
            evtemp.genjet_pt.push_back( genj_pt[i] );
            evtemp.genjet_phi.push_back( genj_phi[i] );
            evtemp.genjet_eta.push_back( genj_eta[i] );
         }
      }

      // fill event vector
      if( isMC ){
         eventvec_MC.push_back( evtemp );
      }else{
         eventvec_data.push_back( evtemp );
      }

   } // event loop

   std::cout << __LINE__ << std::endl;   
      if(isMC) MatchMCjets();
   std::cout << __LINE__ << std::endl;   

   delete file;

}

void Fitter::MatchMCjets(){

   for(int ev=0; ev < int(eventvec_MC.size()); ev++){

      // loop through reco jets
      for(int ireco=0; ireco < int(eventvec_MC[ev].jet_pt.size()); ireco++){
         
         int matchIndex = -1;
         double dR = 1000;

         // loop through genjets
         for(int igen=0; igen < int(eventvec_MC[ev].genjet_pt.size()); igen++){
            
            double dphi = TVector2::Phi_mpi_pi( eventvec_MC[ev].jet_phi[ireco] 
                  - eventvec_MC[ev].genjet_phi[igen] );
            double deta = eventvec_MC[ev].jet_eta[ireco] 
               - eventvec_MC[ev].genjet_eta[igen];
            double dRtemp = sqrt( deta*deta + dphi*dphi );

            if( dRtemp < dR ){
               dR = dRtemp;
               matchIndex = igen;
            }

         } // genjets

         eventvec_MC[ev].jet_matchIndex.push_back( matchIndex );
         eventvec_MC[ev].jet_matchdR.push_back( dR );

      } // jets

      if( eventvec_MC[ev].jet_pt.size() != eventvec_MC[ev].jet_matchIndex.size() ){
         std::cout << "VECTOR SIZE MISMATCH" << std::endl;
      }

   } // event loop

}

