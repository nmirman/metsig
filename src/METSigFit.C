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

#include <cmath>
#include <iostream>
#include <iomanip>

//#include "Math/Minimizer.h"
//#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"

#include "METSigFit.h"

#define JETPT_LOW 10.0
#define JETPT_HIGH 20.0

const double Fitter::sigmaPt[10][4]={{-0.349206, 0.297831, 0, 0.471121},
   {-0.499735, 0.336391, 0, 0.430689},
   {-0.561649, 0.420293, 0, 0.392398},
   {-1.12329, 0.657891, 0, 0.139595},
   {1.04792, 0.466763, 0, 0.193137},
   {2.56933, 0.305802, 0, 0.398929},
   {2.81978, 0.272373, 0, 0.579396},
   {1.65966, 0.223683, 0, 0.60873},
   {1.41584, 0.209477, 0, 0.588872},
   {1.41584, 0.209477, 0, 0.588872}};

const double Fitter::sigmaPhi[10][5]={{926.978, 2.52747, 0.0304001, -926.224, -1.94117},
   {3.32512e-06,     0.063941,  -0.00387593,  0.301932,    -0.825352},
   {    0.38469,    0.0755727,   -0.0044353,  0.453887,      -1.8947},
   {2.92001e-07,    0.0718389,  -0.00385579,  0.403668,     -0.62698},
   { 0.00336639,    0.0880209,   -0.0023084,  0.214304,    -0.416353},
   {    11.1957,     0.643236,   0.00711422,  -10.7613,     0.280927},
   {     1.9027, -4.56986e-06,    0.0304793,  -1.09124,    -0.136402},
   {    2.11446,     0.203329,   -0.0175832,  -1.67946,  -0.00853474},
   {   0.765787, -3.90638e-06, -4.70224e-08,   0.11831,      -1.4675},
   {    259.189,   0.00132792,    -0.311411,  -258.647,            0}};


void Fitter::ReadNtuple(const char* filename, bool isMC){
   std::cout << "ReadNtuple -->" << std::endl;

   int v_size;
   float met_et;

   int pfj_size=0;
   float pfj_pt[1000];
   float pfj_phi[1000];
   float pfj_eta[1000];

   int mu_size=0;
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
   tree->SetBranchAddress("met_et", &met_et);

   tree->SetBranchAddress("pfj_size", &pfj_size);
   tree->SetBranchAddress("pfj_pt", pfj_pt);
   tree->SetBranchAddress("pfj_phi", pfj_phi);
   tree->SetBranchAddress("pfj_eta", pfj_eta);

   tree->SetBranchAddress("mu_size", &mu_size);
   tree->SetBranchAddress("mu_pt", mu_pt);
   tree->SetBranchAddress("mu_phi", mu_phi);

   tree->SetBranchAddress("pfj_l1", pfj_l1);
   tree->SetBranchAddress("pfj_l1l2l3", pfj_l1l2l3);

   if(isMC){
      tree->SetBranchAddress("genj_size", &genj_size);
      tree->SetBranchAddress("genj_pt", genj_pt);
      tree->SetBranchAddress("genj_phi", genj_phi);
      tree->SetBranchAddress("genj_eta", genj_eta);
   }

   std::cout << "  fill event vector" << std::endl;

   for( int ev=0; ev<tree->GetEntries(); ev++){
      tree->GetEntry(ev);
      if( mu_size != 2 or mu_pt[0] < 25 or mu_pt[1] < 20) continue;

      double pjet_scalpt_temp = 0;
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
      int countjets = 0;
      for( int i=0; i < pfj_size; i++){

         if( pfj_pt[i] > JETPT_LOW ){
            // clustered jets

            evtemp.jet_pt.push_back( pfj_pt[i] );
            evtemp.jet_phi.push_back( pfj_phi[i] );
            evtemp.jet_eta.push_back( pfj_eta[i] );

            if( pfj_pt[i]*pfj_l1l2l3[i] > 10 ){
               evtemp.jet_ptL123.push_back( pfj_pt[i]*pfj_l1l2l3[i] );
               evtemp.jet_ptT1.push_back( pfj_pt[i]*(pfj_l1l2l3[i] + 1 - pfj_l1[i]) );
            }else{
               evtemp.jet_ptL123.push_back( pfj_pt[i] );
               evtemp.jet_ptT1.push_back( pfj_pt[i] );
            }

         } else {
            // pseudojet with unclustered energy

            pjet_scalpt_temp += pfj_pt[i];
            pjet_px_temp += pfj_pt[i]*cos(pfj_phi[i]);
            pjet_py_temp += pfj_pt[i]*sin(pfj_phi[i]);

         }
         if(pfj_pt[i]*pfj_l1l2l3[i] > JETPT_HIGH and pfj_eta[i] > 2.3) countjets++;

      } // pfj loop
      //if(countjets > 1) std::cout << ev << ": " << countjets << std::endl;

      // fill pseudojet quantities
      evtemp.pjet_scalpt = pjet_scalpt_temp;
      evtemp.pjet_vectpt = sqrt( pjet_px_temp*pjet_px_temp + pjet_py_temp*pjet_py_temp );
      evtemp.pjet_phi = atan2( pjet_py_temp, pjet_px_temp );

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

   if(isMC) MatchMCjets();

   delete file;
   std::cout << "<--" << std::endl;
}

void Fitter::RunMinimizer(std::vector<event>& eventref_temp){
   std::cout << "RunMinimizer -->" << std::endl;

   gMinuit = new ROOT::Minuit2::Minuit2Minimizer ( ROOT::Minuit2::kMigrad );
   gMinuit->SetTolerance(0.001);
   gMinuit->SetStrategy(0);
   gMinuit->SetPrintLevel(3);

   fFunc = new ROOT::Math::Functor ( this, &Fitter::Min2LL, 12);
   gMinuit->SetFunction( *fFunc );
   gMinuit->SetVariable(0, "a1", 1.5, 0.01);
   gMinuit->SetVariable(1, "a2", 1.5, 0.01);
   gMinuit->SetVariable(2, "a3", 1.5, 0.01);
   gMinuit->SetVariable(3, "a4", 1.5, 0.01);
   gMinuit->SetVariable(4, "a5", 1.5, 0.01);
   gMinuit->SetVariable(5, "k0", 1.0, 0.01);
   gMinuit->SetVariable(6, "k1", 1.0, 0.01);
   gMinuit->SetVariable(7, "k2", 1.0, 0.01);
   gMinuit->SetVariable(8, "N1", 4.0, 0.01);
   gMinuit->SetVariable(9, "S1", 0.5, 0.01);
   gMinuit->SetVariable(10,"N2", 4.0, 0.01);
   gMinuit->SetVariable(11,"S2", 0.5, 0.01);

   // set event vector and minimize
   std::cout << "  minimize, first pass" << std::endl;
   eventvecPnt = &eventref_temp;
   gMinuit->Minimize();

   // new event vector with core of significance
   std::cout << "  build core sig vector" << std::endl;
   std::vector<event> eventvec_coretemp;
   for( std::vector<event>::iterator ev = eventvecPnt->begin(); ev < eventvecPnt->end(); ev++){

      if( ev->sig < 9 ){
         eventvec_coretemp.push_back( *ev );
      }

   }

   // minimize core sig
   std::cout << "  minimize, core sig" << std::endl;
   eventvecPnt = &eventvec_coretemp;
   gMinuit->SetStrategy(1);
   gMinuit->Minimize();
   gMinuit->Hesse();

   // load best-fit significance values
   std::cout << "  fill event vec with best-fit significance" << std::endl;
   const double *xmin = gMinuit->X();
   FindSignificance(xmin, eventref_temp);

   std::cout << "<--" << std::endl;
}

double Fitter::Min2LL(const double *x){

   // load significance values into eventref
   FindSignificance(x, *eventvecPnt);

   // event loop
   double m2ll = 0;
   for( std::vector<event>::iterator ev = eventvecPnt->begin(); ev < eventvecPnt->end(); ev++){
      int evIndex = int(ev - eventvecPnt->begin());
      m2ll += ev->sig + log(ev->det);
      if(countmin == -1){
         std::cout << std::setprecision(20)
            << evIndex << " + " << ev->sig + log(ev->det) << " ---> " << m2ll << std::endl;
         std::cout << std::setprecision(6);
      }
   }

   countmin++;
   return m2ll;
}

void Fitter::FindSignificance(const double *x, std::vector<event>& eventref_temp){

   // event loop
   for( std::vector<event>::iterator ev = eventref_temp.begin(); ev < eventref_temp.end(); ev++){

      int evIndex = int( ev - eventref_temp.begin() );
      bool fcout = false;//(countmin == 0 and evIndex < 10);

      double met_x=0;
      double met_y=0;
      double cov_xx=0;
      double cov_xy=0;
      double cov_yy=0;

      // clustered jets
      int countjet = 0;
      for(int i=0; i < int(ev->jet_pt.size()); i++){

         float feta = fabs(ev->jet_eta[i]);
         double cos = std::cos(ev->jet_phi[i]);
         double sin = std::sin(ev->jet_phi[i]);

         met_x -= cos*(ev->jet_ptT1[i]);
         met_y -= sin*(ev->jet_ptT1[i]);

         double dpt=0;
         double dph=0;

         if(false and evIndex == 410){
            std::cout << " ### " << ev->jet_pt[i] << " " << ev->jet_ptL123[i] << std::endl;
         }

         if( ev->jet_ptL123[i] > JETPT_HIGH ){
            int index=-1;
            if(feta<0.5) index=0;
            else if(feta<1.1) index=1;
            else if(feta<1.7) index=2;
            else if(feta<2.3) index=3;
            else{
               index=4;
            }
            countjet++;

            // CMS 2010 Resolutions -- parameterized by L123 corrected pt
            dpt = x[index] * (ev->jet_ptL123[i]) * dpt_(ev->jet_ptL123[i], ev->jet_eta[i]);
            dph =            (ev->jet_ptL123[i]) * dph_(ev->jet_ptL123[i], ev->jet_eta[i]);
         }
         else if( ev->jet_ptL123[i] > 0 ){
            int index=-1;
            if(feta<2.4) index=0;
            else if(feta<3) index=1;
            else index=2;

            // parameterized by T1 corrected pt
            dpt = x[5+index]*sqrt(ev->jet_ptT1[i]);
            dph = 0;
         }else{
            dpt = 0;
            dph = 0;
         }

         double dtt = dpt*dpt;
         double dff = dph*dph;
         cov_xx += dtt*cos*cos + dff*sin*sin;
         cov_xy += (dtt-dff)*cos*sin;
         cov_yy += dff*cos*cos + dtt*sin*sin;

         if( fcout and i == 0 ){
            std::cout << evIndex << " ########## PT: " << ev->jet_pt[i] << " "
               << ev->jet_ptT1[i] << " " << ev->jet_ptT1[i]*cos << " " 
               << ev->jet_ptT1[i]*sin << std::endl;
            std::cout << evIndex << " ########## DPT: " << dpt << " " << dph << " " << std::endl;
            std::cout << evIndex << " ########## COV: " << cov_xx << " " << cov_xy << " "
               << cov_yy << std::endl;
         }

      }

      // muons -- assume zero resolutions
      met_x -= cos(ev->muon_phi[0])*(ev->muon_pt[0]);
      met_y -= sin(ev->muon_phi[0])*(ev->muon_pt[0]);
      met_x -= cos(ev->muon_phi[1])*(ev->muon_pt[1]);
      met_y -= sin(ev->muon_phi[1])*(ev->muon_pt[1]);

      // unclustered energy -- parameterize by scalar sum of ET
      double cos = std::cos(ev->pjet_phi);
      double sin = std::sin(ev->pjet_phi);

      met_x -= cos*(ev->pjet_vectpt);
      met_y -= sin*(ev->pjet_vectpt);

      double ctt = x[8]*x[8] + x[9]*x[9]*(ev->pjet_scalpt);
      double cff = x[10]*x[10] + x[11]*x[11]*(ev->pjet_scalpt);

      cov_xx += ctt*cos*cos + cff*sin*sin;
      cov_xy += (ctt-cff)*cos*sin;
      cov_yy += cff*cos*cos + ctt*sin*sin;  

      if( fcout ){
         std::cout << evIndex << " ########## MU: " << ev->muon_pt[0] << " " << ev->muon_pt[1] << " "
            << ev->muon_phi[0] << " " << ev->muon_phi[1] << std::endl;
         std::cout << evIndex << " ########## PJET: " << ev->pjet_scalpt << " " 
            << ev->pjet_vectpt << " " << ev->pjet_phi << std::endl;
      }

      double det = cov_xx*cov_yy - cov_xy*cov_xy;

      double ncov_xx = cov_yy / det;
      double ncov_xy = -cov_xy / det;
      double ncov_yy = cov_xx / det;

      double sig = met_x*met_x*ncov_xx + 2*met_x*met_y*ncov_xy + met_y*met_y*ncov_yy;

      if( fcout ){
         std::cout << evIndex << " ########## SIG: " << sig << " " << det << std::endl;
      }

      ev->sig = sig;
      ev->det = det;
      if( false and countjet > 2 ) std::cout << evIndex << ": " << countjet << " "
         << ev->sig + log(ev->det) << std::endl;
      
      //std::cout << evIndex << ": " << countjet << std::endl;
   }
}

void Fitter::MatchMCjets(){
   std::cout << "  match MC jets" << std::endl;

   for( std::vector<event>::iterator ev = eventvec_MC.begin(); ev < eventvec_MC.end(); ev++){

      // loop through reco jets
      for(int ireco=0; ireco < int(ev->jet_pt.size()); ireco++){

         int matchIndex = -1;
         double dR = 1000;

         // loop through genjets
         for(int igen=0; igen < int(ev->genjet_pt.size()); igen++){

            double dphi = TVector2::Phi_mpi_pi( ev->jet_phi[ireco] 
                  - ev->genjet_phi[igen] );
            double deta = ev->jet_eta[ireco] 
               - ev->genjet_eta[igen];
            double dRtemp = sqrt( deta*deta + dphi*dphi );

            if( dRtemp < dR ){
               dR = dRtemp;
               matchIndex = igen;
            }

         } // genjets

         ev->jet_matchIndex.push_back( matchIndex );
         ev->jet_matchdR.push_back( dR );

      } // jets

      if( ev->jet_pt.size() != ev->jet_matchIndex.size() ){
         std::cout << "VECTOR SIZE MISMATCH" << std::endl;
      }

   } // event loop

}

