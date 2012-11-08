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
#include "TRandom3.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <list>

#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"

#include "METSigFit.h"
using namespace std;

//
// constructor and destructor
//

Fitter::Fitter(){
   // MINUIT variables
   gMinuit = 0;
   fFunc = 0;

   // jet pt bounds
   jetfitLOW = 10.0;
   jetfitHIGH = 30.0;
   jetcorrMIN = 10.0;

}

Fitter::~Fitter(){
   if (gMinuit) delete gMinuit;
   if (fFunc) delete fFunc;
}

//
// member definitions
//

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

void Fitter::ReadNtuple(const char* filename, vector<event>& eventref_temp, const int maxevents,
      const bool isMC){
   cout << "---> ReadNtuple" << endl;

   int v_size;
   float met_et;

   int pfj_size=0;
   float pfj_pt[1000];
   float pfj_phi[1000];
   float pfj_eta[1000];

   int mu_size=0;
   float mu_pt[100];
   float mu_px[100];
   float mu_py[100];
   float mu_pz[100];
   float mu_e[100];
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
   tree->SetBranchAddress("mu_px", mu_px);
   tree->SetBranchAddress("mu_py", mu_py);
   tree->SetBranchAddress("mu_pz", mu_pz);
   tree->SetBranchAddress("mu_e", mu_e);
   tree->SetBranchAddress("mu_phi", mu_phi);

   tree->SetBranchAddress("pfj_l1", pfj_l1);
   tree->SetBranchAddress("pfj_l1l2l3", pfj_l1l2l3);

   if(isMC){
      tree->SetBranchAddress("genj_size", &genj_size);
      tree->SetBranchAddress("genj_pt", genj_pt);
      tree->SetBranchAddress("genj_phi", genj_phi);
      tree->SetBranchAddress("genj_eta", genj_eta);
   }

   cout << " -----> fill event vector" << endl;

   int countev=0;
   for( int ev=0; ev<tree->GetEntries() and countev < maxevents; ev++){

      tree->GetEntry(ev);
      if( ev % 100000 == 0 and ev > 0) cout << "    -----> getting entry " << ev << endl;

      // ####################### Z PEAK FILTER #######################
      TLorentzVector mu1temp( mu_px[0], mu_py[0], mu_pz[0], mu_e[0] );
      TLorentzVector mu2temp( mu_px[1], mu_py[1], mu_pz[1], mu_e[1] );
      if( (mu1temp+mu2temp).M() < 86 or (mu1temp+mu2temp).M() > 96 ) continue;
      // #############################################################

      countev++;

      int pjet_size_temp = 0;
      double pjet_scalpt_temp = 0;
      double pjet_px_temp = 0;
      double pjet_py_temp = 0;

      event evtemp;

      // vertices
      evtemp.nvertices = v_size;

      // muons
      for( int i=0; i < mu_size; i++){
         TLorentzVector ptemp( mu_px[i], mu_py[i], mu_pz[i], mu_e[i] );
         evtemp.muon_pt.push_back( mu_pt[i] );
         evtemp.muon_phi.push_back( mu_phi[i] );
         evtemp.muon_4vect.push_back( ptemp );
      }

      // jets
      for( int i=0; i < pfj_size; i++){

         if( pfj_pt[i] > jetfitLOW ){
            // clustered jets

            evtemp.jet_pt.push_back( pfj_pt[i] );
            evtemp.jet_phi.push_back( pfj_phi[i] );
            evtemp.jet_eta.push_back( pfj_eta[i] );

            if( pfj_pt[i]*pfj_l1l2l3[i] > jetcorrMIN ){
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
            pjet_size_temp++;

         }

      } // pfj loop

      // pseudojet
      evtemp.pjet_scalpt = pjet_scalpt_temp;
      evtemp.pjet_vectpt = sqrt( pjet_px_temp*pjet_px_temp + pjet_py_temp*pjet_py_temp );
      evtemp.pjet_phi = atan2( pjet_py_temp, pjet_px_temp );
      evtemp.pjet_size = pjet_size_temp;

      // genjets
      if(isMC){
         for( int i=0; i < genj_size; i++){
            evtemp.genjet_pt.push_back( genj_pt[i] );
            evtemp.genjet_phi.push_back( genj_phi[i] );
            evtemp.genjet_eta.push_back( genj_eta[i] );
         }
      }

      // fill event vector
      eventref_temp.push_back( evtemp );

   } // event loop

   delete file;
}

void Fitter::RunMinimizer(vector<event>& eventref_temp){
   cout << "---> RunMinimizer" << endl;

   gMinuit = new ROOT::Minuit2::Minuit2Minimizer ( ROOT::Minuit2::kMigrad );
   gMinuit->SetTolerance(0.001);
   gMinuit->SetStrategy(0);
   gMinuit->SetPrintLevel(2);

   fFunc = new ROOT::Math::Functor ( this, &Fitter::Min2LL, 12 );
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
   cout << " -----> minimize, first pass" << endl;
   eventvecPnt = &eventref_temp;
   gMinuit->Minimize();

   // new event vector with core of significance
   cout << " -----> build core sig vector" << endl;
   vector<event> eventvec_coretemp;
   for( vector<event>::iterator ev = eventvecPnt->begin(); ev < eventvecPnt->end(); ev++){

      if( ev->sig < 9 ){
         eventvec_coretemp.push_back( *ev );
      }

   }

   // minimize core sig
   cout << " -----> minimize, core sig" << endl;
   eventvecPnt = &eventvec_coretemp;
   gMinuit->SetStrategy(1);
   gMinuit->Minimize();
   gMinuit->Hesse();

   // load best-fit significance values
   cout << " -----> fill event vec with best-fit significance" << endl;
   const double *xmin = gMinuit->X();
   FindSignificance(xmin, eventref_temp);

}

double Fitter::Min2LL(const double *x){

   // load significance values into eventvec
   FindSignificance(x, *eventvecPnt);

   // event loop
   double m2ll = 0;
   for( vector<event>::iterator ev = eventvecPnt->begin(); ev < eventvecPnt->end(); ev++){
      m2ll += ev->sig + log(ev->det);
   }

   return m2ll;
}

void Fitter::FindSignificance(const double *x, vector<event>& eventref_temp){

   // event loop
   for( vector<event>::iterator ev = eventref_temp.begin(); ev < eventref_temp.end(); ev++){

      double met_x=0;
      double met_y=0;
      double cov_xx=0;
      double cov_xy=0;
      double cov_yy=0;

      // clustered jets
      for(int i=0; i < int(ev->jet_pt.size()); i++){

         float feta = fabs(ev->jet_eta[i]);
         double c = cos(ev->jet_phi[i]);
         double s = sin(ev->jet_phi[i]);

         met_x -= c*(ev->jet_ptT1[i]);
         met_y -= s*(ev->jet_ptT1[i]);

         double dpt=0;
         double dph=0;

         // resolutions for two jet categories
         if( ev->jet_pt[i] > jetfitHIGH ){

            int index=-1;
            if(feta<0.5) index=0;
            else if(feta<1.1) index=1;
            else if(feta<1.7) index=2;
            else if(feta<2.3) index=3;
            else{
               index=4;
            }

            // CMS 2010 Resolutions -- parameterized by L123 corrected pt
            dpt = x[index] * (ev->jet_ptL123[i]) * dpt_(ev->jet_ptL123[i], ev->jet_eta[i]);
            dph =            (ev->jet_ptL123[i]) * dph_(ev->jet_ptL123[i], ev->jet_eta[i]);

         }
         else if( ev->jet_pt[i] > jetfitLOW ){

            int index=-1;
            if(feta<2.4) index=0;
            else if(feta<3) index=1;
            else index=2;

            // parameterized by T1 corrected pt
            dpt = x[5+index]*sqrt(ev->jet_ptT1[i]);
            dph = 0;

         }else{
            cout << "ERROR: JET PT OUT OF RANGE" << endl;
         }

         double dtt = dpt*dpt;
         double dff = dph*dph;
         cov_xx += dtt*c*c + dff*s*s;
         cov_xy += (dtt-dff)*c*s;
         cov_yy += dff*c*c + dtt*s*s;

      }

      // muons -- assume zero resolutions
      met_x -= cos(ev->muon_phi[0])*(ev->muon_pt[0]);
      met_y -= sin(ev->muon_phi[0])*(ev->muon_pt[0]);
      met_x -= cos(ev->muon_phi[1])*(ev->muon_pt[1]);
      met_y -= sin(ev->muon_phi[1])*(ev->muon_pt[1]);

      // unclustered energy -- parameterize by scalar sum of ET
      double c = cos(ev->pjet_phi);
      double s = sin(ev->pjet_phi);

      met_x -= c*(ev->pjet_vectpt);
      met_y -= s*(ev->pjet_vectpt);

      double ctt = x[8]*x[8] + x[9]*x[9]*(ev->pjet_scalpt);
      double cff = x[10]*x[10] + x[11]*x[11]*(ev->pjet_scalpt);

      cov_xx += ctt*c*c + cff*s*s;
      cov_xy += (ctt-cff)*c*s;
      cov_yy += cff*c*c + ctt*s*s;  

      double det = cov_xx*cov_yy - cov_xy*cov_xy;

      double ncov_xx = cov_yy / det;
      double ncov_xy = -cov_xy / det;
      double ncov_yy = cov_xx / det;

      double sig = met_x*met_x*ncov_xx + 2*met_x*met_y*ncov_xy + met_y*met_y*ncov_yy;

      ev->sig = sig;
      ev->det = det;

      ev->cov_xx = cov_xx;
      ev->cov_xy = cov_xy;
      ev->cov_yy = cov_yy;

      // fill qt, ut
      double qt_x = ev->muon_pt[0]*cos(ev->muon_phi[0])
         + ev->muon_pt[1]*cos(ev->muon_phi[1]);
      double qt_y = ev->muon_pt[0]*sin(ev->muon_phi[0])
         + ev->muon_pt[1]*sin(ev->muon_phi[1]);
      double qt = sqrt( qt_x*qt_x + qt_y*qt_y );

      double ut_x = -met_x - qt_x;
      double ut_y = -met_y - qt_y;
      double ut = sqrt( ut_x*ut_x + ut_y*ut_y );

      double ut_par = (ut_x*qt_x + ut_y*qt_y)/qt;
      double ut_perp = (ut_y*qt_x - qt_y*ut_x)/qt;

      ev->qt = qt;
      ev->ut = ut;
      ev->ut_par = ut_par;
      ev->ut_perp = ut_perp;

   }
}

void Fitter::MatchMCjets(vector<event>& eventref_temp){
   cout << "---> MatchMCjets" << endl;

   for( vector<event>::iterator ev = eventref_temp.begin(); ev < eventref_temp.end(); ev++){

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
         cout << "ERROR: VECTOR SIZE MISMATCH" << endl;
      }

   } // event loop

}

void Fitter::PlotsDataMC(vector<event>& eventref_MC, vector<event>& eventref_data, 
      const char* filename){

   // pileup reweighting
   int verts_MC [50] = {0};
   int verts_Data [50] = {0};
   double weights_MC [50];
   double weights_Data [50];
   fill_n(weights_MC,50,1.0);
   fill_n(weights_Data,50,1.0);
   for( vector<event>::iterator ev = eventref_MC.begin(); ev < eventref_MC.end(); ev++){
      verts_MC[ ev->nvertices ] += 1;
   }
   for( vector<event>::iterator ev = eventref_data.begin(); ev < eventref_data.end(); ev++){
      verts_Data[ ev->nvertices ] += 1;
   }
   for(int i=0; i < 50; i++){
      if( verts_MC[i] != 0){
         weights_MC[i] = (1.0*eventref_MC.size()/eventref_data.size())
            * (1.0*verts_Data[i]/verts_MC[i]);
      }
   }

   // histograms
   map<string, TH1*> histsData_;
   map<string, TH1*> histsMC_;
   map<string, TProfile*> profsData_;
   map<string, TProfile*> profsMC_;

   // data hists
   histsData_["muon_pt"] = new TH1D("muon_pt_Data", "Muon p_{T}", 100, 0, 200);
   histsData_["muon_invmass"] = new TH1D("muon_invmass_Data", "M_{#mu#mu}", 100, 0, 200);
   histsData_["njets"  ] = new TH1D("njets_Data", "N jets", 100, 0, 100);
   histsData_["jet_pt" ] = new TH1D("jet_pt_Data", "Jet p_{T}", 100, 0, 200);
   histsData_["jet1_pt"] = new TH1D("jet1_pt_Data", "Jet 1 p_{T}", 100, 0, 200);
   histsData_["pjet_size"  ] = new TH1D("pjet_size_Data", "N jets", 100, 0, 500);
   histsData_["pjet_scalpt"] = new TH1D("pjet_scalpt_Data", "Pseudojet Scalar p_{T}", 100, 0, 500);
   histsData_["pjet_vectpt"] = new TH1D("pjet_vectpt_Data", "Pseudojet Scalar p_{T}", 100, 0, 200);
   histsData_["qt"] = new TH1D("qt_Data", "q_{T}", 100, 0, 200);
   histsData_["ut_par"] = new TH1D("ut_par_Data", "|u_{T}|_{#parallel}", 100, 0, 100);
   histsData_["nvert"] = new TH1D("nvert_Data", "N Vertices", 100, 0, 100);
   histsData_["cov_xx"] = new TH1D("cov_xx_Data", "Cov_{xx}", 100, 0, 500);
   histsData_["cov_xy"] = new TH1D("cov_xy_Data", "Cov_{xy}", 100, -150, 150);
   histsData_["cov_yy"] = new TH1D("cov_yy_Data", "Cov_{yy}", 100, 0, 500);
   histsData_["sig"] = new TH1D("sig_Data", "Significance", 100, 0, 50);
   histsData_["det"] = new TH1D("det_Data", "Determinant", 100, 0, 100000);
   histsData_["pchi2"] = new TH1D("pchi2_Data", "P(#chi^{2})", 100, 0, 1);

   // profile histograms
   profsData_["psig_vert"] = new TProfile("psig_vert_Data",
         "Significance vs. N Vertices;N Vertices;<S_{E}>", 15, 0, 15);
   profsData_["psig_qt"] = new TProfile("psig_qt_Data",
         "Significance vs. q_{T};q_{T} (GeV);<S_{E}>", 15, 0, 100);
   profsData_["presp_qt"] = new TProfile("presp_qt_Data",
         "Response = |<u_{#parallel}>|/q_{T} vs. q_{T};q_{T} (GeV);Response", 25, 0, 100);

   // clone data hists for MC
   for(map<string,TH1*>::const_iterator it = histsData_.begin();
         it != histsData_.end(); it++){

      string hname = it->first;
      TH1D *hist = (TH1D*)it->second;
      histsMC_[hname] = (TH1D*)hist->Clone((char*)hname.c_str());

   }
   for(map<string,TProfile*>::const_iterator it = profsData_.begin();
         it != profsData_.end(); it++){

      string pname = it->first;
      TProfile *prof = (TProfile*)it->second;
      profsMC_[pname] = (TProfile*)prof->Clone((char*)pname.c_str());

   }


   // fill hists
   for( int i=0; i < 2; i++ ){

      map<string, TH1*> hists_;
      map<string, TProfile*> profs_;
      vector<event>::iterator iter_begin;
      vector<event>::iterator iter_end;
      double *weights;

      if( i==0 ){
         hists_ = histsMC_;
         profs_ = profsMC_;
         iter_begin = eventref_MC.begin();
         iter_end = eventref_MC.end();
         weights = weights_MC;
      }
      if( i==1 ){
         hists_ = histsData_;
         profs_ = profsData_;
         iter_begin = eventref_data.begin();
         iter_end = eventref_data.end();
         weights = weights_Data;
      }

      for( vector<event>::iterator ev = iter_begin; ev < iter_end; ev++ ){
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
         if( ev->jet_ptL123.size() > 0 ){
            hists_["jet1_pt"]->Fill( ev->jet_ptL123[0] , weights[nvert] );
         }

         // pseudojet
         hists_["pjet_scalpt"]->Fill( ev->pjet_scalpt , weights[nvert] );
         hists_["pjet_vectpt"]->Fill( ev->pjet_vectpt , weights[nvert] );
         hists_["pjet_size"]->Fill( ev->pjet_size , weights[nvert] );

         // other observables
         hists_["nvert"]->Fill( nvert , weights[nvert] );

         hists_["qt"]->Fill( ev->qt, weights[nvert] );
         hists_["ut_par"]->Fill( fabs(ev->ut_par), weights[nvert] );

         hists_["cov_xx"]->Fill( ev->cov_xx, weights[nvert] );
         hists_["cov_xy"]->Fill( ev->cov_xy, weights[nvert] );
         hists_["cov_yy"]->Fill( ev->cov_yy, weights[nvert] );
         hists_["sig"]->Fill( ev->sig, weights[nvert] );
         hists_["det"]->Fill( ev->det, weights[nvert] );
         hists_["pchi2"]->Fill( TMath::Prob(ev->sig,2), weights[nvert] );

         // profiles
         profs_["psig_vert"]->Fill( ev->nvertices, ev->sig );
         profs_["psig_qt"]->Fill( ev->qt, ev->sig );
         profs_["presp_qt"]->Fill( ev->qt, -(ev->ut_par)/(ev->qt) );
      }
   }

   // draw hists and write to file
   TFile *file = new TFile(filename,"RECREATE");
   file->cd();

   for(map<string,TH1*>::const_iterator it = histsData_.begin();
         it != histsData_.end(); it++){

      string hname = it->first;
      TH1D *histData = (TH1D*)it->second;
      TH1D *histMC = (TH1D*)histsMC_[hname];

      TCanvas *canvas  = new TCanvas( (char*)hname.c_str(), (char*)hname.c_str(), 700, 700 );
      canvas->cd();

      histMC->SetLineColor(2);
      histMC->Scale( double(eventref_data.size()) / eventref_MC.size() );
      histData->SetLineColor(1);
      histData->SetMarkerStyle(20);

      histMC->SetMaximum( 1.1*max(histMC->GetMaximum(), histData->GetMaximum()) );

      histMC->Draw();
      histData->Draw("EP same");

      canvas->Write();
   }
   for(map<string,TProfile*>::const_iterator it = profsData_.begin();
         it != profsData_.end(); it++){

      string pname = it->first;
      TProfile *profData = (TProfile*)it->second;
      TProfile *profMC = (TProfile*)profsMC_[pname];

      TCanvas *canvas  = new TCanvas( (char*)pname.c_str(), (char*)pname.c_str(), 700, 700 );
      canvas->cd();

      profMC->SetLineColor(2);
      profData->SetLineColor(1);
      profData->SetMarkerStyle(20);

      profMC->SetMaximum( 1.1*max(profMC->GetMaximum(), profData->GetMaximum()) );
      profMC->SetMinimum( 0.8*min(profMC->GetMinimum(), profData->GetMinimum()) );

      profMC->Draw("HIST");
      profData->Draw("EP same");

      canvas->Write();
   }

}
