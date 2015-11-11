#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TFileCollection.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TRandom3.h"
#include "TProfile.h"
#include "TString.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "TVirtualFFT.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <list>

#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"

#include "Math/GSLRndmEngines.h"
#include "Math/Random.h"

#include "METSigFit.h"
using namespace std;

//
// constructor and destructor
//

Fitter::Fitter(double b /*=1*/){
   // MINUIT variables
   gMinuit = 0;
   fFunc = 0;

   // jet pt bounds
   jetbinpt = 20.0;
   jetcorrpt = 10.0;

   // significance cut for minimization
   significance_cut = false;

   // MET type to use in FindSignificance
   met_type = 4;

   // compute weighted errors
   TH1::SetDefaultSumw2();

   // data hists for plotting
   rebin = b;
   histsData_["lepton_pt"] = new TH1D("lepton_pt_Data", "Muon p_{T}", 50/b, 0, 200);
   histsData_["lepton_eta"] = new TH1D("lepton_eta_Data", "Muon #eta", 50/b, -5, 5);
   histsData_["muon_invmass"] = new TH1D("muon_invmass_Data", "M_{#mu#mu}", 50/b, 0, 200);
   histsData_["njets"  ] = new TH1D("njets_Data", "N jets", 10/b, 0, 10);
   histsData_["jet_pt" ] = new TH1D("jet_pt_Data", "Jet p_{T}", 50/b, 0, 200);
   histsData_["jet1_pt"] = new TH1D("jet1_pt_Data", "Jet 1 p_{T}", 50/b, 0, 500);
   histsData_["jet_eta" ] = new TH1D("jet_eta_Data", "Jet #eta, p_{T} > 30 GeV", 50/b, -5, 5);
   histsData_["pjet_size"  ] = new TH1D("pjet_size_Data", "N jets", 50/b, 0, 500);
   histsData_["pjet_scalpt"] = new TH1D("pjet_scalpt_Data", "Pseudojet Scalar p_{T} ", 200/b, 0, 2000);
   histsData_["pjet_scalptL123"] = new TH1D("pjet_scalptL123_Data", "Pseudojet Scalar p_{T} (L123-corrected)", 200/b, 0, 2000);
   histsData_["pjet_vectptL123"] = new TH1D("pjet_vectptL123_Data", "Pseudojet Scalar p_{T} (L123-corrected)", 50/b, 0, 200);
   histsData_["pjet_scalptT1"] = new TH1D("pjet_scalptT1_Data", "Pseudojet Scalar p_{T} (T1-corrected)", 200/b, 0, 2000);
   histsData_["pjet_vectptT1"] = new TH1D("pjet_vectptT1_Data", "Pseudojet Scalar p_{T} (T1-corrected)", 50/b, 0, 200);
   histsData_["pjet_phi"] = new TH1D("pjet_phi_Data", "Pseudojet Phi", 50/b, -3.5, 3.5);
   histsData_["qt"] = new TH1D("qt_Data", "q_{T}", 50/b, 0, 500);
   histsData_["qx"] = new TH1D("qx_Data", "q_{x}", 50/b, -50, 50);
   histsData_["ut_par"] = new TH1D("ut_par_Data", "|u_{T}|_{#parallel}", 50/b, 0, 100);
   histsData_["nvert"] = new TH1D("nvert_Data", "N Vertices", 50/b, 0, 50);
   histsData_["cov_xx"] = new TH1D("cov_xx_Data", "Cov_{xx}", 50/b, 0, 500);
   histsData_["cov_xy"] = new TH1D("cov_xy_Data", "Cov_{xy}", 50/b, -150, 150);
   histsData_["cov_yy"] = new TH1D("cov_yy_Data", "Cov_{yy}", 50/b, 0, 500);
   histsData_["met"] = new TH1D("met_Data", "Missing E_{T}", 50/b, 0, 100);
   histsData_["met_200"] = new TH1D("met_200_Data", "Missing E_{T}", 50/b, 0, 200);
   histsData_["met_300"] = new TH1D("met_300_Data", "Missing E_{T}", 30/b, 0, 300);
   histsData_["sig"] = new TH1D("sig_Data", "Significance", 50/b, 0, 500);
   histsData_["sig_100"] = new TH1D("sig_100_Data", "Significance", 50/b, 0, 100);
   histsData_["sig_15"] = new TH1D("sig_15_Data", "Significance", 50/b, 0, 15);
   histsData_["sig_old"] = new TH1D("sig_old_Data", "Old Significance", 50/b, 0, 500);
   histsData_["det"] = new TH1D("det_Data", "Determinant", 50/b, 0, 100000);
   histsData_["pchi2"] = new TH1D("pchi2_Data", "P(#chi^{2})", 50/b, 0, 1);
   histsData_["pchi2_old"] = new TH1D("pchi2_old-Data", "P(#chi^{2}) from Old Sig", 50/b, 0, 1);
   histsData_["logpchi2"] = new TH1D("logpchi2_Data", "log(P(#chi^{2}))", 50/b, -50, 0);
   histsData_["cov_xx_highpt"] = new TH1D("cov_xx_highpt_Data", "Cov_{xx} High-p_{T} Jets", 50/b, 0, 500);
   histsData_["cov_xx_pjet"] = new TH1D("cov_xx_pjet_Data", "Cov_{xx} Pseudojet", 50/b, 0, 500);
   histsData_["cov_xx_ratio"] = new TH1D("cov_xx_ratio_Data", "Cov_{xx} High-p_{T}/Total", 50/b, 0, 1 );
   histsData_["met_varx"] = new TH1D("met_varx_Data","#sigma_{x}^{2} for MET smearing",50/b,-100,100);
   histsData_["met_vary"] = new TH1D("met_vary_Data","#sigma_{y}^{2} for MET smearing",50/b,-100,100);
   histsData_["met_rho"] = new TH1D("met_rho_Data", "#rho for MET smearing", 50/b, -10, 10);
   histsData_["pjet_axesratio"] = new TH1D("pjet_axesratio_Data",
         "Pseudo-jet Ratio of Major/Semi-Major Axes", 50/b, 0, 1);
   histsData_["pjet_tiltangle"] = new TH1D("pjet_tiltangle_Data", "Tilt Angle of Pseudo-jet",
         50/b, -2, 2);
   histsData_["pjet_tiltangle_rel"] = new TH1D("pjet_tiltangle_rel_Data",
      "Tilt Angle of Pseudo-jet relative to Pseudo-jet Momentum", 50/b, -2, 2);

   // profile histograms
   profsData_["psig_nvert"] = new TH2D("psig_nvert_Data",
         "Significance vs. N Vertices;N Vertices;<S_{E}>", 30/b, 0, 30, 50/b, 0, 50);
   profsData_["psig_qt"] = new TH2D("psig_qt_Data",
         "Significance vs. q_{T};q_{T} (GeV);<S_{E}>", 15/b, 0, 50, 100/b, 0, 50);
   profsData_["presp_qt"] = new TH2D("presp_qt_Data",
         "Response = |<u_{#parallel}>|/q_{T} vs. q_{T};q_{T} (GeV);Response", 25/b, 0, 50, 100/b, -100, 100);
   profsData_["pMET_nvert"] = new TH2D("pMET_nvert_Data",
         "MET vs. N Vertices;N Vertices;<MET>", 30/b, 0, 30, 50/b, 0, 100);
   profsData_["pjet_scalptL123_nvert"] = new TH2D("pjet_scalpt_nvert_Data",
         "Pseudojet Scalar p_{T} vs. N Vertices", 30/b, 0, 30, 200/b, 0, 2000);
   profsData_["jet_pt_nvert"] = new TH2D("jet_pt_nvert_Data",
         "Jet p_{T} vs. N Vertices", 30/b, 0, 30, 50/b, 0, 200);
   profsData_["njets_nvert"] = new TH2D("njets_nvert_Data",
         "N jets vs. N Vertices", 30/b, 0, 30, 50/b, 0, 100);
   profsData_["sig_met"] = new TH2D("sig_met_Data",
         "Significance vs. MET", 500/b, 0, 100, 500/b, 0, 30);

   // clone data hists for MC
   for(map<string,TH1*>::const_iterator it = histsData_.begin();
         it != histsData_.end(); it++){

      string hname = it->first;
      TH1D *hist = (TH1D*)it->second;
      histsMC_[hname] = (TH1D*)hist->Clone( (hname+"_MC").c_str() );

      histsMC_signal_[hname] = (TH1D*)hist->Clone( (hname+"_MC_signal").c_str() );
      histsMC_top_[hname] = (TH1D*)hist->Clone( (hname+"_MC_top").c_str() );
      histsMC_top_dileptonic_[hname] = (TH1D*)hist->Clone( (hname+"_MC_top_dileptonic").c_str() );
      histsMC_top_hadronic_[hname] = (TH1D*)hist->Clone( (hname+"_MC_top_hadronic").c_str() );
      histsMC_top_single_[hname] = (TH1D*)hist->Clone( (hname+"_MC_top_single").c_str() );
      histsMC_EWK_[hname] = (TH1D*)hist->Clone( (hname+"_MC_EWK").c_str() );
      histsMC_QCD_[hname] = (TH1D*)hist->Clone( (hname+"_MC_QCD").c_str() );
      histsMC_gamma_[hname] = (TH1D*)hist->Clone( (hname+"_MC_gamma").c_str() );
      histsMC_DY_[hname] = (TH1D*)hist->Clone( (hname+"_MC_DY").c_str() );

   }
   for(map<string,TH2*>::const_iterator it = profsData_.begin();
         it != profsData_.end(); it++){

      string pname = it->first;
      TH2 *prof = (TH2*)it->second;
      profsMC_[pname] = (TH2*)prof->Clone((char*)pname.c_str());

   }
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

void Fitter::ReadNtuple(string path, vector<string>& filenames, 
      vector<event>& eventref_temp, const double fracevents,
      const bool isMC, string process, const bool do_resp_correction, 
      const int start_evt_num, const int end_evt_num , const double jvar){
   cout << "---> ReadNtuple " << process << endl;

   std::vector<double> *lep_pt=0, *lep_energy=0, *lep_phi=0, *lep_eta=0;
   std::vector<double> *jet_pt=0, *jet_energy=0, *jet_phi=0, *jet_eta=0;
   std::vector<double> *jet_sigmapt=0, *jet_sigmaphi=0;
   std::vector<bool> *jet_passid=0;
   double met_pt=0, met_energy=0, met_phi=0, met_eta=0, met_sumpt=0;
   int nvertices=0;
   double mcweight = 1.0;
   int run = 0;

   TChain *tree = new TChain("events");
   for( vector<string>::iterator file = filenames.begin(); file != filenames.end(); file++){
      TString fn = path + "/" + (*file);
      tree->Add( fn );
   }

   tree->SetBranchAddress("run", &run);
   tree->SetBranchAddress("mcweight", &mcweight);

   tree->SetBranchAddress("muon_pt", &lep_pt);
   tree->SetBranchAddress("muon_energy", &lep_energy);
   tree->SetBranchAddress("muon_phi", &lep_phi);
   tree->SetBranchAddress("muon_eta", &lep_eta);

   tree->SetBranchAddress("jet_pt", &jet_pt);
   tree->SetBranchAddress("jet_passid", &jet_passid);
   tree->SetBranchAddress("jet_energy", &jet_energy);
   tree->SetBranchAddress("jet_phi", &jet_phi);
   tree->SetBranchAddress("jet_eta", &jet_eta);
   tree->SetBranchAddress("jet_sigmapt", &jet_sigmapt);
   tree->SetBranchAddress("jet_sigmaphi", &jet_sigmaphi);

   tree->SetBranchAddress("met_pt", &met_pt);
   tree->SetBranchAddress("met_energy", &met_energy);
   tree->SetBranchAddress("met_phi", &met_phi);
   tree->SetBranchAddress("met_eta", &met_eta);
   tree->SetBranchAddress("met_sumpt", &met_sumpt);

   tree->SetBranchAddress("nvertices", &nvertices);

   cout << " -----> fill event vector" << endl;

   int numevents = (end_evt_num == -1) ? fracevents*tree->GetEntries(): end_evt_num-start_evt_num;
   int start_evt=0, end_evt=0;
   if( end_evt_num == -1 ){
      start_evt = 0;
      end_evt = numevents;
   }else{
      start_evt = start_evt_num;
      end_evt = end_evt_num;
   }
   eventref_temp.reserve( numevents );

   for( int ev=start_evt; ev<end_evt; ev++){

      tree->GetEntry(ev);
      if( ev % 100000 == 0 and ev > 0) cout << "    -----> getting entry " << ev << endl;

      int sgn_weight = mcweight < 0 ? -1.0 : 1.0;

      event evtemp;

      evtemp.process = process;
      evtemp.weight = 1.0*sgn_weight;

      // MC cross-sections (pb)
      // obtained at https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
      // and http://arxiv.org/pdf/1105.0020.pdf for WZ, ZZ.
      double xsec_dyjetstoll           = 2008.4;
      double xsec_ttjets               = 831.76;
      /*
      double xsec_ttjets_0lept         = 0.46*831.76;
      double xsec_ttjets_1lept         = 0.45*831.76;
      double xsec_ttjets_2lept         = 0.09*831.76;
      */
      double xsec_ww                   = 118.7;
      double xsec_wz                   = 46.74; 
      double xsec_zz                   = 15.99;
      double xsec_wjetstolnu           = 20508.90;

      // MC sample size (events), obtained from DAS.
      // miniAOODv2
      /*
      //int nevts_dyjetstoll          = 9042031; //madgraph
      int nevts_dyjetstoll          = 28747969; //amc@nlo
      //int nevts_dyjetstoll          = 28825132; //miniAODv1
      int nevts_ttjets              = 11344206;
      int nevts_wjetstolnu          = 24184766;
      int nevts_ww                  = 989608;
      int nevts_wz                  = 996920;
      int nevts_zz                  = 998848;
      */
      // miniAODv1
      int nevts_dyjetstoll          = 28825132;
      int nevts_ttjets              = 42730273;
      int nevts_wjetstolnu          = 24151270;
      int nevts_ww                  = 994416;
      int nevts_wz                  = 991232;
      int nevts_zz                  = 996168;

      if( isMC ){

         //evtemp.weight *= 616.9; // pb in dataset

         if( process.compare("DYJetsToLL") == 0 ) evtemp.weight *= xsec_dyjetstoll/nevts_dyjetstoll;
         else if( process.compare("TTJets") == 0 ) evtemp.weight *= xsec_ttjets/nevts_ttjets;
         else if( process.compare("WW") == 0 ) evtemp.weight *= xsec_ww/nevts_ww;
         else if( process.compare("WZ") == 0 ) evtemp.weight *= xsec_wz/nevts_wz;
         else if( process.compare("ZZ") == 0 ) evtemp.weight *= xsec_zz/nevts_zz;
         else if( process.compare("WJetsToLNu") == 0 )
            evtemp.weight *= xsec_wjetstolnu/nevts_wjetstolnu;
         else cout << "No Xsection for process " << process << endl;

      }
      // nvertices
      evtemp.nvertices = nvertices;

      // leptons
      evtemp.lepton_pt = *lep_pt;
      evtemp.lepton_phi = *lep_phi;
      evtemp.lepton_eta = *lep_eta;
      evtemp.lepton_energy = *lep_energy;

      // jets
      evtemp.jet_pt = *jet_pt;
      evtemp.jet_passid = *jet_passid;
      evtemp.jet_phi = *jet_phi;
      evtemp.jet_eta = *jet_eta;
      evtemp.jet_sigmapt = *jet_sigmapt;
      evtemp.jet_sigmaphi = *jet_sigmaphi;

      // met
      evtemp.met_pt = met_pt;
      evtemp.met_phi = met_phi;

      // pseudo-jet
      evtemp.pjet_scalpt = met_sumpt > 0 ? met_sumpt : 0.0;

      int njets = 0;
      for( int j=0; j < int(evtemp.jet_pt.size()); j++){
         if( fabs(evtemp.jet_eta[j]) < 2.5 and evtemp.jet_pt[j] > 30 ) njets++;
      }
      if( njets >=2 ){
         eventref_temp.push_back( evtemp );
      }

   } // event loop

   return;
}

void Fitter::RunMinimizer(vector<event>& eventref_temp){
   cout << "---> RunMinimizer" << endl;

   gMinuit = new ROOT::Minuit2::Minuit2Minimizer ( ROOT::Minuit2::kMigrad );
   gMinuit->SetTolerance(1000.0);
   gMinuit->SetStrategy(0);
   gMinuit->SetPrintLevel(2);

   fFunc = new ROOT::Math::Functor ( this, &Fitter::Min2LL, 7 );
   gMinuit->SetFunction( *fFunc );
   gMinuit->SetVariable(0, "a1", 1.0, 0.05);
   gMinuit->SetVariable(1, "a2", 1.0, 0.05);
   gMinuit->SetVariable(2, "a3", 1.0, 0.05);
   gMinuit->SetVariable(3, "a4", 1.0, 0.05);
   gMinuit->SetVariable(4, "a5", 1.0, 0.05);
   gMinuit->SetVariable(5, "N1", 0.0, 0.05);
   gMinuit->SetVariable(6, "S1", 0.5, 0.05);

   // set event vector and minimize
   cout << " -----> minimize, first pass" << endl;
   eventvecPnt = &eventref_temp;
   gMinuit->Minimize();

   // significance for cut
   for( vector<event>::iterator ev = eventvecPnt->begin(); ev < eventvecPnt->end(); ev++){
      ev->sig_init = ev->sig;
   }

   // minimize core sig
   cout << " -----> minimize, core sig" << endl;
   significance_cut = true;
   gMinuit->SetStrategy(1);
   gMinuit->Minimize();
   gMinuit->Hesse();

   significance_cut = false;

   // load best-fit significance values
   cout << " -----> fill event vec with best-fit significance" << endl;
   xmin = gMinuit->X();
   FindSignificance(xmin, eventref_temp);

}

double Fitter::Min2LL(const double *x){

   // load significance values into eventvec
   FindSignificance(x, *eventvecPnt);

   // event loop
   double m2ll = 0;
   for( vector<event>::iterator ev = eventvecPnt->begin(); ev < eventvecPnt->end(); ev++){
      if( !(significance_cut and ev->sig_init > 9) ){
         m2ll += ev->weight*(ev->sig + log(ev->det));
      }
   }

   return m2ll;
}

void Fitter::FindSignificance(const double *x, vector<event>& eventref_temp){

   for( vector<event>::iterator ev = eventref_temp.begin(); ev < eventref_temp.end(); ev++){

      // metsig covariance
      double cov_xx = 0;
      double cov_xy = 0;
      double cov_yy = 0;

      // jets
      for(unsigned int i=0; i < ev->jet_pt.size(); i++){
         double jpt = ev->jet_pt[i];
         double jeta = ev->jet_eta[i];
         double feta = fabs(jeta);
         double c = cos(ev->jet_phi[i]);
         double s = sin(ev->jet_phi[i]);

         int index=-1;
         if(feta<0.5) index=0;
         else if(feta<1.1) index=1;
         else if(feta<1.7) index=2;
         else if(feta<2.3) index=3;
         else{
            index=4;
         }

         double dpt = x[index]*jpt*ev->jet_sigmapt[i];
         double dph =          jpt*ev->jet_sigmaphi[i];

         double dtt = dpt*dpt;
         double dff = dph*dph;
         cov_xx += dtt*c*c + dff*s*s;
         cov_xy += (dtt-dff)*c*s;
         cov_yy += dff*c*c + dtt*s*s;
      }

      // pseudo-jet
      double ctt = x[5]*x[5] + x[6]*x[6]*(ev->pjet_scalpt);
      cov_xx += ctt;
      cov_yy += ctt;

      // compute significance
      double met_x = ev->met_pt * cos(ev->met_phi);
      double met_y = ev->met_pt * sin(ev->met_phi);

      double det = cov_xx*cov_yy - cov_xy*cov_xy;
      double ncov_xx = cov_yy / det;
      double ncov_xy = -cov_xy / det;
      double ncov_yy = cov_xx / det;

      double sig = met_x*met_x*ncov_xx + 2*met_x*met_y*ncov_xy + met_y*met_y*ncov_yy;

      // load into eventvec
      ev->sig = sig;
      ev->det = det;

      ev->cov_xx = cov_xx;
      ev->cov_xy = cov_xy;
      ev->cov_yy = cov_yy;

      //if( ev->pjet_scalpt <= 0 ) cout << ev->pjet_scalpt << endl;

      // fill qt, ut
      double qt_x=0, qt_y=0;
      for(int i=0; i < int(ev->lepton_pt.size()); i++){
         qt_x += ev->lepton_pt[i]*cos(ev->lepton_phi[i]);
         qt_y += ev->lepton_pt[i]*sin(ev->lepton_phi[i]);
      }
      double qt = sqrt( qt_x*qt_x + qt_y*qt_y );

      double ut_x = -met_x - qt_x;
      double ut_y = -met_y - qt_y;
      double ut = sqrt( ut_x*ut_x + ut_y*ut_y );

      double ut_par = (ut_x*qt_x + ut_y*qt_y)/qt;
      double ut_perp = (ut_y*qt_x - qt_y*ut_x)/qt;

      TLorentzVector mu1, mu2;
      mu1.SetPtEtaPhiE( ev->lepton_pt[0], ev->lepton_eta[0], ev->lepton_phi[0], ev->lepton_energy[0] );
      mu2.SetPtEtaPhiE( ev->lepton_pt[1], ev->lepton_eta[1], ev->lepton_phi[1], ev->lepton_energy[1] );
      ev->invmass = (mu1+mu2).M();

      ev->qt = qt;
      ev->qx = qt_x;
      ev->ut = ut;
      ev->ut_par = ut_par;
      ev->ut_perp = ut_perp;

   }

   return;
}

//
// Yimin's full jet resolution shapes
//
void Fitter::ComplexMult(int entries, double *re1, double *im1, double *re2, double *im2, double *reF, double *imF){
   for (int i=0; i<entries; i++){
      double tempre1=re1[i];
      double tempim1=im1[i];
      double tempre2=re2[i];
      double tempim2=im2[i];
      reF[i]=tempre1*tempre2-tempim1*tempim2;
      imF[i]=tempre1*tempim2+tempre2*tempim1;
   }
}

void Fitter::FillHists(vector<event>& eventref, string stackmode){
   if( eventref.size() == 0 ) return;

   vector<event>::iterator iter_begin = eventref.begin();
   vector<event>::iterator iter_end = eventref.end();

   map<string, TH1*> hists_;
   map<string, TH2*> profs_;

   if( iter_begin->process.compare("Data") == 0 ){
      hists_ = histsData_;
      profs_ = profsData_;
   } else {
      profs_ = profsMC_;
      if( stackmode.compare("Zmumu") == 0 ){

         if( iter_begin->process.compare("DYJetsToLL") == 0 ) hists_ = histsMC_signal_;
         else if( iter_begin->process.compare("TTJets") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->process.compare("Tbar_tW-channel") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->process.compare("T_tW-channel") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->process.compare("WJetsToLNu") == 0 ) hists_ = histsMC_EWK_;
         else if( iter_begin->process.compare("WW") == 0 ) hists_ = histsMC_EWK_;
         else if( iter_begin->process.compare("WZ") == 0 ) hists_ = histsMC_EWK_;
         else if( iter_begin->process.compare("ZZ") == 0 ) hists_ = histsMC_EWK_;
         else cout << "Histogram fill error, process " << iter_begin->process << endl;

      }
      if( stackmode.compare("Wenu") == 0 or stackmode.compare("Wenu_loose") == 0 ){

         if( iter_begin->process.compare("WJetsToLNu") == 0 ) hists_ = histsMC_signal_;
         else if( iter_begin->process.compare("DYJetsToLL") == 0 ) hists_ = histsMC_DY_;
         else if( iter_begin->process.compare("DYJetsToLL_M-50") == 0 ) hists_ = histsMC_DY_;
         else if( iter_begin->process.compare("DYJetsToLL_M-10To50") == 0 ) hists_ = histsMC_DY_;
         else if( iter_begin->process.compare("TTJets") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->process.compare("Tbar_tW-channel") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->process.compare("T_tW-channel") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->process.compare("QCD_EMEnriched_20_30") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_EMEnriched_30_80") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_EMEnriched_80_170") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_EMEnriched_170_250") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_EMEnriched_250_350") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_EMEnriched_350") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_BCtoE_20_30") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_BCtoE_30_80") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_BCtoE_80_170") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_BCtoE_170_250") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_BCtoE_250_350") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("Gamma_0_15") == 0 ) hists_ = histsMC_gamma_;
         else if( iter_begin->process.compare("Gamma_15_30") == 0 ) hists_ = histsMC_gamma_;
         else if( iter_begin->process.compare("Gamma_30_50") == 0 ) hists_ = histsMC_gamma_;
         else if( iter_begin->process.compare("Gamma_50_80") == 0 )  hists_ = histsMC_gamma_;
         else if( iter_begin->process.compare("Gamma_80_120") == 0 ) hists_ = histsMC_gamma_;
         else if( iter_begin->process.compare("Gamma_120_170") == 0 ) hists_ = histsMC_gamma_;  
         else if( iter_begin->process.compare("Gamma_170_300") == 0 ) hists_ = histsMC_gamma_;  
         else if( iter_begin->process.compare("Gamma_300_470") == 0 ) hists_ = histsMC_gamma_;  
         else if( iter_begin->process.compare("WW") == 0 ) hists_ = histsMC_EWK_;
         else if( iter_begin->process.compare("WZ") == 0 ) hists_ = histsMC_EWK_;
         else if( iter_begin->process.compare("ZZ") == 0 ) hists_ = histsMC_EWK_;
         else cout << "Histogram fill error, process " << iter_begin->process << endl;

      }
      if( stackmode.compare("Dijet") == 0 ){

         if( iter_begin->process.find("QCD") != string::npos ) hists_ = histsMC_signal_;
         else if( iter_begin->process.compare("DYJetsToLL") == 0 ) hists_ = histsMC_DY_;
         else if( iter_begin->process.compare("TTJets") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->process.compare("Tbar_tW-channel") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->process.compare("T_tW-channel") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->process.compare("WJetsToLNu") == 0 ) hists_ = histsMC_EWK_;
         else if( iter_begin->process.compare("WW") == 0 ) hists_ = histsMC_EWK_;
         else if( iter_begin->process.compare("WZ") == 0 ) hists_ = histsMC_EWK_;
         else if( iter_begin->process.compare("ZZ") == 0 ) hists_ = histsMC_EWK_;
         else cout << "Histogram fill error, process " << iter_begin->process << endl;

      }
      if( stackmode.compare("Ttbar0lept") == 0){

         if( iter_begin->process.compare("TTJets_Hadronic") == 0 ) hists_ = histsMC_signal_;
         else if( iter_begin->process.compare("TTJets_FullLept") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->process.compare("TTJets_SemiLept") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->process.compare("DYJetsToLL") == 0 ) hists_ = histsMC_DY_;
         else if( iter_begin->process.compare("Tbar_tW-channel") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->process.compare("T_tW-channel") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->process.compare("WJetsToLNu") == 0 ) hists_ = histsMC_EWK_;
         else if( iter_begin->process.compare("QCD_15_30") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_30_50") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_50_80") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_80_120") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_120_170") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_170_300") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_300_470") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_470_600") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_600_800") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_800_1000") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_1000_1400") == 0 ) hists_ = histsMC_QCD_;
         else cout << "Histogram fill error, process " << iter_begin->process << endl;

      }
      if( stackmode.compare("Ttbar1lept") == 0 ){

         if( iter_begin->process.compare("TTJets_SemiLept") == 0 ) hists_ = histsMC_signal_;
         else if( iter_begin->process.compare("TTJets_FullLept") == 0 ) hists_ = histsMC_top_dileptonic_;
         else if( iter_begin->process.compare("TTJets_Hadronic") == 0 ) hists_ = histsMC_top_hadronic_;
         else if( iter_begin->process.compare("Tbar_tW") == 0 ) hists_ = histsMC_top_single_;
         else if( iter_begin->process.compare("T_tW") == 0 ) hists_ = histsMC_top_single_;
         else if( iter_begin->process.compare("WJetsToLNu") == 0 ) hists_ = histsMC_EWK_;
         else if( iter_begin->process.compare("DYJetsToLL") == 0 ) hists_ = histsMC_DY_;
         else if( iter_begin->process.compare("DYJetsToLL_M-50") == 0 ) hists_ = histsMC_DY_;
         else if( iter_begin->process.compare("DYJetsToLL_M-10To50") == 0 ) hists_ = histsMC_DY_;
         else if( iter_begin->process.compare("Tbar_tW-channel") == 0 ) hists_ = histsMC_top_single_;
         else if( iter_begin->process.compare("T_tW-channel") == 0 ) hists_ = histsMC_top_single_;
         else if( iter_begin->process.compare("QCD_EMEnriched_20_30") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_EMEnriched_30_80") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_EMEnriched_80_170") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_EMEnriched_170_250") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_EMEnriched_250_350") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_EMEnriched_350") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_BCtoE_20_30") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_BCtoE_30_80") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_BCtoE_80_170") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_BCtoE_170_250") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("QCD_BCtoE_250_350") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->process.compare("Gamma_0_15") == 0 ) hists_ = histsMC_gamma_;
         else if( iter_begin->process.compare("Gamma_15_30") == 0 ) hists_ = histsMC_gamma_;
         else if( iter_begin->process.compare("Gamma_30_50") == 0 ) hists_ = histsMC_gamma_;
         else if( iter_begin->process.compare("Gamma_50_80") == 0 )  hists_ = histsMC_gamma_;
         else if( iter_begin->process.compare("Gamma_80_120") == 0 ) hists_ = histsMC_gamma_;
         else if( iter_begin->process.compare("Gamma_120_170") == 0 ) hists_ = histsMC_gamma_;  
         else if( iter_begin->process.compare("Gamma_170_300") == 0 ) hists_ = histsMC_gamma_;  
         else if( iter_begin->process.compare("Gamma_300_470") == 0 ) hists_ = histsMC_gamma_;  
         else if( iter_begin->process.compare("WW") == 0 ) hists_ = histsMC_EWK_;
         else if( iter_begin->process.compare("WZ") == 0 ) hists_ = histsMC_EWK_;
         else if( iter_begin->process.compare("ZZ") == 0 ) hists_ = histsMC_EWK_;
         else cout << "Histogram fill error, process " << iter_begin->process << endl;

      }

   }

   for( vector<event>::iterator ev = iter_begin; ev < iter_end; ev++ ){

      // test selections
      //if( ev->jet_pt.size() < 2 ) continue;
      //if( ev->jet_pt[0] < 35 or ev->jet_pt[1] < 35 ) continue;

      // this is a hack, need to clean up!
      if( ev->process.compare("DYJetsToLL") == 0 ) hists_ = histsMC_signal_;
      else if( ev->process.compare("TTJets") == 0 ) hists_ = histsMC_top_;
      else if( ev->process.compare("Tbar_tW-channel") == 0 ) hists_ = histsMC_top_;
      else if( ev->process.compare("T_tW-channel") == 0 ) hists_ = histsMC_top_;
      else if( ev->process.compare("WJetsToLNu") == 0 ) hists_ = histsMC_EWK_;
      else if( ev->process.compare("WW") == 0 ) hists_ = histsMC_EWK_;
      else if( ev->process.compare("WZ") == 0 ) hists_ = histsMC_EWK_;
      else if( ev->process.compare("ZZ") == 0 ) hists_ = histsMC_EWK_;

      // jets
      int njets = 0;
      for( int j=0; j < int(ev->jet_pt.size()); j++){
         //if( fabs(ev->jet_eta[j]) < 2.5 and ev->jet_pt[j] > 30 ){
            hists_["jet_pt"]->Fill( ev->jet_pt[j], ev->weight );
            hists_["jet_eta"]->Fill( ev->jet_eta[j], ev->weight );
            if( j == 0 ) hists_["jet1_pt"]->Fill( ev->jet_pt[j] , ev->weight );
            njets++;
         //}
      }
      hists_["njets"]->Fill( njets, ev->weight );

      // leptons
      for( int j=0; j < int(ev->lepton_pt.size()); j++){
         hists_["lepton_pt"]->Fill( ev->lepton_pt[j] , ev->weight );
         hists_["lepton_eta"]->Fill( ev->lepton_eta[j] , ev->weight );
      }
      hists_["muon_invmass"]->Fill( ev->invmass, ev->weight );
      // other observables
      //if( njets >= 2 ){
         hists_["met"]->Fill( ev->met_pt, ev->weight );
         hists_["met_200"]->Fill( ev->met_pt, ev->weight );
         hists_["met_300"]->Fill( ev->met_pt, ev->weight );
         hists_["sig"]->Fill( ev->sig, ev->weight );
         hists_["sig_100"]->Fill( ev->sig, ev->weight );
         hists_["sig_15"]->Fill( ev->sig, ev->weight );
         hists_["pchi2"]->Fill( TMath::Prob(ev->sig,2), ev->weight );
      //}
      hists_["nvert"]->Fill( ev->nvertices , ev->weight );
      hists_["pjet_scalpt"]->Fill( ev->pjet_scalpt , ev->weight );

      /*
      // pseudojet
      hists_["pjet_scalptL123"]->Fill( ev->pjet_scalptL123 , ev->weight );
      hists_["pjet_vectptL123"]->Fill( ev->pjet_vectptL123 , ev->weight );
      hists_["pjet_scalptT1"]->Fill( ev->pjet_scalptT1 , ev->weight );
      hists_["pjet_vectptT1"]->Fill( ev->pjet_vectptT1 , ev->weight );
      hists_["pjet_size"]->Fill( ev->pjet_size , ev->weight );
      hists_["pjet_phi"]->Fill( ev->pjet_phiL123, ev->weight );
      */

      hists_["qt"]->Fill( ev->qt, ev->weight );
      hists_["qx"]->Fill( ev->qx, ev->weight );
      hists_["ut_par"]->Fill( fabs(ev->ut_par), ev->weight );

      hists_["cov_xx"]->Fill( ev->cov_xx, ev->weight );
      hists_["cov_xy"]->Fill( ev->cov_xy, ev->weight );
      hists_["cov_yy"]->Fill( ev->cov_yy, ev->weight );
      /*
      hists_["met0_px"]->Fill( ev->pfmet_px[0], ev->weight );
      hists_["met1_px"]->Fill( ev->pfmet_px[1], ev->weight );
      hists_["met2_px"]->Fill( ev->pfmet_px[2], ev->weight );
      hists_["met3_px"]->Fill( ev->pfmet_px[3], ev->weight );
      hists_["met4_px"]->Fill( ev->pfmet_px[4], ev->weight );
      hists_["met5_px"]->Fill( ev->pfmet_px[5], ev->weight );
      */
      hists_["sig_old"]->Fill( ev->metsig2011, ev->weight );
      hists_["det"]->Fill( ev->det, ev->weight );
      hists_["pchi2_old"]->Fill( TMath::Prob(ev->metsig2011,2), ev->weight );
      hists_["logpchi2"]->Fill( TMath::Log(TMath::Prob(ev->sig,2)), ev->weight );

      hists_["cov_xx_highpt"]->Fill( ev->cov_xx_highpt, ev->weight );
      hists_["cov_xx_pjet"]->Fill( ev->cov_xx_pjet, ev->weight );
      hists_["cov_xx_ratio"]->Fill( ev->cov_xx_highpt/ev->cov_xx, ev->weight );

      /*
      hists_["met_varx"]->Fill( ev->met_varx, ev->weight );
      hists_["met_vary"]->Fill( ev->met_vary, ev->weight );
      hists_["met_rho"]->Fill( ev->met_rho, ev->weight );
      */

      // profiles
      profs_["psig_nvert"]->Fill( ev->nvertices, ev->sig, ev->weight );
      profs_["psig_qt"]->Fill( ev->qt, ev->sig, ev->weight );
      profs_["presp_qt"]->Fill( ev->qt, -(ev->ut_par)/(ev->qt), ev->weight );
      profs_["pMET_nvert"]->Fill( ev->nvertices, ev->met_pt, ev->weight );
      //profs_["pjet_scalptL123_nvert"]->Fill( ev->nvertices, ev->pjet_scalptL123, ev->weight );

      profs_["njets_nvert"]->Fill( ev->nvertices, ev->jet_pt.size(), ev->weight );
      for( int j=0; j < int(ev->jet_pt.size()); j++){
         profs_["jet_pt_nvert"]->Fill( ev->nvertices, ev->jet_pt[j], ev->weight );
      }

      profs_["sig_met"]->Fill( ev->met_pt, ev->sig, ev->weight );
   }

}

void Fitter::PrintHists( const char* filename, string stackmode, bool overflow ){

   // draw hists and write to file
   gROOT->ProcessLineSync(".L tdrstyle.C");
   gROOT->ProcessLineSync("setTDRStyle()");

   TFile *file = new TFile(filename,"RECREATE");
   TDirectory* th2=file->mkdir("2D Hists");
   file->cd();

   //
   // determine mc hist normalization
   //

   // initial normalization
   TH1D *histData_temp = (TH1D*)histsData_["met"]->Clone("histData_temp");
   TH1D *histMC_signal_temp = (TH1D*)histsMC_signal_["met"]->Clone("histMC_signal_temp");
   TH1D *histMC_top_temp = (TH1D*)histsMC_top_["met"]->Clone("histMC_top_temp");
   TH1D *histMC_top_dileptonic_temp = (TH1D*)histsMC_top_dileptonic_["met"]->Clone("histMC_top_dileptonic_temp");
   TH1D *histMC_top_hadronic_temp = (TH1D*)histsMC_top_hadronic_["met"]->Clone("histMC_top_hadronic_temp");
   TH1D *histMC_top_single_temp = (TH1D*)histsMC_top_single_["met"]->Clone("histMC_top_single_temp");
   TH1D *histMC_EWK_temp = (TH1D*)histsMC_EWK_["met"]->Clone("histMC_EWK_temp");
   TH1D *histMC_QCD_temp = (TH1D*)histsMC_QCD_["met"]->Clone("histMC_QCD_temp");
   TH1D *histMC_gamma_temp = (TH1D*)histsMC_gamma_["met"]->Clone("histMC_gamma_temp");
   TH1D *histMC_DY_temp = (TH1D*)histsMC_DY_["met"]->Clone("histMC_DY_temp");

   // rescale QCD & gamma+jets numerically (approximate)
   double chi2 = -1;
   double histnorm = histData_temp->Integral("width") / (histMC_signal_temp->Integral("width")+histMC_top_temp->Integral("width")+histMC_top_dileptonic_temp->Integral("width")+histMC_top_hadronic_temp->Integral("width")+histMC_top_single_temp->Integral("width")+histMC_EWK_temp->Integral("width")+histMC_QCD_temp->Integral("width")+histMC_gamma_temp->Integral("width")+histMC_DY_temp->Integral("width"));
   double scaleQCD = 1;
   if( stackmode.compare("Wenu") == 0 or stackmode.compare("Wenu_loose") == 0 
        /* or stackmode.compare("Ttbar0lept") == 0 */){ 
      for(double s = 0; s < 5; s += 0.01){
         TH1D *histMC_temp = new TH1D( "histMC_temp", "histMC_temp",
               histData_temp->GetNbinsX(), histData_temp->GetBinLowEdge(1),
               histData_temp->GetBinLowEdge(histData_temp->GetNbinsX())
               + histData_temp->GetBinWidth(histData_temp->GetNbinsX()) );

         histMC_temp->Add( histMC_signal_temp );
         histMC_temp->Add( histMC_top_temp );
         histMC_temp->Add( histMC_top_dileptonic_temp );
         histMC_temp->Add( histMC_top_hadronic_temp );
         histMC_temp->Add( histMC_top_single_temp );
         histMC_temp->Add( histMC_EWK_temp );
         histMC_temp->Add( histMC_QCD_temp, s );
         histMC_temp->Add( histMC_gamma_temp, s );
         histMC_temp->Add( histMC_DY_temp );

         double histnorm_temp = histData_temp->Integral("width") / histMC_temp->Integral("width");
         histMC_temp->Scale( histnorm_temp );
         double chi2_temp = histData_temp->Chi2Test( histMC_temp, "CHI2" );

         if( chi2 < 0 or chi2_temp < chi2 ){
            chi2 = chi2_temp;
            scaleQCD = s;
            histnorm = histnorm_temp;
         }

         delete histMC_temp;
      }
      cout << "Found scale factor = " << scaleQCD << " with chi2 = " << chi2 << endl;
   }

   //
   // fill histograms
   //
   for(map<string,TH1*>::const_iterator it = histsData_.begin();
         it != histsData_.end(); it++){

      string hname = it->first;
      TH1D *histData = (TH1D*)it->second;
      TH1D *histMC = (TH1D*)histsMC_[hname];

      TH1D *histMC_signal = (TH1D*)histsMC_signal_[hname];
      TH1D *histMC_top = (TH1D*)histsMC_top_[hname];
      TH1D *histMC_top_dileptonic = (TH1D*)histsMC_top_dileptonic_[hname];
      TH1D *histMC_top_hadronic = (TH1D*)histsMC_top_hadronic_[hname];
      TH1D *histMC_top_single = (TH1D*)histsMC_top_single_[hname];
      TH1D *histMC_EWK = (TH1D*)histsMC_EWK_[hname];
      TH1D *histMC_QCD = (TH1D*)histsMC_QCD_[hname];
      TH1D *histMC_gamma = (TH1D*)histsMC_gamma_[hname];
      TH1D *histMC_DY = (TH1D*)histsMC_DY_[hname];

      // scale QCD
      histMC_QCD->Scale( scaleQCD );
      histMC_gamma->Scale( scaleQCD );

      // add overflow bin
      if( overflow ){
         AddOverflow(histData);
         AddOverflow(histMC);
         AddOverflow(histMC_signal);
         AddOverflow(histMC_top);
         AddOverflow(histMC_top_dileptonic);
         AddOverflow(histMC_top_hadronic);
         AddOverflow(histMC_top_single);
         AddOverflow(histMC_EWK);
         AddOverflow(histMC_QCD);
         AddOverflow(histMC_gamma);
         AddOverflow(histMC_DY);
      }

      // get total MC histogram
      histMC->Add( histMC_signal );
      histMC->Add( histMC_top );
      histMC->Add( histMC_top_dileptonic );
      histMC->Add( histMC_top_hadronic );
      histMC->Add( histMC_top_single );
      histMC->Add( histMC_EWK );
      histMC->Add( histMC_QCD );
      histMC->Add( histMC_gamma );
      histMC->Add( histMC_DY );

      histMC->Scale( histnorm );
      histMC_signal->Scale( histnorm );
      histMC_QCD->Scale( histnorm );
      histMC_gamma->Scale( histnorm );
      histMC_DY->Scale( histnorm );
      histMC_top->Scale( histnorm );
      histMC_top_dileptonic->Scale( histnorm );
      histMC_top_hadronic->Scale( histnorm );
      histMC_top_single->Scale( histnorm );
      histMC_EWK->Scale( histnorm );

      bool log_axis = (hname != "jet_eta") and (hname != "pchi2") and (hname != "pchi2_old");

      TCanvas *canvas  = new TCanvas( (char*)hname.c_str(), (char*)hname.c_str(), 800, 800 );
      canvas->SetFillColor(0);
      TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
      TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
      pad1->SetTopMargin(0.1);
      pad1->SetBottomMargin(0.01);
      pad1->SetRightMargin(0.1);
      pad1->SetFillColor(0);
      if( log_axis ) pad1->SetLogy();
      pad2->SetTopMargin(0.01);
      pad2->SetBottomMargin(0.3);
      pad2->SetRightMargin(0.1);
      pad2->SetFillColor(0);
      pad1->Draw();
      pad2->Draw();

      pad1->cd();

      histMC->SetLineColor(1);
      histMC->SetFillColor(kAzure-9);
      histData->SetLineColor(1);
      histData->SetMarkerStyle(20);

      histMC->SetMaximum( 1.1*max(histMC->GetMaximum(), histData->GetMaximum()) );
      histMC->SetMinimum( 0.1*min(histMC->GetMinimum(), histData->GetMinimum()) );
      if( log_axis and histMC->GetMinimum() < 1.0e-06 ) histMC->SetMinimum( 1.0e-06 );
      histMC->SetLabelSize(0.0);
      histMC->GetXaxis()->SetTitleSize(0.00);
      histMC->GetYaxis()->SetLabelSize(0.07);
      histMC->GetYaxis()->SetTitleSize(0.08);
      histMC->GetYaxis()->SetTitleOffset(0.76);
      histMC->GetXaxis()->SetLabelFont(42);
      histMC->GetYaxis()->SetLabelFont(42);
      histMC->GetXaxis()->SetTitleFont(42);
      histMC->GetYaxis()->SetTitleFont(42);

      histMC->Draw("HIST");
      if( stackmode.compare("Zmumu") == 0 ){
         histMC_top->SetLineColor(1);
         histMC_top->SetFillColor(kYellow-9);
         histMC_EWK->SetLineColor(1);
         histMC_EWK->SetFillColor(kRed-10);

         histMC_top->Add( histMC_EWK );

         histMC_top->Draw("same HIST");
         histMC_EWK->Draw("same HIST");
      }
      if( stackmode.compare("Wenu") == 0 or stackmode.compare("Wenu_loose") == 0 ){
         histMC_QCD->SetLineColor(1);
         histMC_QCD->SetFillColor(kYellow-9);
         histMC_gamma->SetLineColor(1);
         histMC_gamma->SetFillColor(kYellow-10);
         histMC_DY->SetLineColor(1);
         histMC_DY->SetFillColor(kRed-10);
         histMC_top->SetLineColor(1);
         histMC_top->SetFillColor(kCyan-10);
         histMC_EWK->SetLineColor(1);
         histMC_EWK->SetFillColor(kGreen-10);

         histMC_QCD->Add( histMC_gamma );
         histMC_QCD->Add( histMC_DY );
         histMC_QCD->Add( histMC_top );
         histMC_QCD->Add( histMC_EWK );
         histMC_gamma->Add( histMC_DY );
         histMC_gamma->Add( histMC_top );
         histMC_gamma->Add( histMC_EWK );
         histMC_DY->Add( histMC_top );
         histMC_DY->Add( histMC_EWK );
         histMC_top->Add( histMC_EWK );

         histMC_QCD->Draw("same HIST");
         histMC_gamma->Draw("same HIST");
         histMC_DY->Draw("same HIST");
         histMC_top->Draw("same HIST");
         histMC_EWK->Draw("same HIST");
      }
      if( stackmode.compare("Ttbar1lept") == 0 ){
         histMC_top_dileptonic->SetLineColor(1);
         histMC_top_dileptonic->SetFillColor(kYellow-9);
         histMC_top_hadronic->SetLineColor(1);
         histMC_top_hadronic->SetFillColor(kYellow-3);
         histMC_top_single->SetLineColor(1);
         histMC_top_single->SetFillColor(kOrange-4);
         histMC_DY->SetLineColor(1);
         histMC_DY->SetFillColor(kRed-10);
         histMC_EWK->SetLineColor(1);
         histMC_EWK->SetFillColor(kCyan-10);

         histMC_top_dileptonic->Add( histMC_top_hadronic );
         histMC_top_dileptonic->Add( histMC_top_single );
         histMC_top_dileptonic->Add( histMC_DY );
         histMC_top_dileptonic->Add( histMC_EWK );
         histMC_top_hadronic->Add( histMC_top_single );
         histMC_top_hadronic->Add( histMC_DY );
         histMC_top_hadronic->Add( histMC_EWK );
         histMC_top_single->Add( histMC_DY );
         histMC_top_single->Add( histMC_EWK );
         histMC_DY->Add( histMC_EWK );

         histMC_top_dileptonic->Draw("same HIST");
         histMC_top_hadronic->Draw("same HIST");
         histMC_top_single->Draw("same HIST");
         histMC_DY->Draw("same HIST");
         histMC_EWK->Draw("same HIST");
      }
      if( stackmode.compare("Ttbar0lept") == 0 ){
         histMC_top->SetLineColor(1);
         histMC_top->SetFillColor(kYellow-9);
         histMC_QCD->SetLineColor(1);
         histMC_QCD->SetFillColor(kRed-10);
         histMC_DY->SetLineColor(1);
         histMC_DY->SetFillColor(kCyan-10);

         histMC_top->Add( histMC_QCD );
         histMC_top->Add( histMC_DY );
         histMC_QCD->Add( histMC_DY );

         histMC_top->Draw("same HIST");
         histMC_QCD->Draw("same HIST");
         histMC_DY->Draw("same HIST");
      }
      if( stackmode.compare("Dijet") == 0 ){
         histMC_top->SetLineColor(1);
         histMC_top->SetFillColor(kYellow-9);
         histMC_EWK->SetLineColor(1);
         histMC_EWK->SetFillColor(kRed-10);
         histMC_DY->SetLineColor(1);
         histMC_DY->SetFillColor(kCyan-10);

         histMC_top->Add( histMC_EWK );
         histMC_top->Add( histMC_DY );
         histMC_EWK->Add( histMC_DY );

         histMC_top->Draw("same HIST");
         histMC_EWK->Draw("same HIST");
         histMC_DY->Draw("same HIST");
      }

      TH1D *histMCerror = (TH1D*)histMC->Clone("histMCerror");
      histMCerror->SetFillColor(1);
      histMCerror->SetFillStyle(3005);
      histMCerror->Draw("E2 same");

      histData->Draw("EP same");

      TLegend *leg = new TLegend(0.605528,0.655866,0.866834,0.816333);
      leg->AddEntry(histData, "data");
      if( stackmode.compare("Zmumu") == 0 ){
         leg->AddEntry(histMC, "Z #rightarrow #mu #mu", "f");
         leg->AddEntry(histMC_top, "t #bar{t}", "f");
         leg->AddEntry(histMC_EWK, "EWK", "f");
      }
      if( stackmode.compare("Wenu") == 0 or stackmode.compare("Wenu_loose") == 0 ){
         leg->AddEntry(histMC, "W #rightarrow e #nu", "f");
         leg->AddEntry(histMC_QCD, "QCD", "f");
         leg->AddEntry(histMC_gamma, "#gamma+jets", "f");
         leg->AddEntry(histMC_DY,  "DY", "f");
         leg->AddEntry(histMC_top, "t #bar{t}", "f");
         leg->AddEntry(histMC_EWK, "EWK", "f");
      }
      if( stackmode.compare("Dijet") == 0 ){
         leg->AddEntry(histMC, "QCD", "f");
         leg->AddEntry(histMC_top, "top", "f");
         leg->AddEntry(histMC_EWK, "EWK", "f");
         leg->AddEntry(histMC_DY, "DY", "f");
      }
      if( stackmode.compare("Ttbar1lept") == 0 ){
         leg->AddEntry(histMC, "t#bar{t} Semi-Leptonic", "f");
         leg->AddEntry(histMC_top_dileptonic, "Dileptonic t#bar{t}", "f");
         leg->AddEntry(histMC_top_hadronic, "Hadronic t#bar{t}", "f");
         leg->AddEntry(histMC_top_single, "Single Top", "f");
         leg->AddEntry(histMC_DY, "DY", "f");
         leg->AddEntry(histMC_EWK, "EWK", "f");
      }
      if( stackmode.compare("Ttbar0lept") == 0 ){
         leg->AddEntry(histMC, "t#bar{t} Hadronic", "f");
         leg->AddEntry(histMC_top, "Other Top", "f");
         leg->AddEntry(histMC_QCD, "QCD", "f");
         leg->AddEntry(histMC_DY,  "DY", "f");
      }
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->Draw("same");

      pad2->cd();

      TH1D *hratio = (TH1D*)histData->Clone("hratio");
      hratio->Divide( histMC );
      hratio->SetTitle(";"+TString(histData->GetXaxis()->GetTitle()));
      hratio->SetStats(0);

      hratio->GetXaxis()->SetTitleSize(0.14);
      hratio->GetXaxis()->SetLabelSize(0.14);
      hratio->GetYaxis()->SetLabelSize(0.11);
      hratio->GetYaxis()->SetTitleSize(0.14);
      hratio->GetYaxis()->SetTitleOffset(0.28);
      hratio->GetXaxis()->SetLabelFont(42);
      hratio->GetYaxis()->SetLabelFont(42);
      hratio->GetXaxis()->SetTitleFont(42);
      hratio->GetYaxis()->SetTitleFont(42);
      hratio->SetMaximum( 1.6 );
      hratio->SetMinimum( 0.4 );
      hratio->GetYaxis()->SetNdivisions(505);
      hratio->Draw("EP");

      TF1 *func = new TF1("func","[0]",-10E6,10E6);
      func->SetParameter(0,1.0);
      func->SetLineWidth(1);
      func->SetLineStyle(7);
      func->SetLineColor(1);
      func->Draw("same");

      canvas->cd();
      canvas->Write();

      delete pad1;
      delete pad2;
      delete leg;
      delete func;
      delete canvas;
   }
   cout << "\nCORRELATION COEFFICIENTS: data, mc" << endl;
   for(map<string,TH2*>::const_iterator it = profsData_.begin();
         it != profsData_.end(); it++){

      string pname = it->first;
      TH2 *histData = (TH2*)it->second;
      TH2 *histMC = (TH2*)profsMC_[pname];

      th2->cd();
      histData->Write();
      histMC->Write();
      file->cd();

      //preserve y-axis labels after converting to profile
      const char* labelData=histData->GetYaxis()->GetTitle();
      const char* labelMC=histMC->GetYaxis()->GetTitle();

      TProfile *profData = (TProfile*)histData->ProfileX();
      TProfile *profMC = (TProfile*)histMC->ProfileX();

      profData->GetYaxis()->SetTitle(labelData);
      profMC->GetYaxis()->SetTitle(labelMC);

      TCanvas *canvas  = new TCanvas( (char*)pname.c_str(), (char*)pname.c_str(), 800, 800 );
      canvas->cd();

      profMC->SetLineColor(2);
      profData->SetLineColor(1);
      profData->SetMarkerStyle(20);

      profMC->SetMaximum( 1.1*max(profMC->GetMaximum(), profData->GetMaximum()) );
      profMC->SetMinimum( 0.8*min(profMC->GetMinimum(), profData->GetMinimum()) );

      profMC->GetXaxis()->SetLabelFont(42);
      profMC->GetXaxis()->SetLabelSize(0.05);
      profMC->GetXaxis()->SetTitleFont(42);
      profMC->GetXaxis()->SetTitleSize(0.06);

      profMC->Draw("HIST");
      profData->Draw("EP same");

      TLegend *leg = new TLegend(0.605528,0.655866,0.866834,0.816333);
      leg->AddEntry(profData, "data");
      leg->AddEntry(profMC, "MC");
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->Draw("same");

      canvas->Write();

      cout << pname << ": " << histData->GetCorrelationFactor() << ", "
         << histMC->GetCorrelationFactor() << endl;

      delete leg;
      delete canvas;
   }
   TCanvas *canvas2 = new TCanvas( "sig_met_2D", "sig_met_2D", 800, 800 );
   canvas2->cd();
   profsData_["sig_met"]->Draw("colz");
   canvas2->Write();

   TF1 *pchi2_left = new TF1("pchi2_left","pol1",0.0,0.03);
   histsData_["pchi2"]->Fit("pchi2_left","QR");
   pchi2slope_left = pchi2_left->GetParameter(1);

   TF1 *pchi2_right = new TF1("pchi2_right","pol1",0.5,1.0);
   histsData_["pchi2"]->Fit("pchi2_right","QR");
   pchi2slope_right = pchi2_right->GetParameter(1);

   cout << "\nP(CHI2) SLOPES: " << pchi2slope_left << " L, " << pchi2slope_right << " R\n" << endl;

   psig_nvert_corr = profsData_["psig_nvert"]->GetCorrelationFactor();
   psig_qt_corr = profsData_["psig_qt"]->GetCorrelationFactor();

   file->Close();

   delete canvas2;
   delete pchi2_left;
   delete pchi2_right;
   delete file;
}

void Fitter::AddOverflow( TH1* hist ){

   int nbins = hist->GetNbinsX();
   double overflow = hist->GetBinContent(nbins+1);

   hist->SetBinContent( nbins, hist->GetBinContent(nbins) + overflow );

}
