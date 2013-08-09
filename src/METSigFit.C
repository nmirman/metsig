#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TVector2.h"
#include "TRandom3.h"
#include "TProfile.h"
#include "TString.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"

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

Fitter::Fitter(){
   // MINUIT variables
   gMinuit = 0;
   fFunc = 0;

   // jet pt bounds
   jetbinpt = 20.0;
   jetcorrpt = 10.0;

   // significance cut for minimization
   significance_cut = false;

   // MET type to use in FindSignificance
   met_type = -1;

   // compute weighted errors
   TH1::SetDefaultSumw2();

   // data hists for plotting
   histsData_["muon_pt"] = new TH1D("muon_pt_Data", "Muon p_{T}", 50, 0, 200);
   histsData_["muon_invmass"] = new TH1D("muon_invmass_Data", "M_{#mu#mu}", 50, 0, 200);
   histsData_["njets"  ] = new TH1D("njets_Data", "N jets", 100, 0, 100);
   histsData_["jet_pt" ] = new TH1D("jet_pt_Data", "Jet p_{T}", 50, 0, 200);
   histsData_["jet1_pt"] = new TH1D("jet1_pt_Data", "Jet 1 p_{T}", 50, 0, 200);
   histsData_["jet_eta" ] = new TH1D("jet_eta_Data", "Jet #eta, p_{T} > 30 GeV", 50, -5, 5);
   histsData_["pjet_size"  ] = new TH1D("pjet_size_Data", "N jets", 50, 0, 500);
   histsData_["pjet_scalptL123"] = new TH1D("pjet_scalptL123_Data", "Pseudojet Scalar p_{T} (L123-corrected)", 200, 0, 2000);
   histsData_["pjet_vectptL123"] = new TH1D("pjet_vectptL123_Data", "Pseudojet Scalar p_{T} (L123-corrected)", 50, 0, 200);
   histsData_["pjet_scalptT1"] = new TH1D("pjet_scalptT1_Data", "Pseudojet Scalar p_{T} (T1-corrected)", 200, 0, 2000);
   histsData_["pjet_vectptT1"] = new TH1D("pjet_vectptT1_Data", "Pseudojet Scalar p_{T} (T1-corrected)", 50, 0, 200);
   histsData_["pjet_phi"] = new TH1D("pjet_phi_Data", "Pseudojet Phi", 50, -3.5, 3.5);
   histsData_["qt"] = new TH1D("qt_Data", "q_{T}", 50, 0, 200);
   histsData_["ut_par"] = new TH1D("ut_par_Data", "|u_{T}|_{#parallel}", 50, 0, 100);
   histsData_["nvert"] = new TH1D("nvert_Data", "N Vertices", 50, 0, 100);
   histsData_["cov_xx"] = new TH1D("cov_xx_Data", "Cov_{xx}", 50, 0, 500);
   histsData_["cov_xy"] = new TH1D("cov_xy_Data", "Cov_{xy}", 50, -150, 150);
   histsData_["cov_yy"] = new TH1D("cov_yy_Data", "Cov_{yy}", 50, 0, 500);
   histsData_["met"] = new TH1D("met_Data", "Missing E_{T}", 50, 0, 100);
   histsData_["met_200"] = new TH1D("met_200_Data", "Missing E_{T}", 50, 0, 200);
   histsData_["sig"] = new TH1D("sig_Data", "Significance", 50, 0, 500);
   histsData_["sig_100"] = new TH1D("sig_100_Data", "Significance", 50, 0, 100);
   histsData_["sig_15"] = new TH1D("sig_15_Data", "Significance", 50, 0, 15);
   histsData_["sig_old"] = new TH1D("sig_old_Data", "Old Significance", 50, 0, 500);
   histsData_["det"] = new TH1D("det_Data", "Determinant", 50, 0, 100000);
   histsData_["pchi2"] = new TH1D("pchi2_Data", "P(#chi^{2})", 50, 0, 1);
   histsData_["pchi2_old"] = new TH1D("pchi2_old-Data", "P(#chi^{2}) from Old Sig", 50, 0, 1);
   histsData_["logpchi2"] = new TH1D("logpchi2_Data", "log(P(#chi^{2}))", 50, -50, 0);
   histsData_["cov_xx_highpt"] = new TH1D("cov_xx_highpt_Data", "Cov_{xx} High-p_{T} Jets", 50, 0, 500);
   histsData_["cov_xx_pjet"] = new TH1D("cov_xx_pjet_Data", "Cov_{xx} Pseudojet", 50, 0, 500);
   histsData_["cov_xx_ratio"] = new TH1D("cov_xx_ratio_Data", "Cov_{xx} High-p_{T}/Total", 50, 0, 1 );
   histsData_["met_varx"] = new TH1D("met_varx_Data","#sigma_{x}^{2} for MET smearing",50,-10,30);
   histsData_["met_vary"] = new TH1D("met_vary_Data","#sigma_{y}^{2} for MET smearing",50,-10,30);
   histsData_["met_rho"] = new TH1D("met_rho_Data", "#rho for MET smearing", 50, -10, 10);
   histsData_["pjet_axesratio"] = new TH1D("pjet_axesratio_Data",
         "Pseudo-jet Ratio of Major/Semi-Major Axes", 50, 0, 1);
   histsData_["pjet_tiltangle"] = new TH1D("pjet_tiltangle_Data", "Tilt Angle of Pseudo-jet",
         50, -2, 2);
   histsData_["pjet_tiltangle_rel"] = new TH1D("pjet_tiltangle_rel_Data",
      "Tilt Angle of Pseudo-jet relative to Pseudo-jet Momentum", 50, -2, 2);

   // profile histograms
   profsData_["psig_nvert"] = new TH2D("psig_nvert_Data",
         "Significance vs. N Vertices;N Vertices;<S_{E}>", 30, 0, 30, 50, 0, 50);
   profsData_["psig_qt"] = new TH2D("psig_qt_Data",
         "Significance vs. q_{T};q_{T} (GeV);<S_{E}>", 15, 0, 50, 100, 0, 50);
   profsData_["presp_qt"] = new TH2D("presp_qt_Data",
         "Response = |<u_{#parallel}>|/q_{T} vs. q_{T};q_{T} (GeV);Response", 25, 0, 50, 100, -100, 100);
   profsData_["pMET_nvert"] = new TH2D("pMET_nvert_Data",
         "MET vs. N Vertices;N Vertices;<MET>", 30, 0, 30, 50, 0, 100);
   profsData_["pjet_scalptL123_nvert"] = new TH2D("pjet_scalpt_nvert_Data",
         "Pseudojet Scalar p_{T} vs. N Vertices", 30, 0, 30, 200, 0, 2000);
   profsData_["jet_pt_nvert"] = new TH2D("jet_pt_nvert_Data",
         "Jet p_{T} vs. N Vertices", 30, 0, 30, 50, 0, 200);
   profsData_["njets_nvert"] = new TH2D("njets_nvert_Data",
         "N jets vs. N Vertices", 30, 0, 30, 50, 0, 100);
   profsData_["sig_met"] = new TH2D("sig_met_Data",
         "Significance vs. MET", 500, 0, 100, 500, 0, 30);

   // clone data hists for MC
   for(map<string,TH1*>::const_iterator it = histsData_.begin();
         it != histsData_.end(); it++){

      string hname = it->first;
      TH1D *hist = (TH1D*)it->second;
      histsMC_[hname] = (TH1D*)hist->Clone( (hname+"_MC").c_str() );

      histsMC_signal_[hname] = (TH1D*)hist->Clone( (hname+"_MC_signal").c_str() );
      histsMC_top_[hname] = (TH1D*)hist->Clone( (hname+"_MC_top").c_str() );
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

void Fitter::ReadNtuple(const char* filename, vector<event>& eventref_temp, const double fracevents,
      const bool isMC, string channel, const bool do_resp_correction, 
      const int start_evt_num, const int end_evt_num ){
   cout << "---> ReadNtuple " << channel << endl;

   float gi_xsec=0;
   float gi_eff=0;

   int v_size;
   float met_pt[5];
   float met_px[5];
   float met_py[5];
   float metsig2011[5];

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
   float mu_eta[100];
   int mu_isGlobal[100];
   float mu_chi2[100];
   int mu_muonHits[100];
   int mu_nMatches[100];
   float mu_dxy[100];
   int mu_pixelHits[100];
   int mu_numberOfValidTrackerLayers[100];
   float mu_dr03TkSumPt[100];
   float mu_dr03EcalRecHitSumEt[100];
   float mu_dr03HcalTowerSumEt[100];

   int elec_size=0;
   float elec_pt[100];
   float elec_px[100];
   float elec_py[100];
   float elec_pz[100];
   float elec_e[100];
   float elec_phi[100];
   float elec_eta[100];

   float pfj_l1[1000];
   float pfj_l1l2l3[1000];

   float puMyWeight = 1;

   int genj_size=0;
   float genj_pt[1000];
   float genj_phi[1000];
   float genj_eta[1000];
   float genj_energy[1000];
   float genj_emEnergy[1000];
   float genj_hadEnergy[1000];
   float genj_invEnergy[1000];

   //TFile *file = TFile::Open(filename);
   //if( !file ){ return; }

   //TTree *tree = (TTree*)file->Get("events");
   TChain *tree = new TChain("events");
   tree->Add( filename );

   tree->SetBranchAddress("v_size", &v_size);
   tree->SetBranchAddress("met_pt", met_pt);
   tree->SetBranchAddress("met_px", met_px);
   tree->SetBranchAddress("met_py", met_py);
   tree->SetBranchAddress("met_sig", metsig2011);

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
   tree->SetBranchAddress("mu_eta", mu_eta);

   tree->SetBranchAddress("mu_isGlobal", &mu_isGlobal);
   tree->SetBranchAddress("mu_chi2", &mu_chi2);
   tree->SetBranchAddress("mu_muonHits", &mu_muonHits);
   tree->SetBranchAddress("mu_nMatches", &mu_nMatches);
   tree->SetBranchAddress("mu_dxy", &mu_dxy);
   tree->SetBranchAddress("mu_pixelHits", &mu_pixelHits);
   tree->SetBranchAddress("mu_numberOfValidTrackerLayers", &mu_numberOfValidTrackerLayers);
   tree->SetBranchAddress("mu_dr03TkSumPt", &mu_dr03TkSumPt);
   tree->SetBranchAddress("mu_dr03EcalRecHitSumEt", &mu_dr03EcalRecHitSumEt);
   tree->SetBranchAddress("mu_dr03HcalTowerSumEt", &mu_dr03HcalTowerSumEt);

   tree->SetBranchAddress("elec_size", &elec_size);
   tree->SetBranchAddress("elec_pt", elec_pt);
   tree->SetBranchAddress("elec_px", elec_px);
   tree->SetBranchAddress("elec_py", elec_py);
   tree->SetBranchAddress("elec_pz", elec_pz);
   tree->SetBranchAddress("elec_e", elec_e);
   tree->SetBranchAddress("elec_phi", elec_phi);
   tree->SetBranchAddress("elec_eta", elec_eta);

   tree->SetBranchAddress("pfj_l1", pfj_l1);
   tree->SetBranchAddress("pfj_l1l2l3", pfj_l1l2l3);

   if(isMC){
      tree->SetBranchAddress("gi_xsec", &gi_xsec);
      tree->SetBranchAddress("gi_eff", &gi_eff);
      tree->SetBranchAddress("puMyWeight", &puMyWeight);
      tree->SetBranchAddress("genj_size", &genj_size);
      tree->SetBranchAddress("genj_pt", genj_pt);
      tree->SetBranchAddress("genj_phi", genj_phi);
      tree->SetBranchAddress("genj_eta", genj_eta);
      tree->SetBranchAddress("genj_energy", genj_energy);
      tree->SetBranchAddress("genj_emEnergy", genj_emEnergy);
      tree->SetBranchAddress("genj_hadEnergy", genj_hadEnergy);
      tree->SetBranchAddress("genj_invEnergy", genj_invEnergy);
   }

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

      int pjet_size_temp = 0;
      double pjet_scalptL123_temp = 0;
      double pjet_pxL123_temp = 0;
      double pjet_pyL123_temp = 0;
      double pjet_scalptT1_temp = 0;
      double pjet_pxT1_temp = 0;
      double pjet_pyT1_temp = 0;
      double pjet_pxpx_temp = 0;
      double pjet_pypy_temp = 0;
      double pjet_pxpy_temp = 0;

      event evtemp;

      evtemp.channel = channel;
      evtemp.weight = 1.0;

      // MC cross-sections (pb)
      // obtained at https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV
      // and AN-12-333 on MET in Zmumu channel (8 TeV).
      double xsec_dyjetstoll           = 3351.97;
      double xsec_dyjetstoll_m10to50   = 860.5;
      double xsec_ttjets               = 234;
      double xsec_ttjets_0lept         = 0.46*234;
      double xsec_ttjets_1lept         = 0.45*234;
      double xsec_ttjets_2lept         = 0.09*234;
      double xsec_ww                   = 54.838;
      double xsec_wz                   = 33.21;
      double xsec_zz                   = 8.059;
      double xsec_tbar_tw              = 11.1;
      double xsec_t_tw                 = 11.1;
      double xsec_wjetstolnu           = 37509.0;

      // MC sample size (events), obtained from DAS.
      int nevts_dyjetstoll          = 30459503;
      int nevts_dyjetstoll_m10to50  = 7132223;
      int nevts_ttjets              = 6923750;
      int nevts_ttjets_0lept        = 31223821;
      int nevts_ttjets_1lept        = 25424818;
      int nevts_ttjets_2lept        = 12119013;
      int nevts_qcd_15_30           = 9987968;
      int nevts_qcd_30_50           = 413184;
      int nevts_qcd_50_80           = 385398;
      int nevts_qcd_80_120          = 845280;
      int nevts_qcd_120_170         = 885762;
      int nevts_qcd_170_300         = 296948;
      int nevts_qcd_300_470         = 317900;
      int nevts_qcd_470_600         = 328952;
      int nevts_qcd_600_800         = 269784;
      int nevts_qcd_800_1000        = 318785;
      int nevts_qcd_1000_1400       = 1964088;
      int nevts_qcd_em_20_30        = 35040695;
      int nevts_qcd_em_30_80        = 33088888;
      int nevts_qcd_em_80_170       = 34542763;
      int nevts_qcd_em_170_250      = 31697066;
      int nevts_qcd_em_250_350      = 34611322;
      int nevts_qcd_em_350          = 34080562;
      int nevts_qcd_bc_20_30        = 1740229;
      int nevts_qcd_bc_30_80        = 2048152;
      int nevts_qcd_bc_80_170       = 1945525;
      int nevts_qcd_bc_170_250      = 1948112;
      int nevts_qcd_bc_250_350      = 2026521;
      int nevts_gamma_0_15          = 2000488;
      int nevts_gamma_15_30         = 1970745;
      int nevts_gamma_30_50         = 1993325;
      int nevts_gamma_50_80         = 1995062;
      int nevts_gamma_80_120        = 1992627;
      int nevts_gamma_120_170       = 2000043;
      int nevts_gamma_170_300       = 2000069;
      int nevts_gamma_300_470       = 2000130;
      int nevts_ww                  = 10000431;
      int nevts_wz                  = 10000283;
      int nevts_zz                  = 9799908;
      int nevts_tbar_tw             = 493460;
      int nevts_t_tw                = 497658;
      int nevts_wjetstolnu          = 57709905;

      double xsec = gi_xsec*gi_eff;
      double scale = 1.0;
      if( isMC ){

         if( channel.compare("DYJetsToLL") == 0 ) evtemp.weight *= xsec_dyjetstoll/nevts_dyjetstoll;
         else if( channel.compare("DYJetsToLL_M10To50") == 0 )
            evtemp.weight *= xsec_dyjetstoll_m10to50/nevts_dyjetstoll_m10to50;
         else if( channel.compare("TTJets") == 0 ) evtemp.weight *= xsec_ttjets/nevts_ttjets;
         else if( channel.compare("TTJets_Hadronic") == 0 )
            evtemp.weight *= xsec_ttjets_0lept/nevts_ttjets_0lept;
         else if( channel.compare("TTJets_SemiLept") == 0 )
            evtemp.weight *= xsec_ttjets_1lept/nevts_ttjets_1lept;
         else if( channel.compare("TTJets_FullLept") == 0 )
            evtemp.weight *= xsec_ttjets_2lept/nevts_ttjets_2lept;
         else if( channel.compare("QCD_15_30") == 0 )
            evtemp.weight *= scale*xsec/nevts_qcd_15_30;
         else if( channel.compare("QCD_30_50") == 0 )
            evtemp.weight *= scale*xsec/nevts_qcd_30_50;
         else if( channel.compare("QCD_50_80") == 0 )
            evtemp.weight *= scale*xsec/nevts_qcd_50_80;
         else if( channel.compare("QCD_80_120") == 0 )
            evtemp.weight *= scale*xsec/nevts_qcd_80_120;
         else if( channel.compare("QCD_120_170") == 0 )
            evtemp.weight *= scale*xsec/nevts_qcd_120_170;
         else if( channel.compare("QCD_170_300") == 0 )
            evtemp.weight *= scale*xsec/nevts_qcd_170_300;
         else if( channel.compare("QCD_300_470") == 0 )
            evtemp.weight *= scale*xsec/nevts_qcd_300_470;
         else if( channel.compare("QCD_470_600") == 0 )
            evtemp.weight *= scale*xsec/nevts_qcd_470_600;
         else if( channel.compare("QCD_600_800") == 0 )
            evtemp.weight *= scale*xsec/nevts_qcd_600_800;
         else if( channel.compare("QCD_800_1000") == 0 )
            evtemp.weight *= scale*xsec/nevts_qcd_800_1000;
         else if( channel.compare("QCD_1000_1400") == 0 )
            evtemp.weight *= scale*xsec/nevts_qcd_1000_1400;
         else if( channel.compare("QCD_EMEnriched_20_30") == 0 )
            evtemp.weight *= scale*xsec/nevts_qcd_em_20_30;
         else if( channel.compare("QCD_EMEnriched_30_80") == 0 )
            evtemp.weight *= scale*xsec/nevts_qcd_em_30_80;
         else if( channel.compare("QCD_EMEnriched_80_170") == 0 )
            evtemp.weight *= scale*xsec/nevts_qcd_em_80_170;
         else if( channel.compare("QCD_EMEnriched_170_250") == 0 )
            evtemp.weight *= scale*xsec/nevts_qcd_em_170_250;
         else if( channel.compare("QCD_EMEnriched_250_350") == 0 )
            evtemp.weight *= scale*xsec/nevts_qcd_em_250_350;
         else if( channel.compare("QCD_EMEnriched_350") == 0 )
            evtemp.weight *= scale*xsec/nevts_qcd_em_350;
         else if( channel.compare("QCD_BCtoE_20_30") == 0 )
            evtemp.weight *= scale*xsec/nevts_qcd_bc_20_30;
         else if( channel.compare("QCD_BCtoE_30_80") == 0 )
            evtemp.weight *= scale*xsec/nevts_qcd_bc_30_80;
         else if( channel.compare("QCD_BCtoE_80_170") == 0 )
            evtemp.weight *= scale*xsec/nevts_qcd_bc_80_170;
         else if( channel.compare("QCD_BCtoE_170_250") == 0 )
            evtemp.weight *= scale*xsec/nevts_qcd_bc_170_250;
         else if( channel.compare("QCD_BCtoE_250_350") == 0 )
            evtemp.weight *= scale*xsec/nevts_qcd_bc_250_350;
         else if( channel.compare("Gamma_0_15") == 0 )
            evtemp.weight *= scale*gi_xsec/nevts_gamma_0_15;
         else if( channel.compare("Gamma_15_30") == 0 )
            evtemp.weight *= scale*gi_xsec/nevts_gamma_15_30;
         else if( channel.compare("Gamma_30_50") == 0 )
            evtemp.weight *= scale*gi_xsec/nevts_gamma_30_50;
         else if( channel.compare("Gamma_50_80") == 0 )
            evtemp.weight *= scale*gi_xsec/nevts_gamma_50_80;
         else if( channel.compare("Gamma_80_120") == 0 )
            evtemp.weight *= scale*gi_xsec/nevts_gamma_80_120;
         else if( channel.compare("Gamma_120_170") == 0 )
            evtemp.weight *= scale*gi_xsec/nevts_gamma_120_170;
         else if( channel.compare("Gamma_170_300") == 0 )
            evtemp.weight *= scale*gi_xsec/nevts_gamma_170_300;
         else if( channel.compare("Gamma_300_470") == 0 )
            evtemp.weight *= scale*gi_xsec/nevts_gamma_300_470;

         else if( channel.compare("WW") == 0 ) evtemp.weight *= xsec_ww/nevts_ww;
         else if( channel.compare("WZ") == 0 ) evtemp.weight *= xsec_wz/nevts_wz;
         else if( channel.compare("ZZ") == 0 ) evtemp.weight *= xsec_zz/nevts_zz;
         else if( channel.compare("Tbar_tW") == 0 ) evtemp.weight *= xsec_tbar_tw/nevts_tbar_tw;
         else if( channel.compare("T_tW") == 0 ) evtemp.weight *= xsec_t_tw/nevts_t_tw;
         else if( channel.compare("WJetsToLNu") == 0 )
            evtemp.weight *= xsec_wjetstolnu/nevts_wjetstolnu;
         else cout << "No Xsection for channel " << channel << endl;

         evtemp.weight *= double(nevts_dyjetstoll)/xsec_dyjetstoll;
         evtemp.weight *= puMyWeight;
      }

      // vertices
      evtemp.nvertices = v_size;

      // old metsig
      evtemp.metsig2011 = metsig2011[0];

      // muons
      for( int i=0; i < mu_size; i++){
         evtemp.muon_pt.push_back( mu_pt[i] );
         evtemp.muon_phi.push_back( mu_phi[i] );
      }

      // electrons
      for( int i=0; i < elec_size; i++){
         evtemp.electron_pt.push_back( elec_pt[i] );
         evtemp.electron_phi.push_back( elec_phi[i] );
      }

      // jets
      for( int i=0; i < pfj_size; i++){

         double jet_ptL123_temp = (pfj_pt[i]*pfj_l1l2l3[i] > jetcorrpt)
            ? pfj_pt[i]*pfj_l1l2l3[i] : pfj_pt[i];
         double jet_ptT1_temp = (pfj_pt[i]*pfj_l1l2l3[i] > jetcorrpt)
            ? pfj_pt[i]*(pfj_l1l2l3[i] + 1 - pfj_l1[i]) : pfj_pt[i];

         if( jet_ptL123_temp > jetbinpt ){
            // clustered jets

            evtemp.jet_phi.push_back( pfj_phi[i] );
            evtemp.jet_eta.push_back( pfj_eta[i] );
            evtemp.jet_ptUncor.push_back( pfj_pt[i] );

            evtemp.jet_ptL123.push_back( jet_ptL123_temp );
            evtemp.jet_ptT1.push_back( jet_ptT1_temp );

         } else {
            // pseudojet with unclustered energy

            pjet_scalptL123_temp += jet_ptL123_temp;
            pjet_pxL123_temp += jet_ptL123_temp*cos(pfj_phi[i]);
            pjet_pyL123_temp += jet_ptL123_temp*sin(pfj_phi[i]);

            pjet_scalptT1_temp += jet_ptT1_temp;
            pjet_pxT1_temp += jet_ptT1_temp*cos(pfj_phi[i]);
            pjet_pyT1_temp += jet_ptT1_temp*sin(pfj_phi[i]);

            pjet_size_temp++;

            pjet_pxpx_temp += pow(jet_ptL123_temp*cos(pfj_phi[i]),2);
            pjet_pypy_temp += pow(jet_ptL123_temp*sin(pfj_phi[i]),2);
            pjet_pxpy_temp += pow(jet_ptL123_temp,2)*cos(pfj_phi[i])*sin(pfj_phi[i]);

         }

      } // pfj loop

      // pseudojet
      evtemp.pjet_phiL123 = atan2( pjet_pyL123_temp, pjet_pxL123_temp );
      evtemp.pjet_scalptL123 = pjet_scalptL123_temp;
      evtemp.pjet_vectptL123 = sqrt( pow(pjet_pxL123_temp,2) + pow(pjet_pyL123_temp,2) );

      evtemp.pjet_phiT1 = atan2( pjet_pyT1_temp, pjet_pxT1_temp );
      evtemp.pjet_scalptT1 = pjet_scalptT1_temp;
      evtemp.pjet_vectptT1 = sqrt( pow(pjet_pxT1_temp,2) + pow(pjet_pyT1_temp,2) );

      evtemp.pjet_size = pjet_size_temp;

      // met
      for(int i=0; i < 5; i++){
         evtemp.pfmet_px[i] = met_px[i];
         evtemp.pfmet_py[i] = met_py[i];
      }
      
      // met smearing
      evtemp.met_varx = 0.0;
      evtemp.met_vary = 0.0;
      evtemp.met_rho = 0.0;

      // fill event vector
      eventref_temp.push_back( evtemp );

   } // event loop

   if(do_resp_correction) {
	   cout << " -----> met response correction" << endl;
	   double xtemp[20]={0};
	   FindSignificance(xtemp,eventref_temp);
	   ResponseCorrection(eventref_temp);
   }

   return;
}

void Fitter::ResponseCorrection(vector<event>& eventvec) {

	vector<event>* eventref=&eventvec;
	TProfile* presp_qt = new TProfile ("presp_qt_temp",
			"Response = |<u_{#parallel}>|/q_{T} vs. q_{T};q_{T} (GeV);Response", 25, 0, 100);

	for( vector<event>::iterator ev = eventref->begin(); ev < eventref->end(); ev++ ){
		presp_qt->Fill( ev->qt, -(ev->ut_par)/(ev->qt), ev->weight );
	}

	TF1* func=new TF1("response curve fit","[0]+[1]*exp([2]*x)");
	func->SetParName(0,"Offset");
	func->SetParName(1,"Scale");
	func->SetParName(2,"Power");

	func->SetParameter(0,1);
	func->SetParameter(1,-0.4);;

	presp_qt->Fit(func,"NOQ");


	for( vector<event>::iterator ev = eventvec.begin(); ev < eventvec.end(); ev++ ){

		double pt_mult=abs( 1/func->Eval(ev->qt) );
		ev->resp_correction=pt_mult;

		// pseudojet
		ev->pjet_vectptL123*=pt_mult;
		ev->pjet_scalptL123*=pt_mult;
		ev->pjet_vectptT1*=pt_mult;
		ev->pjet_scalptT1*=pt_mult;
	}
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
   const double *xmin = gMinuit->X();
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

   // random number engine for MET-smearing
   ROOT::Math::Random<ROOT::Math::GSLRngMT> GSLr;

   int count = 0;
   // event loop
   for( vector<event>::iterator ev = eventref_temp.begin(); ev < eventref_temp.end(); ev++){

      double met_x=0;
      double met_y=0;
      double cov_xx=0;
      double cov_xy=0;
      double cov_yy=0;
      double cov_xx_highpt=0;
      double cov_xx_pjet=0;

      // clustered jets
      for(int i=0; i < int(ev->jet_ptUncor.size()); i++){

         float feta = fabs(ev->jet_eta[i]);
         double c = cos(ev->jet_phi[i]);
         double s = sin(ev->jet_phi[i]);

         met_x -= c*(ev->jet_ptT1[i]);
         met_y -= s*(ev->jet_ptT1[i]);

         double dpt=0;
         double dph=0;

         // resolutions for two jet categories
         if( true || ev->jet_ptL123[i] > jetbinpt ){ //dummy true

            int index=-1;
            if(feta<0.5) index=0;
            else if(feta<1.1) index=1;
            else if(feta<1.7) index=2;
            else if(feta<2.3) index=3;
            else{
               index=4;
            }

            // CMS 2010 Resolutions -- parameterized by L123 corrected pt
            dpt = x[index] * (ev->jet_ptT1[i]) * dpt_(ev->jet_ptL123[i], ev->jet_eta[i]);
            dph =            (ev->jet_ptT1[i]) * dph_(ev->jet_ptL123[i], ev->jet_eta[i]);

         }else{
            cout << "ERROR: JET PT OUT OF RANGE" << endl;
         }

         double dtt = dpt*dpt;
         double dff = dph*dph;
         cov_xx += dtt*c*c + dff*s*s;
         cov_xy += (dtt-dff)*c*s;
         cov_yy += dff*c*c + dtt*s*s;

         cov_xx_highpt += dtt*c*c + dff*s*s;
      }

      // muons -- assume zero resolutions
      for(int i=0; i < int(ev->muon_pt.size()); i++){
         met_x -= cos(ev->muon_phi[i])*(ev->muon_pt[i]);
         met_y -= sin(ev->muon_phi[i])*(ev->muon_pt[i]);
      }

      // electrons -- assume zero resolutions
      for(int i=0; i < int(ev->electron_pt.size()); i++){
         met_x -= cos(ev->electron_phi[i])*(ev->electron_pt[i]);
         met_y -= sin(ev->electron_phi[i])*(ev->electron_pt[i]);
      }

      // unclustered energy -- parameterize by scalar sum of ET
      double c = cos(ev->pjet_phiT1);
      double s = sin(ev->pjet_phiT1);

      met_x -= c*(ev->pjet_vectptT1);
      met_y -= s*(ev->pjet_vectptT1);

      double ctt = x[5]*x[5] + x[6]*x[6]*(ev->pjet_scalptL123);

      cov_xx += ctt;
      cov_yy += ctt;

      cov_xx_pjet += ctt;

      double det = cov_xx*cov_yy - cov_xy*cov_xy;

      double ncov_xx = cov_yy / det;
      double ncov_xy = -cov_xy / det;
      double ncov_yy = cov_xx / det;

      if( met_type != -1 ){
         met_x = ev->pfmet_px[met_type];
         met_y = ev->pfmet_py[met_type];
      }

      // smear MC MET
      double smear_x = 0;
      double smear_y = 0;
      double sigma_x = 0;
      double sigma_y = 0;
      if( ev->met_varx >= 0.0 and ev->met_vary >= 0.0 and fabs(ev->met_rho) <= 1.0 ){

         sigma_x = sqrt(ev->met_varx);
         sigma_y = sqrt(ev->met_vary);

         GSLr.Gaussian2D( sigma_x, sigma_y, ev->met_rho, smear_x, smear_y );
         met_x += smear_x;
         met_y += smear_y;

      }
      count++;

      double sig = met_x*met_x*ncov_xx + 2*met_x*met_y*ncov_xy + met_y*met_y*ncov_yy;

      ev->met = sqrt( met_x*met_x + met_y*met_y );
      ev->sig = sig;
      ev->det = det;

      ev->cov_xx = cov_xx;
      ev->cov_xy = cov_xy;
      ev->cov_yy = cov_yy;

      ev->cov_xx_highpt = cov_xx_highpt;
      ev->cov_xx_pjet = cov_xx_pjet;

      // fill qt, ut
      double qt_x=0, qt_y=0;
      for(int i=0; i < int(ev->muon_pt.size()); i++){
         qt_x += ev->muon_pt[i]*cos(ev->muon_phi[i]);
         qt_y += ev->muon_pt[i]*sin(ev->muon_phi[i]);
      }
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
/*
void Fitter::MatchMCjets(vector<event>& eventref_temp){
   cout << "---> MatchMCjets" << endl;

   for( vector<event>::iterator ev = eventref_temp.begin(); ev < eventref_temp.end(); ev++){

      // loop through reco jets
      for(int ireco=0; ireco < int(ev->jet_ptUncor.size()); ireco++){

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

      if( ev->jet_ptUncor.size() != ev->jet_matchIndex.size() ){
         cout << "ERROR: VECTOR SIZE MISMATCH" << endl;
      }

   } // event loop

}
*/

void Fitter::FillHists(vector<event>& eventref, string stackmode){
   if( eventref.size() == 0 ) return;

   vector<event>::iterator iter_begin = eventref.begin();
   vector<event>::iterator iter_end = eventref.end();

   map<string, TH1*> hists_;
   map<string, TH2*> profs_;

   if( iter_begin->channel.compare("Data") == 0 ){
      hists_ = histsData_;
      profs_ = profsData_;
   } else {
      profs_ = profsMC_;
      if( stackmode.compare("Zmumu") == 0 ){

         if( iter_begin->channel.compare("DYJetsToLL") == 0 ) hists_ = histsMC_signal_;
         else if( iter_begin->channel.compare("TTJets") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->channel.compare("Tbar_tW") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->channel.compare("T_tW") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->channel.compare("WJetsToLNu") == 0 ) hists_ = histsMC_EWK_;
         else if( iter_begin->channel.compare("WW") == 0 ) hists_ = histsMC_EWK_;
         else if( iter_begin->channel.compare("WZ") == 0 ) hists_ = histsMC_EWK_;
         else if( iter_begin->channel.compare("ZZ") == 0 ) hists_ = histsMC_EWK_;
         else cout << "Histogram fill error, channel " << iter_begin->channel << endl;

      }
      if( stackmode.compare("Wenu") == 0 or stackmode.compare("Wenu_loose") == 0 ){

         if( iter_begin->channel.compare("WJetsToLNu") == 0 ) hists_ = histsMC_signal_;
         else if( iter_begin->channel.compare("DYJetsToLL") == 0 ) hists_ = histsMC_DY_;
         else if( iter_begin->channel.compare("DYJetsToLL_M10To50") == 0 ) hists_ = histsMC_DY_;
         else if( iter_begin->channel.compare("TTJets") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->channel.compare("Tbar_tW") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->channel.compare("T_tW") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->channel.compare("QCD_EMEnriched_20_30") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_EMEnriched_30_80") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_EMEnriched_80_170") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_EMEnriched_170_250") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_EMEnriched_250_350") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_EMEnriched_350") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_BCtoE_20_30") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_BCtoE_30_80") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_BCtoE_80_170") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_BCtoE_170_250") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_BCtoE_250_350") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("Gamma_0_15") == 0 ) hists_ = histsMC_gamma_;
         else if( iter_begin->channel.compare("Gamma_15_30") == 0 ) hists_ = histsMC_gamma_;
         else if( iter_begin->channel.compare("Gamma_30_50") == 0 ) hists_ = histsMC_gamma_;
         else if( iter_begin->channel.compare("Gamma_50_80") == 0 )  hists_ = histsMC_gamma_;
         else if( iter_begin->channel.compare("Gamma_80_120") == 0 ) hists_ = histsMC_gamma_;
         else if( iter_begin->channel.compare("Gamma_120_170") == 0 ) hists_ = histsMC_gamma_;  
         else if( iter_begin->channel.compare("Gamma_170_300") == 0 ) hists_ = histsMC_gamma_;  
         else if( iter_begin->channel.compare("Gamma_300_470") == 0 ) hists_ = histsMC_gamma_;  
         else if( iter_begin->channel.compare("WW") == 0 ) hists_ = histsMC_EWK_;
         else if( iter_begin->channel.compare("WZ") == 0 ) hists_ = histsMC_EWK_;
         else if( iter_begin->channel.compare("ZZ") == 0 ) hists_ = histsMC_EWK_;
         else cout << "Histogram fill error, channel " << iter_begin->channel << endl;

      }
      if( stackmode.compare("Dijet") == 0 ){

         hists_ = histsMC_signal_;

      }
      if( stackmode.compare("Ttbar0lept") == 0){

         if( iter_begin->channel.compare("TTJets_Hadronic") == 0 ) hists_ = histsMC_signal_;
         else if( iter_begin->channel.compare("TTJets_FullLept") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->channel.compare("TTJets_SemiLept") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->channel.compare("DYJetsToLL") == 0 ) hists_ = histsMC_DY_;
         else if( iter_begin->channel.compare("Tbar_tW") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->channel.compare("T_tW") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->channel.compare("WJetsToLNu") == 0 ) hists_ = histsMC_EWK_;
         else if( iter_begin->channel.compare("QCD_15_30") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_30_50") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_50_80") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_80_120") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_120_170") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_170_300") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_300_470") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_470_600") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_600_800") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_800_1000") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_1000_1400") == 0 ) hists_ = histsMC_QCD_;
         else cout << "Histogram fill error, channel " << iter_begin->channel << endl;

      }
      if( stackmode.compare("Ttbar1lept") == 0 ){

         if( iter_begin->channel.compare("TTJets_SemiLept") == 0 ) hists_ = histsMC_signal_;
         else if( iter_begin->channel.compare("TTJets_FullLept") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->channel.compare("TTJets_Hadronic") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->channel.compare("WJetsToLNu") == 0 ) hists_ = histsMC_signal_;
         else if( iter_begin->channel.compare("DYJetsToLL") == 0 ) hists_ = histsMC_DY_;
         else if( iter_begin->channel.compare("DYJetsToLL_M10To50") == 0 ) hists_ = histsMC_DY_;
         else if( iter_begin->channel.compare("Tbar_tW") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->channel.compare("T_tW") == 0 ) hists_ = histsMC_top_;
         else if( iter_begin->channel.compare("QCD_EMEnriched_20_30") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_EMEnriched_30_80") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_EMEnriched_80_170") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_EMEnriched_170_250") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_EMEnriched_250_350") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_EMEnriched_350") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_BCtoE_20_30") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_BCtoE_30_80") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_BCtoE_80_170") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_BCtoE_170_250") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("QCD_BCtoE_250_350") == 0 ) hists_ = histsMC_QCD_;
         else if( iter_begin->channel.compare("Gamma_0_15") == 0 ) hists_ = histsMC_gamma_;
         else if( iter_begin->channel.compare("Gamma_15_30") == 0 ) hists_ = histsMC_gamma_;
         else if( iter_begin->channel.compare("Gamma_30_50") == 0 ) hists_ = histsMC_gamma_;
         else if( iter_begin->channel.compare("Gamma_50_80") == 0 )  hists_ = histsMC_gamma_;
         else if( iter_begin->channel.compare("Gamma_80_120") == 0 ) hists_ = histsMC_gamma_;
         else if( iter_begin->channel.compare("Gamma_120_170") == 0 ) hists_ = histsMC_gamma_;  
         else if( iter_begin->channel.compare("Gamma_170_300") == 0 ) hists_ = histsMC_gamma_;  
         else if( iter_begin->channel.compare("Gamma_300_470") == 0 ) hists_ = histsMC_gamma_;  
         else if( iter_begin->channel.compare("WW") == 0 ) hists_ = histsMC_EWK_;
         else if( iter_begin->channel.compare("WZ") == 0 ) hists_ = histsMC_EWK_;
         else if( iter_begin->channel.compare("ZZ") == 0 ) hists_ = histsMC_EWK_;
         else cout << "Histogram fill error, channel " << iter_begin->channel << endl;

      }

   }

   for( vector<event>::iterator ev = iter_begin; ev < iter_end; ev++ ){

      // muons
      for( int j=0; j < int(ev->muon_pt.size()); j++){
         hists_["muon_pt"]->Fill( ev->muon_pt[j] , ev->weight );
      }

      // jets
      hists_["njets"]->Fill( ev->jet_ptL123.size(), ev->weight );
      for( int j=0; j < int(ev->jet_ptL123.size()); j++){
         hists_["jet_pt"]->Fill( ev->jet_ptL123[j], ev->weight );
         if( ev->jet_ptL123[j] > 30 and ev->jet_ptUncor[j] != ev->jet_ptL123[j] ){
            hists_["jet_eta"]->Fill( ev->jet_eta[j], ev->weight );
         }
      }
      if( ev->jet_ptL123.size() > 0 ){
         hists_["jet1_pt"]->Fill( ev->jet_ptL123[0] , ev->weight );
      }

      // pseudojet
      hists_["pjet_scalptL123"]->Fill( ev->pjet_scalptL123 , ev->weight );
      hists_["pjet_vectptL123"]->Fill( ev->pjet_vectptL123 , ev->weight );
      hists_["pjet_scalptT1"]->Fill( ev->pjet_scalptT1 , ev->weight );
      hists_["pjet_vectptT1"]->Fill( ev->pjet_vectptT1 , ev->weight );
      hists_["pjet_size"]->Fill( ev->pjet_size , ev->weight );
      hists_["pjet_phi"]->Fill( ev->pjet_phiL123, ev->weight );

      // other observables
      hists_["nvert"]->Fill( ev->nvertices , ev->weight );

      hists_["qt"]->Fill( ev->qt, ev->weight );
      hists_["ut_par"]->Fill( fabs(ev->ut_par), ev->weight );

      hists_["cov_xx"]->Fill( ev->cov_xx, ev->weight );
      hists_["cov_xy"]->Fill( ev->cov_xy, ev->weight );
      hists_["cov_yy"]->Fill( ev->cov_yy, ev->weight );
      hists_["met"]->Fill( ev->met, ev->weight );
      hists_["met_200"]->Fill( ev->met, ev->weight );
      hists_["sig"]->Fill( ev->sig, ev->weight );
      hists_["sig_100"]->Fill( ev->sig, ev->weight );
      hists_["sig_15"]->Fill( ev->sig, ev->weight );
      hists_["sig_old"]->Fill( ev->metsig2011, ev->weight );
      hists_["det"]->Fill( ev->det, ev->weight );
      hists_["pchi2"]->Fill( TMath::Prob(ev->sig,2), ev->weight );
      hists_["pchi2_old"]->Fill( TMath::Prob(ev->metsig2011,2), ev->weight );
      hists_["logpchi2"]->Fill( TMath::Log(TMath::Prob(ev->sig,2)), ev->weight );

      hists_["cov_xx_highpt"]->Fill( ev->cov_xx_highpt, ev->weight );
      hists_["cov_xx_pjet"]->Fill( ev->cov_xx_pjet, ev->weight );
      hists_["cov_xx_ratio"]->Fill( ev->cov_xx_highpt/ev->cov_xx, ev->weight );

      hists_["met_varx"]->Fill( ev->met_varx, ev->weight );
      hists_["met_vary"]->Fill( ev->met_vary, ev->weight );
      hists_["met_rho"]->Fill( ev->met_rho, ev->weight );

      // profiles
      profs_["psig_nvert"]->Fill( ev->nvertices, ev->sig, ev->weight );
      profs_["psig_qt"]->Fill( ev->qt, ev->sig, ev->weight );
      profs_["presp_qt"]->Fill( ev->qt, -(ev->ut_par)/(ev->qt), ev->weight );
      profs_["pMET_nvert"]->Fill( ev->nvertices, ev->met, ev->weight );
      profs_["pjet_scalptL123_nvert"]->Fill( ev->nvertices, ev->pjet_scalptL123, ev->weight );

      profs_["njets_nvert"]->Fill( ev->nvertices, ev->jet_ptL123.size(), ev->weight );
      for( int j=0; j < int(ev->jet_ptL123.size()); j++){
         profs_["jet_pt_nvert"]->Fill( ev->nvertices, ev->jet_ptL123[j], ev->weight );
      }

      profs_["sig_met"]->Fill( ev->met, ev->sig, ev->weight );

   }

}

void Fitter::PrintHists( const char* filename, string stackmode ){

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
   TH1D *histMC_EWK_temp = (TH1D*)histsMC_EWK_["met"]->Clone("histMC_EWK_temp");
   TH1D *histMC_QCD_temp = (TH1D*)histsMC_QCD_["met"]->Clone("histMC_QCD_temp");
   TH1D *histMC_gamma_temp = (TH1D*)histsMC_gamma_["met"]->Clone("histMC_gamma_temp");
   TH1D *histMC_DY_temp = (TH1D*)histsMC_DY_["met"]->Clone("histMC_DY_temp");

   // rescale QCD & gamma+jets numerically (approximate)
   double chi2 = -1;
   double histnorm = histData_temp->Integral("width") / (histMC_signal_temp->Integral("width")+histMC_top_temp->Integral("width")+histMC_EWK_temp->Integral("width")+histMC_QCD_temp->Integral("width")+histMC_gamma_temp->Integral("width")+histMC_DY_temp->Integral("width"));
   double scaleQCD = 1;
   if( stackmode.compare("Wenu") == 0 or stackmode.compare("Wenu_loose") == 0 
         or stackmode.compare("Ttbar0lept") == 0 ){ 
      for(double s = 0; s < 5; s += 0.01){
         TH1D *histMC_temp = new TH1D( "histMC_temp", "histMC_temp",
               histData_temp->GetNbinsX(), histData_temp->GetBinLowEdge(1),
               histData_temp->GetBinLowEdge(histData_temp->GetNbinsX())
               + histData_temp->GetBinWidth(histData_temp->GetNbinsX()) );

         histMC_temp->Add( histMC_signal_temp );
         histMC_temp->Add( histMC_top_temp );
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
      TH1D *histMC_EWK = (TH1D*)histsMC_EWK_[hname];
      TH1D *histMC_QCD = (TH1D*)histsMC_QCD_[hname];
      TH1D *histMC_gamma = (TH1D*)histsMC_gamma_[hname];
      TH1D *histMC_DY = (TH1D*)histsMC_DY_[hname];

      histMC_QCD->Scale( scaleQCD );
      histMC_gamma->Scale( scaleQCD );

      // get total MC histogram
      histMC->Add( histMC_signal );
      histMC->Add( histMC_top );
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
         histMC_top->SetLineColor(1);
         histMC_top->SetFillColor(kYellow-9);
         histMC_DY->SetLineColor(1);
         histMC_DY->SetFillColor(kRed-10);
         histMC_EWK->SetLineColor(1);
         histMC_EWK->SetFillColor(kCyan-10);

         histMC_top->Add( histMC_DY );
         histMC_top->Add( histMC_EWK );
         histMC_DY->Add( histMC_EWK );

         histMC_top->Draw("same HIST");
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
         // nothing happens here
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
      }
      if( stackmode.compare("Ttbar1lept") == 0 ){
         leg->AddEntry(histMC, "t#overline{t} Semi-Leptonic", "f");
         leg->AddEntry(histMC_top, "Other Top", "f");
         leg->AddEntry(histMC_DY, "DY", "f");
         leg->AddEntry(histMC_EWK, "EWK", "f");
      }
      if( stackmode.compare("Ttbar0lept") == 0 ){
         leg->AddEntry(histMC, "t#overline{t} Hadronic", "f");
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

}
