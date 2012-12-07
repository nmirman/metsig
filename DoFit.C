#include "METSigFit.h"

#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include <iostream>
#include <unistd.h>
using namespace std;

int main(int argc, char* argv[]){

   // tree with fit results
   double jetfitLOW=0, jetfitHIGH=0, psig_nvert_corr=0, psig_qt_corr=0,
          pchi2slope_left=0, pchi2slope_right=0,
          a1=0, a2=0, a3=0, a4=0, a5=0, k0=0, k1=0, k2=0, N1=0, S1=0;
   int fitStatus=-1;

   TTree *tree = new TTree("FitResults", "FitResults");
   tree->Branch("jetfitLOW", &jetfitLOW);
   tree->Branch("jetfitHIGH", &jetfitHIGH);
   tree->Branch("fitStatus", &fitStatus);
   tree->Branch("a1", &a1);
   tree->Branch("a2", &a2);
   tree->Branch("a3", &a3);
   tree->Branch("a4", &a4);
   tree->Branch("a5", &a5);
   tree->Branch("k0", &k0);
   tree->Branch("k1", &k1);
   tree->Branch("k2", &k2);
   tree->Branch("N1", &N1);
   tree->Branch("S1", &S1);
   tree->Branch("psig_vert_corr", &psig_nvert_corr);
   tree->Branch("psig_qt_corr", &psig_qt_corr);
   tree->Branch("pchi2slope_left", &pchi2slope_left);
   tree->Branch("pchi2slope_right", &pchi2slope_right);

   // declarations
   Fitter fitter;
   vector<event> eventvec_MC;
   vector<event> eventvec_data;


   // option flags
   char c;
   while( (c = getopt(argc, argv, "n:")) != -1 ) {
      switch(c)
      {
         case 'n' :
            fitter.jetfitLOW = atoi(optarg);
            break;

         default :
            continue;
      }
   }
   if( fitter.jetfitLOW > fitter.jetfitHIGH ) return 0;

   // fill eventvecs
   int numevents = 100000;
   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Zmumu/"
         "Zmumu_MC_DYJettoLL_TuneZ2_M-50_7TeV_madgraph_tauola_20121116.root",
         eventvec_MC, numevents, true);
   fitter.MatchMCjets( eventvec_MC );

   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Zmumu/"
         "Zmumu_data_DoubleMu_Run2011A_08Nov2011_v1_20121116.root",
         eventvec_data, numevents/2, false);
   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Zmumu/"
         "Zmumu_data_DoubleMu_Run2011B_19Nov2011_v1_20121116.root",
         eventvec_data, numevents/2, false);
   
   // pile-up reweighting
   fitter.GetPUWeights( eventvec_data, eventvec_MC );

   cout << "\n  MC EVENTS: " << eventvec_MC.size() << endl;
   cout << "DATA EVENTS: " << eventvec_data.size() << endl;

   // minimize
   cout << "\n ############################ " << endl;
   cout << " ###########  MC  ########### " << endl;
   cout << " ############################ \n" << endl;
   fitter.RunMinimizer( eventvec_MC );
   cout << "\n ############################ " << endl;
   cout << " ########### Data ########### " << endl;
   cout << " ############################ \n" << endl;
   fitter.RunMinimizer( eventvec_data );

   fitter.PlotsDataMC( eventvec_data, eventvec_MC, "results/plotsDataMC.root" );

   // fill tree with fit results
   jetfitLOW = fitter.jetfitLOW;
   jetfitHIGH = fitter.jetfitHIGH;
   fitStatus = fitter.gMinuit->Status();
   const double *par = fitter.gMinuit->X();
   a1 = par[0];
   a2 = par[1];
   a3 = par[2];
   a4 = par[3];
   a5 = par[4];
   k0 = par[5];
   k1 = par[6];
   k2 = par[7];
   N1 = par[8];
   S1 = par[9];
   psig_nvert_corr = fitter.psig_nvert_corr;
   psig_qt_corr = fitter.psig_qt_corr;
   pchi2slope_left = fitter.pchi2slope_left;
   pchi2slope_right = fitter.pchi2slope_right;

   tree->Fill();
   
   // set up output file path
   std::string pathstr;
   char* path = std::getenv("WORKING_DIR");
   if (path==NULL) {
      pathstr = "./results";
   }else {
      pathstr = path;
   }
   string outfilename = "/fitresults.root";

   TFile *file = new TFile((pathstr+outfilename).c_str(), "RECREATE");
   file->cd();
   tree->Write();
   file->Write();
   file->Close();

   return 0;
}
