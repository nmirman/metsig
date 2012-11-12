#include "METSigFit.h"

#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include <iostream>
using namespace std;

int main(){

   // declarations
   Fitter fitter;
   vector<event> eventvec_MC;
   vector<event> eventvec_data;

   // fill eventvecs
   int numevents = 100000;
   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Zmumu/"
         "Zmumu_MC_DYJettoLL_TuneZ2_M-50_7TeV_madgraph_tauola_20121107.root",
         eventvec_MC, numevents, true);
   fitter.MatchMCjets( eventvec_MC );

   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Zmumu/"
         "Zmumu_data_DoubleMu_Run2011A_08Nov2011_v1_20121104.root",
         eventvec_data, numevents/2, false);
   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Zmumu/"
         "Zmumu_data_DoubleMu_Run2011B_19Nov2011_v1_20121104.root",
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

   return 0;
}
