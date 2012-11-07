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
   std::vector<event> eventvec_MC;
   std::vector<event> eventvec_data;

   // fill eventvecs
   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Zmumu/"
         "Zmumu_MC_DYJettoLL_TuneZ2_M-50_7TeV_madgraph_tauola_20121104.root",
         eventvec_MC, 100000, true);
   fitter.MatchMCjets( eventvec_MC );

   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Zmumu/"
         "Zmumu_data_DoubleMu_Run2011A_08Nov2011_v1_20121104.root",
         eventvec_data, 50000, false);
   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Zmumu/"
         "Zmumu_data_DoubleMu_Run2011B_19Nov2011_v1_20121104.root",
         eventvec_data, 50000, false);

   std::cout << "MC EVENTS:   " << eventvec_MC.size() << std::endl;
   std::cout << "DATA EVENTS: " << eventvec_data.size() << std::endl;

   // minimize
   std::cout << " ########### MC ########### " << std::endl;
   fitter.RunMinimizer( eventvec_MC );
   std::cout << " ########### Data ########### " << std::endl;
   fitter.RunMinimizer( eventvec_data );

   fitter.PlotsDataMC( eventvec_MC, eventvec_data, "results/plotsDataMC.root" );

   return 0;
}
