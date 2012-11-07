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

   fitter.PlotsDataMC( eventvec_MC, eventvec_data, "plots/DataMC.root" );

   return 0;
}
