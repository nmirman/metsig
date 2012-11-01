#include "METSigFit.h"

int main(){

   Fitter fitter;
   std::vector<event> eventvec_MC;

   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Zmumu/Zmumu_ntuple_20121005.root",
         eventvec_MC, true);
   fitter.MatchMCjets( eventvec_MC );
   fitter.RunMinimizer( eventvec_MC );

   return 0;
}
