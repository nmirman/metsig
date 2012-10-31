#include "METSigFit.h"

int main(){

   Fitter fitter;
   fitter.ReadNtuple("/eos/uscms/store/user/nmirman/Zmumu/Zmumu_ntuple_20121005.root",true);
   fitter.RunMinimizer( fitter.eventvec_MC );

   return 0;
}
