#include <iostream>
#include "METSigFit.h"

int main(){

   std::cout << "Hello World 1" << std::endl;
   
   Fitter fitter;
   fitter.ReadNtuple("/eos/uscms/store/user/nmirman/Zmumu/Zmumu_ntuple_20121005.root",true);

   std::cout << "Hello World 2" << std::endl;

   return 0;
}
