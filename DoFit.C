#include <iostream>
#include "METSigFit.h"

int main(){

   std::cout << "Hello World 1" << std::endl;
   
   Fitter fitter;
   fitter.ReadNtuple("Zmumu_MC_20121005.root",false);

   std::cout << "Hello World 2" << std::endl;

   return 0;
}
