#include <iostream>
#include "METSigFit.h"

#include "TH1.h"
#include "TCanvas.h"

int main(){

   const double parameter [12] = { 1.5, 1.5, 1.5, 1.5, 5.0,
                                   1.0, 1.0, 1.0,
                                   4.0, 0.5, 4.0, 0.5 };

   Fitter fitter;
   fitter.ReadNtuple("/eos/uscms/store/user/nmirman/Zmumu/Zmumu_ntuple_20121005.root",true);
   fitter.RunMinimizer( fitter.eventvec_MC );
   //fitter.FindSignificance( parameter, fitter.eventvec_MC );

   TH1D *hjet_pt = new TH1D("hjet_pt","Jet p_{T}",100,0,200);
   TH1D *hsig = new TH1D("hsig","Significance",100,0,20);

   double m2ll = 0;
   for( std::vector<Fitter::event>::iterator ev = fitter.eventvec_MC.begin(); 
         ev < fitter.eventvec_MC.end(); ev++){
      hsig->Fill( ev->sig );
      for(int i=0; i < ev->jet_pt.size(); i++){
         hjet_pt->Fill( ev->jet_pt[i] );
      }
      m2ll += ev->sig + log(ev->det);
   }
   std::cout << "M2LL = " << m2ll << std::endl;

   TCanvas *canvas = new TCanvas("canvas","canvas",1000,500);
   canvas->Divide(2,1);
   canvas->cd(1);
   hsig->Draw();
   canvas->cd(2);
   hjet_pt->Draw();
   canvas->cd();
   canvas->Print("plots.root");

   return 0;
}
