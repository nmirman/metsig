#include <iostream>
#include "METSigFit.h"

#include "TH1.h"
#include "TCanvas.h"

int main(){

   Fitter fitter;
   fitter.ReadNtuple("/eos/uscms/store/user/nmirman/Zmumu/Zmumu_ntuple_20121005.root",true);
   fitter.RunMinimizer( fitter.eventvec_MC );

   TH1D *hjet_pt = new TH1D("hjet_pt","Jet p_{T}",100,0,200);
   TH1D *hsig = new TH1D("hsig","Significance",100,0,20);

   for( std::vector<Fitter::event>::iterator ev = fitter.eventvec_MC.begin(); 
         ev < fitter.eventvec_MC.end(); ev++){
      hsig->Fill( ev->sig );
      for(int i=0; i < ev->jet_pt.size(); i++){
         hjet_pt->Fill( ev->jet_pt[i] );
      }
   }

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
