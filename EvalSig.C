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

   int numevents = -1;
   bool stackMC = true;
   bool do_resp_correction = false;

   // declarations
   Fitter fitter;
   vector<event> eventvec_MC;
   vector<event> eventvec_data;

   // fill eventvecs
   
   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130501/DYJetsToLL.root",
		   eventvec_MC, numevents, true, "DYJetsToLL", do_resp_correction);
   if( stackMC ){
      fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130501/QCD.root",
            eventvec_MC, numevents, true, "QCD", do_resp_correction);
      fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130501/TTJets.root",
            eventvec_MC, numevents, true, "TTJets", do_resp_correction);
      fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130501/Tbar_tW-channel.root",
            eventvec_MC, numevents, true, "Tbar_tW", do_resp_correction);
      fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130501/T_tW-channel.root",
            eventvec_MC, numevents, true, "T_tW", do_resp_correction);
      fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130501/WJetsToLNu.root",
            eventvec_MC, numevents, true, "WJetsToLNu", do_resp_correction);
      fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130501/WW.root",
            eventvec_MC, numevents, true, "WW", do_resp_correction);
      fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130501/WZ.root",
            eventvec_MC, numevents, true, "WZ", do_resp_correction);
      fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130501/ZZ.root",
            eventvec_MC, numevents, true, "ZZ", do_resp_correction);
   }
   fitter.MatchMCjets( eventvec_MC );

   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130427/Data.root",
         eventvec_data, numevents, false, "Data", false);

   cout << "\n  MC EVENTS: " << eventvec_MC.size() << endl;
   cout << "DATA EVENTS: " << eventvec_data.size() << endl;

   double parMC [] =   {1.12660,1.09322,1.10951,1.17178,1.12164,0.0,0.585145};
   double parData [] = {1.39669,1.32037,1.32047,1.38161,1.51508,0.0,0.639158};

   vector<event> eventvec_sigmaMC = eventvec_MC;
   vector<event> eventvec_sigmaData = eventvec_MC;

   // reweight MC for pseudojet scalar pt
   fitter.PJetReweight( eventvec_MC, eventvec_data );
   fitter.PJetReweight( eventvec_sigmaData, eventvec_data );

   fitter.FindSignificance(parMC, eventvec_sigmaMC);
   fitter.FindSignificance(parData, eventvec_sigmaData);

   for( int i=0; i < int(eventvec_MC.size()); i++ ){
      eventvec_MC[i].met_varx = eventvec_sigmaData[i].cov_xx - eventvec_sigmaMC[i].cov_xx;
      eventvec_MC[i].met_vary = eventvec_sigmaData[i].cov_yy - eventvec_sigmaMC[i].cov_yy;
      eventvec_MC[i].met_rho = (eventvec_sigmaData[i].cov_xy - eventvec_sigmaMC[i].cov_xy)
         / sqrt(eventvec_MC[i].met_varx * eventvec_MC[i].met_vary);
   }

   fitter.FindSignificance(parData, eventvec_MC);
   fitter.FindSignificance(parData, eventvec_data);

   fitter.PlotsDataMC( eventvec_data, eventvec_MC, "results/plotsDataMC.root", stackMC, "Zmumu");

   return 0;
}
