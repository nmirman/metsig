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
   
   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130321/DYJetsToLL.root",
		   eventvec_MC, numevents, true, "DYJetsToLL", do_resp_correction);
   if( stackMC ){
      fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130321/QCD.root",
            eventvec_MC, numevents, true, "QCD", do_resp_correction);
      fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130321/TTJets.root",
            eventvec_MC, numevents, true, "TTJets", do_resp_correction);
      fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130321/Tbar_tW-channel.root",
            eventvec_MC, numevents, true, "Tbar_tW", do_resp_correction);
      fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130321/T_tW-channel.root",
            eventvec_MC, numevents, true, "T_tW", do_resp_correction);
      fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130321/WJetsToLNu.root",
            eventvec_MC, numevents, true, "WJetsToLNu", do_resp_correction);
      fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130321/WW.root",
            eventvec_MC, numevents, true, "WW", do_resp_correction);
      fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130321/WZ.root",
            eventvec_MC, numevents, true, "WZ", do_resp_correction);
      fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130321/ZZ.root",
            eventvec_MC, numevents, true, "ZZ", do_resp_correction);
   }
   fitter.MatchMCjets( eventvec_MC );

   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130321/Data.root",
         eventvec_data, numevents, false, "Data", false);

   cout << "\n  MC EVENTS: " << eventvec_MC.size() << endl;
   cout << "DATA EVENTS: " << eventvec_data.size() << endl;

   double parMC [] =   {1.14271,1.05988,1.10172,1.08334,1.13369,0.000138076,0.606791,-0.907264};
   double parData [] = {1.30118,1.22218,1.24115,1.24444,1.43998,1.89428e-05,0.695186,-2.90158};

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
