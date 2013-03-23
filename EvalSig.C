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

   double parMC [] =   {1.14210,1.06446,1.10291,1.08131,1.14072,0.000615132,0.606726,-0.905279};
   double parData [] = {1.30361,1.22821,1.24193,1.25893,1.43332,-0.00024461,0.695474,-2.88705};

   vector<event> eventvec_sigmaMC = eventvec_MC;
   vector<event> eventvec_sigmaData = eventvec_MC;

   // reweight MC for pseudojet scalar pt
   double pjet_weights [] = {0.0,1.17948,1.12869,1.10327,1.09697,1.08439,1.07407,1.06587,1.05725,
      1.04646,1.0377,1.03198,1.02595,1.02087,1.0134,1.01204,1.0073,1.0009,0.997246,0.994572,0.984984,
      0.98568,0.989608,0.973945,0.999722,1.00649,0.946263,0.912974,1.00385,0.805591};
   //fitter.PJetReweight( eventvec_MC, pjet_weights );
   //fitter.PJetReweight( eventvec_sigmaData, pjet_weights );

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

   fitter.PlotsDataMC( eventvec_data, eventvec_MC, "results/plotsDataMC.root", stackMC );

   return 0;
}
