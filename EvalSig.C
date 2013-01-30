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

   // declarations
   Fitter fitter;
   vector<event> eventvec_MC;
   vector<event> eventvec_data;

   // fill eventvecs
   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Zmumu/"
         "Zmumu_MC_DYJettoLL_TuneZ2_M-50_7TeV_madgraph_tauola_20121221.root",
         eventvec_MC, numevents, true);
   fitter.MatchMCjets( eventvec_MC );

   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Zmumu/"
         "Zmumu_data_DoubleMu_Run2011A_08Nov2011_v1_20121221.root",
         eventvec_data, numevents/2, false);
   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Zmumu/"
         "Zmumu_data_DoubleMu_Run2011B_19Nov2011_v1_20121221.root",
         eventvec_data, numevents/2, false);
   
   cout << "\n  MC EVENTS: " << eventvec_MC.size() << endl;
   cout << "DATA EVENTS: " << eventvec_data.size() << endl;

   double parMC [] =   {1.14753,1.12098,1.18506,1.20834,1.55375,2.56196,0.53392,-0.193797};
   double parData [] = {1.27509,1.24593,1.29233,1.34648,1.86162,2.47271,0.577556,-1.43896};

   vector<event> eventvec_sigmaMC = eventvec_MC;
   vector<event> eventvec_sigmaData = eventvec_MC;

   // reweight MC for pseudojet scalar pt
   double pjet_weights [] = {0.0,1.17948,1.12869,1.10327,1.09697,1.08439,1.07407,1.06587,1.05725,
      1.04646,1.0377,1.03198,1.02595,1.02087,1.0134,1.01204,1.0073,1.0009,0.997246,0.994572,0.984984,
      0.98568,0.989608,0.973945,0.999722,1.00649,0.946263,0.912974,1.00385,0.805591};
   fitter.PJetReweight( eventvec_MC, pjet_weights );
   fitter.PJetReweight( eventvec_sigmaData, pjet_weights );

   fitter.FindSignificance(parMC, eventvec_sigmaMC);
   fitter.FindSignificance(parData, eventvec_sigmaData);

   for( int i=0; i < eventvec_MC.size(); i++ ){
     eventvec_MC[i].met_varx = eventvec_sigmaData[i].cov_xx - eventvec_sigmaMC[i].cov_xx;
     eventvec_MC[i].met_vary = eventvec_sigmaData[i].cov_yy - eventvec_sigmaMC[i].cov_yy;
     eventvec_MC[i].met_rho = (eventvec_sigmaData[i].cov_xy - eventvec_sigmaMC[i].cov_xy)
        / sqrt(eventvec_MC[i].met_varx * eventvec_MC[i].met_vary);
   }

   fitter.FindSignificance(parData, eventvec_MC);
   fitter.FindSignificance(parData, eventvec_data);

   fitter.PlotsDataMC( eventvec_data, eventvec_MC, "results/plotsDataMC.root" );

   return 0;
}
