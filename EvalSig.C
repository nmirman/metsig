#include "METSigFit.h"

#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include <iostream>
#include <unistd.h>
using namespace std;

struct Dataset {
   const char* filename;
   const char* channel;
   bool isMC;
   Dataset( const char* f, const char* c, bool i) : filename(f), channel(c), isMC(i) {}
};

int main(int argc, char* argv[]){

   int numevents = -1;
   bool stackMC = true;
   bool do_resp_correction = false;

   // declarations
   Fitter fitter;
   vector<Dataset> datasets;

   //
   // get all ntuples
   //

   // data
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/Run2012A-13Jul2012.root", "Data", false) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/Run2012A-recover-06Aug2012.root", "Data", false) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/Run2012B-13Jul2012.root", "Data", false) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/Run2012C-24Aug2012.root", "Data", false) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/Run2012C-EcalRecover_11Dec2012.root", "Data", false) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/Run2012C-PromptReco-v2.root", "Data", false) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/Run2012D-PromptReco-v1.root", "Data", false) );

   // mc
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/DYJetsToLL_M-50.root", "DYJetsToLL", true) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/DYJetsToLL_M-10To50.root", "DYJetsToLL_M10To50", true) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/QCD_EMEnriched_20_30.root", "QCD_EMEnriched_20_30", true) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/QCD_EMEnriched_30_80.root", "QCD_EMEnriched_30_80", true) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/QCD_EMEnriched_80_170.root", "QCD_EMEnriched_80_170", true) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/QCD_BCtoE_20_30.root", "QCD_BCtoE_20_30", true) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/QCD_BCtoE_30_80.root", "QCD_BCtoE_30_80", true) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/QCD_BCtoE_80_170.root", "QCD_BCtoE_80_170", true) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/Gamma_0_15.root", "Gamma_0_15", true) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/Gamma_15_30.root", "Gamma_15_30", true) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/Gamma_30_50.root", "Gamma_30_50", true) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/Gamma_50_80.root", "Gamma_50_80", true) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/Gamma_80_120.root", "Gamma_80_120", true) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/Gamma_120_170.root", "Gamma_120_170", true) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/TTJets.root", "TTJets", true) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/Tbar_tW-channel.root", "Tbar_tW", true) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/T_tW-channel.root", "T_tW", true) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/WJetsToLNu.root", "WJetsToLNu", true) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/WW.root", "WW", true) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/WZ.root", "WZ", true) );
   datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130504/ZZ.root", "ZZ", true) );

   //
   // loop through datasets, fill histograms
   //
   double parMC [] =   {1.12660,1.09322,1.10951,1.17178,1.12164,0.0,0.585145};
   double parData [] = {1.39669,1.32037,1.32047,1.38161,1.51508,0.0,0.639158};

   for( vector<Dataset>::iterator data = datasets.begin(); data != datasets.end(); data++ ){

      vector<event> eventvec;
      fitter.ReadNtuple( data->filename, eventvec, numevents,
            data->isMC, data->channel, do_resp_correction );

      vector<event> eventvec_sigmaMC;
      vector<event> eventvec_sigmaData;

      // met smearing for mc datasets
      if( data->isMC ){
         eventvec_sigmaMC = eventvec;
         eventvec_sigmaData = eventvec;
         fitter.FindSignificance(parMC, eventvec_sigmaMC);
         fitter.FindSignificance(parData, eventvec_sigmaData);
         for( int i=0; i < int(eventvec.size()); i++ ){
            eventvec[i].met_varx = eventvec_sigmaData[i].cov_xx - eventvec_sigmaMC[i].cov_xx;
            eventvec[i].met_vary = eventvec_sigmaData[i].cov_yy - eventvec_sigmaMC[i].cov_yy;
            eventvec[i].met_rho = (eventvec_sigmaData[i].cov_xy - eventvec_sigmaMC[i].cov_xy)
               / sqrt(eventvec[i].met_varx * eventvec[i].met_vary);
         }
      }

      // compute significance
      fitter.FindSignificance(parData, eventvec);

      // fill histograms
      fitter.FillHists(eventvec, "Wenu"); 
   }

   // combine all channels, print histograms
   fitter.PrintHists("results/plotsDataMC.root", "Wenu");

   return 0;
}
