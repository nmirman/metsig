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
   int size;
   Dataset( const char* f, const char* c, bool i) : filename(f), channel(c), isMC(i) {}
};

int main(int argc, char* argv[]){

   // declarations
   Fitter fitter;
   vector<Dataset> datasets;

   // option flags
   char c;
   double fracevents = 1;
   bool do_resp_correction = false;
   char* channel = "Zmumu";
   bool smear_met = true;
   char* fileout = "results/plotsDataMC.root";

   while( (c = getopt(argc, argv, "n:j:p:hscobm")) != -1 ) {
      switch(c)
      {
         case 'n' :
            fracevents = atof(optarg);
            break;

         case 'j' :
            fitter.jetbinpt = atof(optarg);
            break;

         case 's' :
            fracevents = 0.1;
            break;

         case 'c' :
            do_resp_correction=true;
            break;

         case 'p':
            channel = optarg;
            break;

         case 'm':
            smear_met = false;
            break;

         case 'f':
            fileout = optarg;
            break;

         case 'h' :
            cout << "Usage: ./EvalSig <flags>\n";
            cout << "Flags: \n";
            cout << "\t-n <number>\t  Fraction of events to fit.  Default at -1.\n";
            cout << "\t-j <number>\t  Jet bin pt threshold.  Default at 20 GeV.\n";
            cout << "\t-s\t          'Short' run, 10%% of events.\n";
            cout << "\t-c\t          Apply response correction.\n";
            cout << "\t-p <string>\t  Physics channel: Zmumu or Wenu.\n";
            cout << "\t-f <string>\t  Filename for Data/MC plots.\n";
            cout << "\t-m\t          Turn off MET smearing.\n";
            cout << "\t-h\t          Display this menu.\n";
            return -1;
            break;

         default :
            continue;
      }
   }


   //
   // get all ntuples
   //

   // data
   if( strcmp(channel,"Wenu") == 0 ){
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/Run2012A-22Jan2013.root", "Data", false));
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/Run2012B-22Jan2013.root", "Data", false));
      //datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/Run2012C-22Jan2013.root", "Data", false));
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/Run2012D-22Jan2013.root", "Data", false));
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/Run2012C-part1-22Jan2013.root", "Data", false));
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/Run2012C-part2-22Jan2013.root", "Data", false));
      //datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/Run2012D-part1-22Jan2013.root", "Data", false));
      //datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/Run2012D-part2-22Jan2013.root", "Data", false));

      // mc
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/DYJetsToLL_M-50.root", "DYJetsToLL", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/DYJetsToLL_M-10To50.root", "DYJetsToLL_M10To50", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/QCD_EMEnriched_20_30.root", "QCD_EMEnriched_20_30", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/QCD_EMEnriched_30_80.root", "QCD_EMEnriched_30_80", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/QCD_EMEnriched_80_170.root", "QCD_EMEnriched_80_170", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/QCD_EMEnriched_170_250.root", "QCD_EMEnriched_170_250", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/QCD_EMEnriched_250_350.root", "QCD_EMEnriched_250_350", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/QCD_EMEnriched_350.root", "QCD_EMEnriched_350", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/QCD_BCtoE_20_30.root", "QCD_BCtoE_20_30", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/QCD_BCtoE_30_80.root", "QCD_BCtoE_30_80", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/QCD_BCtoE_80_170.root", "QCD_BCtoE_80_170", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/QCD_BCtoE_170_250.root", "QCD_BCtoE_170_250", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/QCD_BCtoE_250_350.root", "QCD_BCtoE_250_350", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/Gamma_0_15.root", "Gamma_0_15", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/Gamma_15_30.root", "Gamma_15_30", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/Gamma_30_50.root", "Gamma_30_50", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/Gamma_50_80.root", "Gamma_50_80", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/Gamma_80_120.root", "Gamma_80_120", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/Gamma_120_170.root", "Gamma_120_170", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/Gamma_170_300.root", "Gamma_170_300", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/Gamma_300_470.root", "Gamma_300_470", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/TTJets.root", "TTJets", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/Tbar_tW-channel.root", "Tbar_tW", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/T_tW-channel.root", "T_tW", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/WJetsToLNu.root", "WJetsToLNu", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/WW.root", "WW", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/WZ.root", "WZ", true) );
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Wenu/20130626/ZZ.root", "ZZ", true) );
   }
   else if( strcmp(channel,"Zmumu") == 0 ){

      // data
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130626/Run2012A-22Jan2013.root", "Data", false));
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130626/Run2012B-22Jan2013.root", "Data", false));
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130626/Run2012C-22Jan2013.root", "Data", false));
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130626/Run2012D-22Jan2013.root", "Data", false));

      // mc
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130626/DYJetsToLL.root", "DYJetsToLL", true));
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130626/TTJets.root", "TTJets", true));
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130626/Tbar_tW-channel.root", "Tbar_tW", true));
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130626/T_tW-channel.root", "T_tW", true));
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130626/WJetsToLNu.root", "WJetsToLNu", true));
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130626/WW.root", "WW", true));
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130626/WZ.root", "WZ", true));
      datasets.push_back( Dataset("/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130626/ZZ.root", "ZZ", true));

   }
   else{ cout << "Unknown physics channel.  Use option 'p' to input channel name." << endl; }

   // get number of events in datasets
   cout << "Getting number of events in datasets." << endl;
   for( vector<Dataset>::iterator data = datasets.begin(); data != datasets.end(); data++ ){
      TFile *file = TFile::Open(data->filename);
      if( !file ){
         data->size = 0;
         continue;
      }
      TTree *tree = (TTree*)file->Get("events");
      data->size = tree->GetEntries();
      cout << data->filename << ": " << data->size << " events." << endl;
   }
   cout << endl;

   //
   // loop through datasets, fill histograms
   //
   //double parMC [] =   {1.12660,1.09322,1.10951,1.17178,1.12164,0.0,0.585145};
   //double parData [] = {1.39669,1.32037,1.32047,1.38161,1.51508,0.0,0.639158};
   double parData [] =   {1.29446,1.24207,1.26686,1.34076,1.49548,0.0,0.6117};
   double parMC   [] =   {1.11659,1.06256,1.09741,1.11931,1.17266,0.0,0.569454};

   for( vector<Dataset>::iterator data = datasets.begin(); data != datasets.end(); data++ ){

      int section_size = 2000000;
      int num_events = fracevents*data->size;
      int num_sections = 1 + ((num_events-1)/section_size);
      cout << "Opening dataset " << data->filename << endl;
      cout << "Divide into " << num_sections << " sections..." << endl;

      for(int isec=0; isec < num_sections; isec++){
         int start = isec*section_size;
         int end = (isec == num_sections-1) ? num_events : start + section_size;
         cout << "Begin section [" << start << ", " << end << "]" << endl;

         vector<event> eventvec;
         fitter.ReadNtuple( data->filename, eventvec, 1,
               data->isMC, data->channel, do_resp_correction, start, end );

         vector<event> eventvec_sigmaMC;
         vector<event> eventvec_sigmaData;

         // met smearing for mc datasets
         if( data->isMC and smear_met ){
            eventvec_sigmaMC = eventvec;
            //eventvec_sigmaData = eventvec;
            fitter.FindSignificance(parMC, eventvec_sigmaMC);
            fitter.FindSignificance(parData, eventvec/*_sigmaData*/);
            for( int i=0; i < int(eventvec.size()); i++ ){
               eventvec[i].met_varx = eventvec/*_sigmaData*/[i].cov_xx - eventvec_sigmaMC[i].cov_xx;
               eventvec[i].met_vary = eventvec/*_sigmaData*/[i].cov_yy - eventvec_sigmaMC[i].cov_yy;
               eventvec[i].met_rho = (eventvec/*_sigmaData*/[i].cov_xy - eventvec_sigmaMC[i].cov_xy)
                  / sqrt(eventvec[i].met_varx * eventvec[i].met_vary);
            }
         }

         // compute significance
         if( !(data->isMC) or smear_met ){
            fitter.FindSignificance(parData, eventvec);
         }else{
            fitter.FindSignificance(parMC, eventvec);
         }

         // fill histograms
         fitter.FillHists(eventvec, channel); 
      }

   }

   // combine all channels, print histograms
   fitter.PrintHists("results/plotsDataMC.root", channel);

   return 0;
}
