#include "METSigFit.h"

#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TLegend.h"

#include <iostream>
#include <unistd.h>
using namespace std;

struct ROCPoint {
   double cut;
   double pass;
   double total;
   ROCPoint( double c, double p, double t ) : cut(c), pass(p), total(t) {}
};

struct Dataset {
   // dataset info
   string path;
   string date;
   string channel;
   string filename;
   string process;
   bool isMC;
   int size;

   // for ROC curve
   vector<ROCPoint> ROCmet;
   vector<ROCPoint> ROCmetsig2011;
   vector<ROCPoint> ROCmetsig2012;
   vector<ROCPoint> ROCmetrht;

   // constructor
   Dataset( string f, string p, bool i) : filename(f), process(p), isMC(i) {}
};

int main(int argc, char* argv[]){

   // option flags
   char c;
   double fracevents = 1;
   bool do_resp_correction = false;
   string channel = "Zmumu";
   bool smear_met = true;
   string fileout = "results/plotsDataMC.root";
   bool compute_roc = false;
   bool run_data = true;
   bool run_mc = true;
   int met_type = 4;
   double rebin = 1;
   bool fullshape = false;
   double jec_var = 0;

   while( (c = getopt(argc, argv, "n:p:o:t:b:v:fhscbmrdw")) != -1 ) {
      switch(c)
      {
         case 'n' :
            fracevents = atof(optarg);
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

         case 'o':
            fileout = optarg;
            break;

         case 'r':
            compute_roc = true;
            break;

         case 'd':
            run_data = false;
            break;

         case 'w':
            run_mc = false;
            break;

         case 't':
            met_type = atoi(optarg);
            break;

         case 'b':
            rebin = atof(optarg);
            break;

         case 'f':
            fullshape = true;
            break;

         case 'v':
            jec_var = atof(optarg);
            break;

         case 'h' :
            cout << "Usage: ./EvalSig <flags>\n";
            cout << "Flags: \n";
            cout << "\t-n <number>\t  Fraction of events to fit.  Default at -1.\n";
            cout << "\t-j <number>\t  Jet bin pt threshold.  Default at 20 GeV.\n";
            cout << "\t-s\t          'Short' run, 10%% of events.\n";
            cout << "\t-c\t          Apply response correction.\n";
            cout << "\t-p <string>\t  Physics channel: Zmumu or Wenu.\n";
            cout << "\t-o <string>\t  Filename for Data/MC plots.\n";
            cout << "\t-m\t          Turn off MET smearing.\n";
            cout << "\t-r\t          Compute ROC curve.\n";
            cout << "\t-d\t          Do not run on data.\n";
            cout << "\t-w\t          Do not run on MC.\n";
            cout << "\t-t <number>\t MET type, in range [-1,4].\n";
            cout << "\t-b <number>\t Rebin -- divide bins by number.\n";
            cout << "\t-f\t          Compute Significance with full jet resolution shapes.\n";
            cout << "\t-v\t          Scale up (1) or down (-1) by JEC uncertainty.\n";
            cout << "\t-h\t          Display this menu.\n";
            return -1;
            break;

         default :
            continue;
      }
   }

   // declarations
   Fitter fitter(rebin);
   vector<Dataset> datasets;


   //
   // get all ntuples
   //

   if( channel.compare("Wenu") == 0 or channel.compare("Wenu_loose") == 0 ){

      // data
      if( run_data ){
         datasets.push_back( Dataset("Run2012A-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012B-1-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012B-2-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012B-3-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012B-4-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012B-5-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012B-6-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012B-7-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012B-8-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012B-9-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012B-10-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012C-1-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012C-2-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012C-3-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012C-4-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012C-5-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012C-6-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012C-7-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012C-8-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012C-9-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012C-10-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012D-1-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012D-2-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012D-3-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012D-4-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012D-5-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012D-6-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012D-7-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012D-8-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012D-9-22Jan2013/*.root", "Data", false));
         datasets.push_back( Dataset("Run2012D-10-22Jan2013/*.root", "Data", false));
      }

      // mc
      if( run_mc ){
         datasets.push_back( Dataset("DYJetsToLL_M-50/*.root", "DYJetsToLL", true) );
         datasets.push_back( Dataset("DYJetsToLL_M-10To50/*.root", "DYJetsToLL_M10To50", true) );
         datasets.push_back( Dataset("QCD_EMEnriched_20_30/*.root", "QCD_EMEnriched_20_30", true) );
         datasets.push_back( Dataset("QCD_EMEnriched_30_80/*.root", "QCD_EMEnriched_30_80", true) );
         datasets.push_back( Dataset("QCD_EMEnriched_80_170/*.root", "QCD_EMEnriched_80_170", true) );
         datasets.push_back( Dataset("QCD_EMEnriched_170_250/*.root", "QCD_EMEnriched_170_250", true) );
         datasets.push_back( Dataset("QCD_EMEnriched_250_350/*.root", "QCD_EMEnriched_250_350", true) );
         datasets.push_back( Dataset("QCD_EMEnriched_350/*.root", "QCD_EMEnriched_350", true) );
         datasets.push_back( Dataset("QCD_BCtoE_20_30/*.root", "QCD_BCtoE_20_30", true) );
         datasets.push_back( Dataset("QCD_BCtoE_30_80/*.root", "QCD_BCtoE_30_80", true) );
         datasets.push_back( Dataset("QCD_BCtoE_80_170/*.root", "QCD_BCtoE_80_170", true) );
         datasets.push_back( Dataset("QCD_BCtoE_170_250/*.root", "QCD_BCtoE_170_250", true) );
         datasets.push_back( Dataset("QCD_BCtoE_250_350/*.root", "QCD_BCtoE_250_350", true) );
         datasets.push_back( Dataset("Gamma_0_15/*.root", "Gamma_0_15", true) );
         datasets.push_back( Dataset("Gamma_15_30/*.root", "Gamma_15_30", true) );
         datasets.push_back( Dataset("Gamma_30_50/*.root", "Gamma_30_50", true) );
         datasets.push_back( Dataset("Gamma_50_80/*.root", "Gamma_50_80", true) );
         datasets.push_back( Dataset("Gamma_80_120/*.root", "Gamma_80_120", true) );
         datasets.push_back( Dataset("Gamma_120_170/*.root", "Gamma_120_170", true) );
         datasets.push_back( Dataset("Gamma_170_300/*.root", "Gamma_170_300", true) );
         datasets.push_back( Dataset("Gamma_300_470/*.root", "Gamma_300_470", true) );
         datasets.push_back( Dataset("TTJets/*.root", "TTJets", true) );
         datasets.push_back( Dataset("Tbar_tW-channel/*.root", "Tbar_tW", true) );
         datasets.push_back( Dataset("T_tW-channel/*.root", "T_tW", true) );
         datasets.push_back( Dataset("WJetsToLNu/*.root", "WJetsToLNu", true) );
         datasets.push_back( Dataset("WW/*.root", "WW", true) );
         datasets.push_back( Dataset("WZ/*.root", "WZ", true) );
         datasets.push_back( Dataset("ZZ/*.root", "ZZ", true) );
      }

   }
   else if( channel.compare("Zmumu") == 0 ){

      // data
      if( run_data ){
         datasets.push_back( Dataset("Run2012A-22Jan2013/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012B-part1-22Jan2013/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012B-part2-22Jan2013/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012B-part3-22Jan2013/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012B-part4-22Jan2013/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012C-part1-22Jan2013/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012C-part2-22Jan2013/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012C-part3-22Jan2013/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012C-part4-22Jan2013/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012D-part1-22Jan2013/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012D-part2-22Jan2013/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012D-part3-22Jan2013/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012D-part4-22Jan2013/*.root", "Data", false) );
      }

      // mc
      if( run_mc ){
         datasets.push_back( Dataset("DYJetsToLL/*.root", "DYJetsToLL", true) );
         datasets.push_back( Dataset("TTJets/*.root", "TTJets", true) );
         datasets.push_back( Dataset("Tbar_tW-channel/*.root", "Tbar_tW", true) );
         datasets.push_back( Dataset("T_tW-channel/*.root", "T_tW", true) );
         datasets.push_back( Dataset("WJetsToLNu/*.root", "WJetsToLNu", true) );
         datasets.push_back( Dataset("WW/*.root", "WW", true) );
         datasets.push_back( Dataset("WZ/*.root", "WZ", true) );
         datasets.push_back( Dataset("ZZ/*.root", "ZZ", true) );
      }

   }
   else if( channel.compare("Dijet") == 0 ){

      // data
      if( run_data ){
         datasets.push_back( Dataset("Run2012A-22Jan2013/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012B-part1-22Jan2013/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012B-part2-22Jan2013/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012B-part3-22Jan2013/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012B-part4-22Jan2013/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012C-part1-22Jan2013/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012C-part2-22Jan2013/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012C-part3-22Jan2013/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012C-part4-22Jan2013/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012D-part1-22Jan2013/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012D-part2-22Jan2013/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012D-part3-22Jan2013/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012D-part4-22Jan2013/*.root", "Data", false) );
      }

      // mc
      if( run_mc ){
         datasets.push_back( Dataset("QCD_15_30/*.root", "QCD_15_30", true) );
         datasets.push_back( Dataset("QCD_30_50/*.root", "QCD_30_50", true) );
         datasets.push_back( Dataset("QCD_50_80/*.root", "QCD_50_80", true) );
         datasets.push_back( Dataset("QCD_80_120/*.root", "QCD_80_120", true) );
         datasets.push_back( Dataset("QCD_120_170/*.root", "QCD_120_170", true) );
         datasets.push_back( Dataset("QCD_170_300/*.root", "QCD_170_300", true) );
         datasets.push_back( Dataset("QCD_300_470/*.root", "QCD_300_470", true) );
         datasets.push_back( Dataset("QCD_470_600/*.root", "QCD_470_600", true) );
         datasets.push_back( Dataset("QCD_600_800/*.root", "QCD_600_800", true) );
         datasets.push_back( Dataset("QCD_800_1000/*.root", "QCD_800_1000", true) );
         datasets.push_back( Dataset("QCD_1000_1400/*.root", "QCD_1000_1400", true) );
         datasets.push_back( Dataset("DYJetsToLL/*.root", "DYJetsToLL", true) );
         datasets.push_back( Dataset("TTJets/*.root", "TTJets", true) );
         datasets.push_back( Dataset("Tbar_tW-channel/*.root", "Tbar_tW", true) );
         datasets.push_back( Dataset("T_tW-channel/*.root", "T_tW", true) );
         datasets.push_back( Dataset("WJetsToLNu/*.root", "WJetsToLNu", true) );
         datasets.push_back( Dataset("WW/*.root", "WW", true) );
         datasets.push_back( Dataset("WZ/*.root", "WZ", true) );
         datasets.push_back( Dataset("ZZ/*.root", "ZZ", true) );
      }

   }
   else if( channel.compare("Ttbar0lept") == 0 ){

      // data
      if( run_data ){
         datasets.push_back( Dataset("Run2012A-13Jul2012/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012A-recover-06Aug2012/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012B-13Jul2012/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012C-24Aug2012/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012C-PromptReco/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012D-PromptReco/*.root", "Data", false) );
      }

      // mc
      if( run_mc ){
         datasets.push_back( Dataset("TTJets_Hadronic/*.root", "TTJets_Hadronic", true) );
         datasets.push_back( Dataset("TTJets_FullLept/*.root", "TTJets_FullLept", true) );
         datasets.push_back( Dataset("TTJets_SemiLept/*.root", "TTJets_SemiLept", true) );
         datasets.push_back( Dataset("DYJetsToLL/*.root", "DYJetsToLL", true) );
         datasets.push_back( Dataset("Tbar_tW-channel/*.root", "Tbar_tW", true) );
         datasets.push_back( Dataset("T_tW-channel/*.root", "T_tW", true) );
         datasets.push_back( Dataset("WJetsToLNu/*.root", "WJetsToLNu", true) );
         datasets.push_back( Dataset("QCD_15_30/*.root", "QCD_15_30", true) );
         datasets.push_back( Dataset("QCD_30_50/*.root", "QCD_30_50", true) );
         datasets.push_back( Dataset("QCD_50_80/*.root", "QCD_50_80", true) );
         datasets.push_back( Dataset("QCD_80_120/*.root", "QCD_80_120", true) );
         datasets.push_back( Dataset("QCD_120_170/*.root", "QCD_120_170", true) );
         datasets.push_back( Dataset("QCD_170_300/*.root", "QCD_170_300", true) );
         datasets.push_back( Dataset("QCD_300_470/*.root", "QCD_300_470", true) );
         datasets.push_back( Dataset("QCD_470_600/*.root", "QCD_470_600", true) );
         datasets.push_back( Dataset("QCD_600_800/*.root", "QCD_600_800", true) );
         datasets.push_back( Dataset("QCD_800_1000/*.root", "QCD_800_1000", true) );
         datasets.push_back( Dataset("QCD_1000_1400/*.root", "QCD_1000_1400", true) );
      }

   }
   else if( channel.compare("Ttbar1lept") == 0 ){

      // data
      if( run_data ){
         datasets.push_back( Dataset("Run2012A-22Jan2013-SingleMu/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012A-22Jan2013-SingleElectron/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012B-22Jan2013-SingleMu/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012B-22Jan2013-SingleElectron/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012C-22Jan2013-SingleMu/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012C-22Jan2013-SingleElectron/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012D-22Jan2013-SingleMu/*.root", "Data", false) );
         datasets.push_back( Dataset("Run2012D-22Jan2013-SingleElectron/*.root", "Data", false) );
      }

      // mc
      if( run_mc ){
         datasets.push_back( Dataset("TTJets_FullLept/*.root", "TTJets_FullLept", true) );
         datasets.push_back( Dataset("TTJets_SemiLept/*.root", "TTJets_SemiLept", true) );
         datasets.push_back( Dataset("TTJets_Hadronic/*.root", "TTJets_Hadronic", true) );
         datasets.push_back( Dataset("DYJetsToLL_M-50/*.root", "DYJetsToLL", true) );
         datasets.push_back( Dataset("DYJetsToLL_M-10To50/*.root", "DYJetsToLL_M10To50", true) );
         datasets.push_back( Dataset("Tbar_tW-channel/*.root", "Tbar_tW", true) );
         datasets.push_back( Dataset("T_tW-channel/*.root", "T_tW", true) );
         datasets.push_back( Dataset("WJetsToLNu/*.root", "WJetsToLNu", true) );
         datasets.push_back( Dataset("WW/*.root", "WW", true) );
         datasets.push_back( Dataset("WZ/*.root", "WZ", true) );
         datasets.push_back( Dataset("ZZ/*.root", "ZZ", true) );
         datasets.push_back( Dataset("QCD_EMEnriched_20_30/*.root", "QCD_EMEnriched_20_30", true) );
         datasets.push_back( Dataset("QCD_EMEnriched_30_80/*.root", "QCD_EMEnriched_30_80", true) );
         datasets.push_back( Dataset("QCD_EMEnriched_80_170/*.root", "QCD_EMEnriched_80_170", true) );
         datasets.push_back( Dataset("QCD_EMEnriched_170_250/*.root", "QCD_EMEnriched_170_250", true) );
         datasets.push_back( Dataset("QCD_EMEnriched_250_350/*.root", "QCD_EMEnriched_250_350", true) );
         datasets.push_back( Dataset("QCD_EMEnriched_350/*.root", "QCD_EMEnriched_350", true) );
         datasets.push_back( Dataset("QCD_BCtoE_20_30/*.root", "QCD_BCtoE_20_30", true) );
         datasets.push_back( Dataset("QCD_BCtoE_30_80/*.root", "QCD_BCtoE_30_80", true) );
         datasets.push_back( Dataset("QCD_BCtoE_80_170/*.root", "QCD_BCtoE_80_170", true) );
         datasets.push_back( Dataset("QCD_BCtoE_170_250/*.root", "QCD_BCtoE_170_250", true) );
         datasets.push_back( Dataset("QCD_BCtoE_250_350/*.root", "QCD_BCtoE_250_350", true) );
         datasets.push_back( Dataset("Gamma_0_15/*.root", "Gamma_0_15", true) );
         datasets.push_back( Dataset("Gamma_15_30/*.root", "Gamma_15_30", true) );
         datasets.push_back( Dataset("Gamma_30_50/*.root", "Gamma_30_50", true) );
         datasets.push_back( Dataset("Gamma_50_80/*.root", "Gamma_50_80", true) );
         datasets.push_back( Dataset("Gamma_80_120/*.root", "Gamma_80_120", true) );
         datasets.push_back( Dataset("Gamma_120_170/*.root", "Gamma_120_170", true) );
         datasets.push_back( Dataset("Gamma_170_300/*.root", "Gamma_170_300", true) );
         datasets.push_back( Dataset("Gamma_300_470/*.root", "Gamma_300_470", true) );
      }

   }
   else{ cout << "Unknown physics channel.  Use option 'p' to input channel name." << endl; }

   for( vector<Dataset>::iterator data = datasets.begin(); data != datasets.end(); data++ ){
      data->channel = channel;
      data->path = "/mnt/xrootd/user/nmirman/Ntuples/METsig";
      data->date = "20130830";
      if( channel.compare("Zmumu") == 0 and data->isMC ){
         data->date = "20130913";
      }
      if( channel.compare("Wenu") == 0 ){
         data->date = "20130916";
      }
      if( channel.compare("Dijet") == 0 and data->isMC ){
         data->date = "20130913";
      }
      if( channel.compare("Ttbar0lept") == 0 ){
         data->date = "20130913";
      }
   }

   // get number of events in datasets
   cout << "Getting number of events in datasets." << endl;
   for( vector<Dataset>::iterator data = datasets.begin(); data != datasets.end(); data++ ){
      //TFile *file = TFile::Open(data->filename);
      //if( !file ){
      //   data->size = 0;
      //   continue;
      //}
      //Tree *tree = (TTree*)file->Get("events");
      string fullname = data->path+"/"+data->channel+"/"+data->date+"/"+data->filename;
      TChain tree("events");
      tree.Add( fullname.c_str() );
      data->size = tree.GetEntries();
      cout << fullname << ": " << data->size << " events." << endl;
   }
   cout << endl;

   //
   // loop through datasets, fill histograms
   //
   //double parMC [] =   {1.12660,1.09322,1.10951,1.17178,1.12164,0.0,0.585145};
   //double parData [] = {1.39669,1.32037,1.32047,1.38161,1.51508,0.0,0.639158};
   //double parData [] =   {1.29446,1.24207,1.26686,1.34076,1.49548,0.0,0.6117};
   //double parMC   [] =   {1.11659,1.06256,1.09741,1.11931,1.17266,0.0,0.569454};
   double parData [] = {1.15061,1.07776,1.04204,1.12509,1.56414,0.0,0.548758};
   double parMC [] = {1.05347,0.975375,0.957986,0.97269,1.28106,-1.10982,0.52039};

   double parData_up [] = {1.15061,1.07776,1.04204,1.12509,1.56414,0.0,0.548758};
   //double parMC_up [] = {1.05347,0.975375,0.957986,0.97269,1.28106,-1.10982,0.52039};
   double parMC_up [] = {1.07177,1.00423,0.979847,1.00971,1.37391,-0.0128714,0.558511};

   double parData_down [] = {1.15061,1.07776,1.04204,1.12509,1.56414,0.0,0.548758};
   //double parMC_down [] = {1.05347,0.975375,0.957986,0.97269,1.28106,-1.10982,0.52039};
   double parMC_down [] = {1.03122,0.951483,0.930119,0.92954,1.17063,-2.51767,0.486471};

   fitter.met_type = met_type;

   for( vector<Dataset>::iterator data = datasets.begin(); data != datasets.end(); data++ ){

      // initialize counters for ROC curve
      if( compute_roc and data->isMC ){
         for(int i=0; i < 100; i++){
            data->ROCmet.push_back( ROCPoint(2*i, 0, 0) );
            data->ROCmetsig2011.push_back( ROCPoint(exp(double(i)/5-15), 0, 0) );
            data->ROCmetsig2012.push_back( ROCPoint(exp(double(i)/5-15), 0, 0) );
            data->ROCmetrht.push_back( ROCPoint(exp(double(i)/5-15), 0, 0) );
         }
      }

      //int section_size = 1000000;
      int section_size = 500000;
      int num_events = fracevents*data->size;
      int num_sections = 1 + ((num_events-1)/section_size);
      cout << "Opening dataset " << data->filename << endl;
      cout << "Divide " << num_events << " events into " << num_sections << " sections..." << endl;

      for(int isec=0; isec < num_sections; isec++){
         int start = isec*section_size;
         int end = (isec == num_sections-1) ? num_events : start + section_size;
         cout << "Begin section [" << start << ", " << end << "]" << endl;

         double jec_varmc = data->isMC ? jec_var : 0;

         vector<event> eventvec;
         string fullname = data->path+"/"+data->channel+"/"+data->date+"/"+data->filename;
         fitter.ReadNtuple( fullname.c_str(), eventvec, 1,
               data->isMC, data->process, do_resp_correction, start, end, jec_varmc );

         vector<event> eventvec_sigmaMC;

         // met smearing for mc datasets
         if( data->isMC and smear_met ){
            eventvec_sigmaMC = eventvec;

            if( jec_var == 1 ){
               fitter.FindSignificance(parMC_up, eventvec_sigmaMC);
               fitter.FindSignificance(parData_up, eventvec);
            } else if ( jec_var == -1 ){
               fitter.FindSignificance(parMC_down, eventvec_sigmaMC);
               fitter.FindSignificance(parData_down, eventvec);
            } else {
               fitter.FindSignificance(parMC, eventvec_sigmaMC);
               fitter.FindSignificance(parData, eventvec);
            }

            for( int i=0; i < int(eventvec.size()); i++ ){
               eventvec[i].met_varx = eventvec[i].cov_xx - eventvec_sigmaMC[i].cov_xx;
               eventvec[i].met_vary = eventvec[i].cov_yy - eventvec_sigmaMC[i].cov_yy;
               eventvec[i].met_rho = (eventvec[i].cov_xy - eventvec_sigmaMC[i].cov_xy)
                  / sqrt(eventvec[i].met_varx * eventvec[i].met_vary);
            }
         }

         // compute significance
         if( !fullshape ){
            if( !(data->isMC) or smear_met ){
               if( jec_var == 1 ){
                  fitter.FindSignificance(parData_up, eventvec);
               } else if ( jec_var == -1 ){
                  fitter.FindSignificance(parData_down, eventvec);
               } else {
                  fitter.FindSignificance(parData, eventvec);
               }
            }else{
               if( jec_var == 1 ){
                  fitter.FindSignificance(parMC_up, eventvec);
               } else if ( jec_var == -1 ){
                  fitter.FindSignificance(parMC_down, eventvec);
               } else {
                  fitter.FindSignificance(parMC, eventvec);
               }
            }
         }else{
            if( !(data->isMC) or smear_met ){
               fitter.FullShapeSig(parData, eventvec);
            }else{
               fitter.FullShapeSig(parMC, eventvec);
            }
         }

         // compute average met
         double norm = 0;
         double met = 0;
         for( int i=0; i < int(eventvec.size()); i++ ){
            met += eventvec[i].met * eventvec[i].weight;
            norm += eventvec[i].weight;
         }
         cout << "Avg MET = " << met << "/" << norm << " = " << met/norm << endl;

         // fill histograms
         fitter.FillHists(eventvec, channel); 

         // ROC curve
         if( compute_roc and data->isMC ){
            for(vector<event>::iterator ev = eventvec.begin(); ev < eventvec.end(); ev++){
               // met
               for( int i = 0; i < int(data->ROCmet.size()); i++ ){
                  data->ROCmet[i].total += ev->weight;
                  if( data->channel.compare("Ttbar1lept") == 0
                        or data->channel.compare("Ttbar0lept") == 0
                        or data->channel.compare("Dijet") == 0 ){
                     if( ev->met < data->ROCmet[i].cut ) data->ROCmet[i].pass += ev->weight;
                  }else{
                     if( ev->met > data->ROCmet[i].cut ) data->ROCmet[i].pass += ev->weight;
                  }
               }
               // metsig2011
               for( int i = 0; i < int(data->ROCmetsig2011.size()); i++ ){
                  data->ROCmetsig2011[i].total += ev->weight;
                  if( data->channel.compare("Ttbar1lept") == 0
                        or data->channel.compare("Ttbar0lept") == 0
                        or data->channel.compare("Dijet") == 0 ){
                     if(ev->metsig2011 < data->ROCmetsig2011[i].cut)
                        data->ROCmetsig2011[i].pass += ev->weight;
                  }else{
                     if(ev->metsig2011 > data->ROCmetsig2011[i].cut)
                        data->ROCmetsig2011[i].pass += ev->weight;
                  }
               }
               // metsig2012
               for( int i = 0; i < int(data->ROCmetsig2012.size()); i++ ){
                  data->ROCmetsig2012[i].total += ev->weight;
                  if( data->channel.compare("Ttbar1lept") == 0
                        or data->channel.compare("Ttbar0lept") == 0
                        or data->channel.compare("Dijet") == 0 ){
                     if(ev->sig < data->ROCmetsig2012[i].cut) data->ROCmetsig2012[i].pass += ev->weight;
                  }else{
                     if(ev->sig > data->ROCmetsig2012[i].cut) data->ROCmetsig2012[i].pass += ev->weight;
                  }
               }
               // dumb metsig
               double ht = ev->pjet_scalptL123;
               for(int i=0; i < int(ev->jet_ptUncor.size()); i++){
                  ht += ev->jet_ptL123[i];
               }
               double metrht = pow(ev->met,2)/ht;
               for( int i = 0; i < int(data->ROCmetrht.size()); i++ ){
                  data->ROCmetrht[i].total += ev->weight;
                  if( data->channel.compare("Ttbar1lept") == 0
                        or data->channel.compare("Ttbar0lept") == 0
                        or data->channel.compare("Dijet") == 0 ){
                     if(metrht < data->ROCmetrht[i].cut) data->ROCmetrht[i].pass += ev->weight;
                  }else{
                     if(metrht > data->ROCmetrht[i].cut) data->ROCmetrht[i].pass += ev->weight;
                  }
               }

            }
         } // ROC curve

      }

   }

   // combine all channels, print histograms
   fitter.PrintHists(fileout.c_str(), channel);

   // ROC plots
   if( compute_roc ){
      TGraph* gROCmet = new TGraph();
      TGraph* gROCmetsig2011 = new TGraph();
      TGraph* gROCmetsig2012 = new TGraph();
      TGraph* gROCmetrht = new TGraph();

      vector<Dataset>::iterator datatemp = datasets.end() - 1;
      vector<double> met_sigpass(datatemp->ROCmet.size(),0);
      vector<double> met_bkgpass(datatemp->ROCmet.size(),0);
      vector<double> met_sigtot(datatemp->ROCmet.size(),0);
      vector<double> met_bkgtot(datatemp->ROCmet.size(),0);

      vector<double> metsig2011_sigpass(datatemp->ROCmetsig2011.size(),0);
      vector<double> metsig2011_bkgpass(datatemp->ROCmetsig2011.size(),0);
      vector<double> metsig2011_sigtot(datatemp->ROCmetsig2011.size(),0);
      vector<double> metsig2011_bkgtot(datatemp->ROCmetsig2011.size(),0);

      vector<double> metsig2012_sigpass(datatemp->ROCmetsig2012.size(),0);
      vector<double> metsig2012_bkgpass(datatemp->ROCmetsig2012.size(),0);
      vector<double> metsig2012_sigtot(datatemp->ROCmetsig2012.size(),0);
      vector<double> metsig2012_bkgtot(datatemp->ROCmetsig2012.size(),0);

      vector<double> metrht_sigpass(datatemp->ROCmetrht.size(),0);
      vector<double> metrht_bkgpass(datatemp->ROCmetrht.size(),0);
      vector<double> metrht_sigtot(datatemp->ROCmetrht.size(),0);
      vector<double> metrht_bkgtot(datatemp->ROCmetrht.size(),0);

      for( vector<Dataset>::iterator data = datasets.begin(); data != datasets.end(); data++ ){
         if( data->isMC ){

            if( (data->channel.compare("Wenu") == 0 and data->process.compare("WJetsToLNu") == 0)
                  or (data->channel.compare("Ttbar1lept") == 0 and data->process.compare("TTJets_SemiLept") == 0)
                  or (data->channel.compare("Ttbar0lept") == 0 and data->process.compare("TTJets_Hadronic") == 0)
                  or (data->channel.compare("Dijet") == 0 and data->process.find("QCD") != string::npos)  ) {
               for(int i=0; i < int(data->ROCmet.size()); i++){ // met
                  met_sigpass[i] += data->ROCmet[i].pass;
                  met_sigtot[i] += data->ROCmet[i].total;
               }
               for(int i=0; i < int(data->ROCmetsig2011.size()); i++){ // metsig2011
                  metsig2011_sigpass[i] += data->ROCmetsig2011[i].pass;
                  metsig2011_sigtot[i] += data->ROCmetsig2011[i].total;
               }
               for(int i=0; i < int(data->ROCmetsig2012.size()); i++){ // metsig2012
                  metsig2012_sigpass[i] += data->ROCmetsig2012[i].pass;
                  metsig2012_sigtot[i] += data->ROCmetsig2012[i].total;
               }
               for(int i=0; i < int(data->ROCmetrht.size()); i++){ // metrht
                  metrht_sigpass[i] += data->ROCmetrht[i].pass;
                  metrht_sigtot[i] += data->ROCmetrht[i].total;
               }
            //}else if( string::npos != string(data->channel).find("DY")
            //      or string::npos != string(data->channel).find("QCD") 
            //      or string::npos != string(data->channel).find("Gamma")
            //      ){
            }else{
               for(int i=0; i < int(data->ROCmet.size()); i++){ // met
                  met_bkgpass[i] += data->ROCmet[i].pass;
                  met_bkgtot[i] += data->ROCmet[i].total;
               }
               for(int i=0; i < int(data->ROCmetsig2011.size()); i++){ // metsig2011
                  metsig2011_bkgpass[i] += data->ROCmetsig2011[i].pass;
                  metsig2011_bkgtot[i] += data->ROCmetsig2011[i].total;
               }
               for(int i=0; i < int(data->ROCmetsig2012.size()); i++){ // metsig2012
                  metsig2012_bkgpass[i] += data->ROCmetsig2012[i].pass;
                  metsig2012_bkgtot[i] += data->ROCmetsig2012[i].total;
               }
               for(int i=0; i < int(data->ROCmetrht.size()); i++){ // metrht
                  metrht_bkgpass[i] += data->ROCmetrht[i].pass;
                  metrht_bkgtot[i] += data->ROCmetrht[i].total;
               }
            }

         }
      } // loop through datasets

      for(int i=0; i < int(met_sigpass.size()); i++){ // met
         gROCmet->SetPoint(i, double(met_bkgpass[i])/met_bkgtot[i],
               double(met_sigpass[i])/met_sigtot[i]);
      }
      for(int i=0; i < int(metsig2011_sigpass.size()); i++){ // metsig2011
         gROCmetsig2011->SetPoint(i, double(metsig2011_bkgpass[i])/metsig2011_bkgtot[i],
               double(metsig2011_sigpass[i])/metsig2011_sigtot[i]);
      }
      for(int i=0; i < int(metsig2012_sigpass.size()); i++){ // metsig2012
         gROCmetsig2012->SetPoint(i, double(metsig2012_bkgpass[i])/metsig2012_bkgtot[i],
               double(metsig2012_sigpass[i])/metsig2012_sigtot[i]);
      }
      for(int i=0; i < int(metrht_sigpass.size()); i++){ // metrht
         gROCmetrht->SetPoint(i, double(metrht_bkgpass[i])/metrht_bkgtot[i],
               double(metrht_sigpass[i])/metrht_sigtot[i]);
      }

      TCanvas* cROC = new TCanvas("cROC","cROC",800,800);
      cROC->cd();

      gROCmet->SetTitle("ROC Curve;Background Efficiency;Signal Efficiency");
      gROCmet->GetXaxis()->SetTitleSize(0.07);
      gROCmet->GetYaxis()->SetTitleSize(0.07);
      gROCmet->GetXaxis()->SetTitleOffset(0.8);
      gROCmet->GetYaxis()->SetTitleOffset(1.0);
      gROCmet->GetXaxis()->SetLimits(0.01,1.05);

      gROCmet->SetMarkerStyle(20);
      gROCmet->SetMarkerSize(0.7);
      gROCmet->SetMarkerColor(1);
      gROCmet->SetLineColor(1);
      gROCmet->SetLineStyle(2);

      gROCmetsig2011->SetMarkerStyle(20);
      gROCmetsig2011->SetMarkerSize(0.7);
      gROCmetsig2011->SetMarkerColor(1);
      gROCmetsig2011->SetLineColor(1);

      gROCmetsig2012->SetMarkerStyle(20);
      gROCmetsig2012->SetMarkerSize(0.7);
      gROCmetsig2012->SetMarkerColor(2);
      gROCmetsig2012->SetLineColor(2);

      gROCmetrht->SetMarkerStyle(20);
      gROCmetrht->SetMarkerSize(0.7);
      gROCmetrht->SetMarkerColor(4);
      gROCmetrht->SetLineColor(4);

      gROCmet->Draw("ACP");
      gROCmetsig2012->Draw("CP");
      gROCmetsig2011->Draw("CP");
      gROCmetrht->Draw("CP");

      TF1* fline = new TF1("fline", "x", 0, 1);
      fline->SetLineColor(1);
      fline->SetLineStyle(7);
      fline->Draw("same");

      TLegend *lROC = new TLegend(0.605528,0.655866,0.866834,0.816333);
      lROC->AddEntry(gROCmet,"met","lp");
      lROC->AddEntry(gROCmetsig2011,"metsig2011","lp");
      lROC->AddEntry(gROCmetsig2012,"metsig2012","lp");
      lROC->AddEntry(gROCmetrht,"met/#sqrt{H_{T}}","lp");
      lROC->Draw("same");

      TFile *file = new TFile(fileout.c_str(),"UPDATE");
      file->cd();
      cROC->Write();
      file->Close();

      delete cROC;
      delete lROC;
      delete file;
   } // ROC plots

   return 0;
}
