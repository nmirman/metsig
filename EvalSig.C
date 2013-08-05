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
   int pass;
   int total;
   ROCPoint( double c, int p, int t ) : cut(c), pass(p), total(t) {}
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

   // constructor
   Dataset( string f, string p, bool i) : filename(f), process(p), isMC(i) {}
};

int main(int argc, char* argv[]){

   // declarations
   Fitter fitter;
   vector<Dataset> datasets;

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

   while( (c = getopt(argc, argv, "n:j:p:f:hscobmrdw")) != -1 ) {
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
      }

   }
   else if( channel.compare("Ttbar0lept") == 0 ){

      // data
      if( run_data ){
         datasets.push_back( Dataset("Run2012A-22Jan2013/*.root", "Data", false) );
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
      data->date = "20130723";
      if( channel.compare("Zmumu") == 0 or ((channel.compare("Dijet") == 0) and !(data->isMC)) ){
         data->date = "20130728";
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
   double parData [] =   {1.29446,1.24207,1.26686,1.34076,1.49548,0.0,0.6117};
   double parMC   [] =   {1.11659,1.06256,1.09741,1.11931,1.17266,0.0,0.569454};

   for( vector<Dataset>::iterator data = datasets.begin(); data != datasets.end(); data++ ){

      // initialize counters for ROC curve
      if( compute_roc and data->isMC ){
         for(int i=0; i < 100; i++){
            data->ROCmet.push_back( ROCPoint(i, 0, 0) );
            data->ROCmetsig2011.push_back( ROCPoint(exp(double(i)/5-10), 0, 0) );
            data->ROCmetsig2012.push_back( ROCPoint(exp(double(i)/5-10), 0, 0) );
         }
      }

      int section_size = 1000000;
      int num_events = fracevents*data->size;
      int num_sections = 1 + ((num_events-1)/section_size);
      cout << "Opening dataset " << data->filename << endl;
      cout << "Divide " << num_events << " events into " << num_sections << " sections..." << endl;

      for(int isec=0; isec < num_sections; isec++){
         int start = isec*section_size;
         int end = (isec == num_sections-1) ? num_events : start + section_size;
         cout << "Begin section [" << start << ", " << end << "]" << endl;

         vector<event> eventvec;
         string fullname = data->path+"/"+data->channel+"/"+data->date+"/"+data->filename;
         fitter.ReadNtuple( fullname.c_str(), eventvec, 1,
               data->isMC, data->process, do_resp_correction, start, end );

         vector<event> eventvec_sigmaMC;

         // met smearing for mc datasets
         if( data->isMC and smear_met ){
            eventvec_sigmaMC = eventvec;

            fitter.FindSignificance(parMC, eventvec_sigmaMC);
            fitter.FindSignificance(parData, eventvec);
            for( int i=0; i < int(eventvec.size()); i++ ){
               eventvec[i].met_varx = eventvec[i].cov_xx - eventvec_sigmaMC[i].cov_xx;
               eventvec[i].met_vary = eventvec[i].cov_yy - eventvec_sigmaMC[i].cov_yy;
               eventvec[i].met_rho = (eventvec[i].cov_xy - eventvec_sigmaMC[i].cov_xy)
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

         // ROC curve
         if( compute_roc and data->isMC ){
            for(vector<event>::iterator ev = eventvec.begin(); ev < eventvec.end(); ev++){
               // met
               for( int i = 0; i < int(data->ROCmet.size()); i++ ){
                  data->ROCmet[i].total += 1;
                  if( ev->met > data->ROCmet[i].cut ) data->ROCmet[i].pass += 1;
               }
               // metsig2011
               for( int i = 0; i < int(data->ROCmetsig2011.size()); i++ ){
                  data->ROCmetsig2011[i].total += 1;
                  if(ev->metsig2011 > data->ROCmetsig2011[i].cut) data->ROCmetsig2011[i].pass += 1;
               }
               // metsig2012
               for( int i = 0; i < int(data->ROCmetsig2012.size()); i++ ){
                  data->ROCmetsig2012[i].total += 1;
                  if(ev->sig > data->ROCmetsig2012[i].cut) data->ROCmetsig2012[i].pass += 1;
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

      vector<int> met_sigpass(datasets[0].ROCmet.size(),0);
      vector<int> met_bkgpass(datasets[0].ROCmet.size(),0);
      vector<int> met_sigtot(datasets[0].ROCmet.size(),0);
      vector<int> met_bkgtot(datasets[0].ROCmet.size(),0);

      vector<int> metsig2011_sigpass(datasets[0].ROCmetsig2011.size(),0);
      vector<int> metsig2011_bkgpass(datasets[0].ROCmetsig2011.size(),0);
      vector<int> metsig2011_sigtot(datasets[0].ROCmetsig2011.size(),0);
      vector<int> metsig2011_bkgtot(datasets[0].ROCmetsig2011.size(),0);

      vector<int> metsig2012_sigpass(datasets[0].ROCmetsig2012.size(),0);
      vector<int> metsig2012_bkgpass(datasets[0].ROCmetsig2012.size(),0);
      vector<int> metsig2012_sigtot(datasets[0].ROCmetsig2012.size(),0);
      vector<int> metsig2012_bkgtot(datasets[0].ROCmetsig2012.size(),0);

      for( vector<Dataset>::iterator data = datasets.begin(); data != datasets.end(); data++ ){
         if( data->isMC ){

            if( data->channel.compare("WJetsToLNu") == 0 ){
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

      TCanvas* cROC = new TCanvas("cROC","cROC",800,800);
      cROC->cd();

      gROCmet->SetTitle("ROC Curve;Background Efficiency;Signal Efficiency");

      gROCmet->SetMarkerStyle(20);
      gROCmet->SetMarkerColor(1);
      gROCmet->Draw("AP");

      gROCmetsig2011->SetMarkerStyle(20);
      gROCmetsig2011->SetMarkerColor(2);
      gROCmetsig2011->Draw("P");

      gROCmetsig2012->SetMarkerStyle(20);
      gROCmetsig2012->SetMarkerColor(3);
      gROCmetsig2012->Draw("P");

      TF1* fline = new TF1("fline", "x", 0, 1);
      fline->SetLineColor(1);
      fline->SetLineStyle(7);
      fline->Draw("same");

      TLegend *lROC = new TLegend(0.605528,0.655866,0.866834,0.816333);
      lROC->AddEntry(gROCmet,"met","p");
      lROC->AddEntry(gROCmetsig2011,"metsig2011","p");
      lROC->AddEntry(gROCmetsig2012,"metsig2012","p");
      lROC->Draw("same");

      cROC->Print("results/ROCplot.root");

      delete cROC;
      delete lROC;
   } // ROC plots

   return 0;
}
