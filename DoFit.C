#include "METSigFit.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include <iostream>
#include <unistd.h>
using namespace std;

struct Dataset {
   // dataset info
   string path;
   string date;
   string channel;
   string filename;
   string process;
   bool isMC;
   int size;

   // constructor
   Dataset( string f, string p, bool i) : filename(f), process(p), isMC(i) {}
};

int main(int argc, char* argv[]){

   // setup fit results tree
   double jetbinpt=0, psig_nvert_corr=0, psig_qt_corr=0,
          pchi2slope_left=0, pchi2slope_right=0,
          a1=0, a2=0, a3=0, a4=0, a5=0, N1=0, S1=0;
   int fitStatus=-1;

   TTree *tree = new TTree("FitResults", "FitResults");
   tree->Branch("fitStatus", &fitStatus);
   tree->Branch("jetbinpt", &jetbinpt);
   tree->Branch("psig_vert_corr", &psig_nvert_corr);
   tree->Branch("psig_qt_corr", &psig_qt_corr);
   tree->Branch("pchi2slope_left", &pchi2slope_left);
   tree->Branch("pchi2slope_right", &pchi2slope_right);
   tree->Branch("a1", &a1);
   tree->Branch("a2", &a2);
   tree->Branch("a3", &a3);
   tree->Branch("a4", &a4);
   tree->Branch("a5", &a5);
   tree->Branch("N1", &N1);
   tree->Branch("S1", &S1);

   // declarations
   Fitter fitter;
   vector<Dataset> datasets;

   // option flags
   char c;
   double numevents = 1;
   bool do_resp_correction=false;
   int met_type = 4;
   int jec_var = 0;

   string file_results = "/fitresults.root";
   string file_plots = "results/plotsDataMC.root";
   while( (c = getopt(argc, argv, "n:j:m:r:p:v:hscob")) != -1 ) {
      switch(c)
      {
         case 'n' :
            numevents = atof(optarg);
            break;

         case 'j' :
            fitter.jetbinpt = atof(optarg);
            break;

         case 's' :
            numevents = 0.1;
            break;

         case 'c' :
            do_resp_correction=true;
            break;

         case 'm' :
            met_type = atoi(optarg);
            break;

         case 'r' :
            file_results = optarg;
            break;

         case 'p' :
            file_plots = optarg;
            break;

         case 'v':
            jec_var = atoi(optarg);
            break;

         case 'h' :
            cout << "Usage: ./DoFit <flags>\n";
            cout << "Flags: \n";
            cout << "\t-n <number>\t Fraction of events to fit.  Default at -1.\n";
            cout << "\t-j <number>\t Jet bin pt threshold.  Default at 20 GeV.\n";
            cout << "\t-m <number>\t Type of MET to use.  Default at -1.\n";
            cout << "\t-s\t          'Short' run, 10%% of events.\n";
            cout << "\t-c\t          Apply response correction.\n";
            cout << "\t-r <string>\t Output fit results to file.\n";
            cout << "\t-p <string>\t Output plots to file.\n";
            cout << "\t-v\t          Scale up (1) or down (-1) by JEC uncertainty.\n";
            cout << "\t-h\t          Display this menu.\n";
            return -1;
            break;

         default :
            continue;
      }
   }

   //
   // ######################### BEGIN FIT #########################
   //

   //
   // get all ntuples
   //

   // data
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

   // mc
   datasets.push_back( Dataset("DYJetsToLL/*.root", "DYJetsToLL", true) );
   datasets.push_back( Dataset("TTJets/*.root", "TTJets", true) );
   datasets.push_back( Dataset("Tbar_tW-channel/*.root", "Tbar_tW", true) );
   datasets.push_back( Dataset("T_tW-channel/*.root", "T_tW", true) );
   datasets.push_back( Dataset("WJetsToLNu/*.root", "WJetsToLNu", true) );
   datasets.push_back( Dataset("WW/*.root", "WW", true) );
   datasets.push_back( Dataset("WZ/*.root", "WZ", true) );
   datasets.push_back( Dataset("ZZ/*.root", "ZZ", true) );

   for( vector<Dataset>::iterator data = datasets.begin(); data != datasets.end(); data++ ){
      data->path = "/mnt/xrootd/user/nmirman/Ntuples/METsig";
      data->date = "20130830";
      data->channel = "Zmumu";
      if( data->isMC ) data->date = "20130913";
   }

   //
   // loop through data and mc, run fit, fill histograms
   //
   for(int i=1; i < 2; i++){

      if(i==0){
         cout << "\n ############################ " << endl;
         cout << " ########### Data ########### " << endl;
         cout << " ############################ \n" << endl;

      } else if(i==1){
         cout << "\n ############################ " << endl;
         cout << " ###########  MC  ########### " << endl;
         cout << " ############################ \n" << endl;
      }

      vector<event> eventvec;

      // read ntuples
      for( vector<Dataset>::iterator data = datasets.begin(); data != datasets.end(); data++ ){
         if( i==0 and data->isMC ) continue;
         if( i==1 and !(data->isMC) ) continue;

         string fullname = data->path+"/"+data->channel+"/"+data->date+"/"+data->filename;
         fitter.ReadNtuple( fullname.c_str(), eventvec, numevents,
               data->isMC, data->process, do_resp_correction, -1, -1, jec_var );
      }

      cout << "\nDATASET SIZE: " << eventvec.size() << " EVENTS\n" << endl;

      // set type of MET
      fitter.met_type = met_type;

      // run fit
      fitter.RunMinimizer( eventvec );

      // fill histograms
      fitter.FillHists( eventvec, "Zmumu" ); 

   }
   // print histograms
   fitter.PrintHists( file_plots.c_str(), "Zmumu" );

   //
   // ######################### END FIT #########################
   //


   // fill tree with fit results
   fitStatus = fitter.gMinuit->Status();
   jetbinpt = fitter.jetbinpt;
   const double *par = fitter.gMinuit->X();
   a1 = par[0];
   a2 = par[1];
   a3 = par[2];
   a4 = par[3];
   a5 = par[4];
   N1 = par[5];
   S1 = par[6];
   psig_nvert_corr = fitter.psig_nvert_corr;
   psig_qt_corr = fitter.psig_qt_corr;
   pchi2slope_left = fitter.pchi2slope_left;
   pchi2slope_right = fitter.pchi2slope_right;

   tree->Fill();

   // set up output file path
   std::string pathstr;
   char* path = std::getenv("WORKING_DIR");
   if (path==NULL) {
      pathstr = "./results";
   }else {
      pathstr = path;
   }

   TFile *file = new TFile((pathstr+file_results).c_str(), "RECREATE");
   file->cd();
   tree->Write();
   file->Write();
   file->Close();

   return 0;
}
