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
   const char* filename;
   const char* channel;
   bool isMC;
   Dataset( const char* f, const char* c, bool i) : filename(f), channel(c), isMC(i) {}
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
   while( (c = getopt(argc, argv, "n:j:hscob")) != -1 ) {
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

         case 'h' :
            cout << "Usage: ./DoFit <flags>\n";
            cout << "Flags: \n";
            cout << "\t-n <number>\t Fraction of events to fit.  Default at -1.\n";
            cout << "\t-j <number>\t Jet bin pt threshold.  Default at 20 GeV.\n";
            cout << "\t-s\t          'Short' run, 10%% of events.\n";
            cout << "\t-c\t          Apply response correction.\n";
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


   //
   // loop through data and mc, run fit, fill histograms
   //
   for(int i=0; i < 2; i++){

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

         fitter.ReadNtuple( data->filename, eventvec, numevents,
               data->isMC, data->channel, do_resp_correction );
      }

      cout << "\nDATASET SIZE: " << eventvec.size() << " EVENTS\n" << endl;

      // run fit
      fitter.RunMinimizer( eventvec );

      // fill histograms
      fitter.FillHists( eventvec, "Zmumu" ); 

   }
   // print histograms
   fitter.PrintHists( "results/plotsDataMC.root", "Zmumu" );

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
   string outfilename = "/fitresults.root";

   TFile *file = new TFile((pathstr+outfilename).c_str(), "RECREATE");
   file->cd();
   tree->Write();
   file->Write();
   file->Close();

   return 0;
}
