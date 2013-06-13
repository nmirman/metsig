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
   vector<event> eventvec_MC;
   vector<event> eventvec_data;
   bool do_resp_correction=false;

   // option flags
   char c;
   double numevents = -1;
   while( (c = getopt(argc, argv, "n:j:hscob")) != -1 ) {
      switch(c)
      {
         case 'n' :
            numevents = atoi(optarg);
            break;

         case 'j' :
            fitter.jetbinpt = atoi(optarg);
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
            cout << "\t-s\t          'Short' run, 10\% of events.\n";
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



   // fill eventvecs

   // mc
   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130501/DYJetsToLL.root",
         eventvec_MC, numevents, true, "DYJetsToLL", do_resp_correction);

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

   fitter.MatchMCjets( eventvec_MC );

   // data
   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Ntuples/Zmumu/20130427/Data.root",
         eventvec_data, numevents, false, "Data", false);

   cout << "\n  MC EVENTS: " << eventvec_MC.size() << endl;
   cout << "DATA EVENTS: " << eventvec_data.size() << endl;

   cout << "\n ############################ " << endl;
   cout << " ###########  MC  ########### " << endl;
   cout << " ############################ \n" << endl;
   fitter.RunMinimizer( eventvec_MC );


   cout << "\n ############################ " << endl;
   cout << " ########### Data ########### " << endl;
   cout << " ############################ \n" << endl;
   fitter.RunMinimizer( eventvec_data );

   fitter.FillHists( eventvec_data, "Zmumu" );
   fitter.FillHists( eventvec_MC, "Zmumu" );
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
