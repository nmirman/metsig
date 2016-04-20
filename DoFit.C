#include "METSigFit.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
using namespace std;

struct Dataset {
   // dataset info
   string path;
   string date;
   string channel;
   string dirname;
   string process;
   bool isMC;
   int size;

   vector<string> filenames;

   // constructor
   Dataset( string f="", string p="", bool i=0) : dirname(f), process(p), isMC(i) {}
};

int main(int argc, char* argv[]){
   
   freopen ("stderr.txt","w",stderr);

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


   ifstream inFile;
   inFile.open("file_catalog.txt"/*file_catalog.c_str()*/);

   while(!inFile.eof()){

      Dataset data;
      string xrdopt = "cache";
      data.path = "Ntuples/";
      string channel = "Zmumu";

      // get line from file
      string line;
      getline(inFile,line);
      stringstream stream(line);

      // get ntuple attributes
      stream >> data.channel;
      stream >> data.date;
      stream >> data.dirname;

      data.process = data.dirname;
      data.isMC = true;
      if( data.dirname.find("Data") != string::npos ){
         data.isMC = false;
         data.process = "Data";
      }

      string date = "20160418";

      // vector of filenames
      string file;
      while( stream >> file ){
         data.filenames.push_back( file );
      }

      // add to datasets
      if( channel.compare(data.channel) == 0 and date.compare(data.date) == 0 )
         datasets.push_back( data );

   }

   // get size of datasets
   cout << endl;
   for( vector<Dataset>::iterator data = datasets.begin(); data != datasets.end(); data++ ){
      TChain chain("events");
      for( vector<string>::iterator file = data->filenames.begin();
            file != data->filenames.end(); file++){
         TString fn = data->path+"/"+data->channel+"/"+data->date+"/"+data->dirname+"/"+(*file);
         chain.Add( fn );
      }
      data->size = chain.GetEntries();
      cout << data->channel << " " << data->date << " " << data->dirname
         << ": " << data->size << " events" << endl;
   }
   cout << endl;


   //
   // loop through data and mc, run fit, fill histograms
   //
   double pileup_dist [2][100] = {};
   double mupt_dist [2][200] = {};
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

         string fullname = data->path+"/"+data->channel+"/"+data->date+"/"+data->dirname;
         //string xrdname = fullname;
         //xrdname.replace(xrdname.begin(),xrdname.begin()+11,
         //      "root://osg-se.cac.cornell.edu//xrootd/path/cms/store");
         fitter.ReadNtuple( fullname.c_str(), data->filenames, eventvec, numevents,
               data->isMC, data->process, do_resp_correction, -1/*start*/, -1/*end*/, 0/*jec_varmc*/ );
      }

      cout << "\nDATASET SIZE: " << eventvec.size() << " EVENTS\n" << endl;

      // pile up reweighting
      for( vector<event>::iterator ev = eventvec.begin(); ev != eventvec.end(); ev++ ){
         if( ev->nvertices < 100 ) pileup_dist[i][ev->nvertices] += ev->weight;
      }
      // normalize
      double norm = 0.0;
      for(int n=0; n < 100; n++) norm += pileup_dist[i][n];
      for(int n=0; n < 100; n++) pileup_dist[i][n] /= norm;
      if( i==1 ){ // MC
         // find minimum weight
         double minx = 100;
         for( int n=0; n < 100; n++){
            if( pileup_dist[1][n] != 0 ){
               double x = pileup_dist[0][n] / pileup_dist[1][n];
               if( x < minx and x > 0 ) minx = x;
            }
         }
         for( vector<event>::iterator ev = eventvec.begin(); ev != eventvec.end(); ev++ ){
            if( pileup_dist[0][ev->nvertices] != 0 ){
               ev->weight *= pileup_dist[0][ev->nvertices] / pileup_dist[1][ev->nvertices];
            }else{
               ev->weight *= minx;
            }
         }
      }
      
      // muon pt reweighting
      /*
      for( vector<event>::iterator ev = eventvec.begin(); ev != eventvec.end(); ev++ ){
         for( int m=0; m < 2; m++ ){
            int ipt = floor(ev->lepton_pt[m]);
            if( ipt > 199 ) ipt = 199;
            mupt_dist[i][ipt] += ev->weight;
         }
      }
      // normalize
      norm = 0.0;
      for(int n=0; n < 200; n++) norm += mupt_dist[i][n];
      for(int n=0; n < 200; n++) mupt_dist[i][n] /= norm;
      if( i==1 ){ // MC
         for( vector<event>::iterator ev = eventvec.begin(); ev != eventvec.end(); ev++ ){
               for( int m=0; m < 2; m++ ){
                  int ipt = floor(ev->lepton_pt[m]);
                  if( ipt > 199 ) ipt = 199;
                  ev->weight *= mupt_dist[0][ipt] / mupt_dist[1][ipt];
               }
         }
      }
      */

      // set type of MET
      fitter.met_type = met_type;

      // run fit
      if( true or i==0 ){
         fitter.RunMinimizer( eventvec );
      }else{
         fitter.FindSignificance(fitter.xmin, eventvec);
      }

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
