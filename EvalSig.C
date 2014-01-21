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
#include "TSystem.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
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
   string dirname;
   string process;
   bool isMC;
   int size;

   vector<string> filenames;

   // for ROC curve
   vector<ROCPoint> ROCmet;
   vector<ROCPoint> ROCmetsig2011;
   vector<ROCPoint> ROCmetsig2012;
   vector<ROCPoint> ROCmetrht;

   // constructor
   Dataset( string f="", string p="", bool i=0) : dirname(f), process(p), isMC(i) {}
};

int main(int argc, char* argv[]){

   freopen ("stderr.txt","w",stderr);

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

   string xrdopt = "path";
   string file_catalog = "file_catalog.txt";

   while( (c = getopt(argc, argv, "n:p:o:t:b:v:q:fhscbmrdw")) != -1 ) {
      switch(c)
      {
         case 'n' :
            fracevents = atof(optarg);
            break;

         case 's' :
            fracevents = 0.1;
            break;

         case 'c' :
            xrdopt = "cache";
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

         case 'q':
            file_catalog = optarg;
            break;

         case 'h' :
            cout << "Usage: ./EvalSig <flags>\n";
            cout << "Flags: \n";
            cout << "\t-n <number>\t  Fraction of events to fit.  Default at -1.\n";
            cout << "\t-j <number>\t  Jet bin pt threshold.  Default at 20 GeV.\n";
            cout << "\t-s\t          'Short' run, 10% of events.\n";
            cout << "\t-c\t          Read from xrootd cache.\n";
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
            cout << "\t-q\t          Filename containing ntuple filenames.\n";
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


   ifstream inFile;
   inFile.open(file_catalog.c_str());

   while(!inFile.eof()){

      Dataset data;
      data.path = "root://osg-se.cac.cornell.edu//xrootd/"
         +xrdopt+"/cms/store/user/nmirman/Ntuples/METsig";

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
      if( data.dirname.find("Run") != string::npos ){
         data.isMC = false;
         data.process = "Data";
      }

      string date = "20130830";
      if( channel.compare("Zmumu") == 0 and data.isMC ){
         date = "20130913";
      }
      if( channel.compare("Wenu") == 0 ){
         date = "20130916";
      }
      if( channel.compare("Dijet") == 0 and data.isMC ){
         date = "20130913";
      }
      if( channel.compare("Ttbar0lept") == 0 ){
         date = "20130913";
      }

      // vector of filenames
      string file;
      while( stream >> file ){
         data.filenames.push_back( file );
      }

      // add to datasets
      if( channel.compare(data.channel) == 0 and date.compare(data.date) == 0 )
         if( (run_data and !data.isMC) or (run_mc and data.isMC) )
            datasets.push_back( data );

   }

   // get size of datasets
   cout << endl;
   for( vector<Dataset>::iterator data = datasets.begin(); data != datasets.end(); data++ ){
      TChain tree("events");
      for( vector<string>::iterator file = data->filenames.begin();
            file != data->filenames.end(); file++){
         TString fn = data->path+"/"+data->channel+"/"+data->date+"/"+data->dirname+"/"+(*file);
         tree.Add( fn );
      }
      data->size = tree.GetEntries();
      cout << data->channel << " " << data->date << " " << data->dirname
         << ": " << data->size << " events" << endl;
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

      int section_size = 1000000;
      int num_events = fracevents*data->size;
      int num_sections = 1 + ((num_events-1)/section_size);
      cout << "Opening dataset " << data->dirname << endl;
      cout << "Divide " << num_events << " events into " << num_sections << " sections..." << endl;

      for(int isec=0; isec < num_sections; isec++){
         int start = isec*section_size;
         int end = (isec == num_sections-1) ? num_events : start + section_size;
         cout << "Begin section [" << start << ", " << end << "]" << endl;

         double jec_varmc = data->isMC ? jec_var : 0;

         vector<event> eventvec;
         string fullname = data->path+"/"+data->channel+"/"+data->date+"/"+data->dirname;
         fitter.ReadNtuple( fullname.c_str(), data->filenames, eventvec, 1,
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
