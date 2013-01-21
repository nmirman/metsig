#include "METSigFit.h"

#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include <iostream>
#include <unistd.h>
using namespace std;


//code to perform met response corrections


TH2D* make_resp_plot( vector<event>& eventref_MC ) {

	TH2D* presp_qt = new TH2D("presp_qt_Data",
			"Response = |<u_{#parallel}>|/q_{T} vs. q_{T};q_{T} (GeV);Response", 25, 0, 100, 100, -100, 100);

	for( vector<event>::iterator ev = eventref_MC.begin(); ev < eventref_MC.end(); ev++ ){
		presp_qt->Fill( ev->qt, -(ev->ut_par)/(ev->qt), ev->weight );
	}

	return presp_qt;
}

void correct_pt( vector<event>& event_source, vector<event>& event_dest, TF1* func ) {

	for( vector<event>::iterator ev = event_source.begin(); ev < event_source.end(); ev++ ){

		event* evtemp=new event(*ev);
		double pt_mult=func->eval(ev->qt)/ev->ut_perp;

		// high pt jets
		evtemp->jet_ptUncor*=pt_mult;
		evtemp->jet_ptL123*=pt_mult;
		evtemp->jet_ptT1*=pt_mult;

		// pseudojet
		evtemp->pjet_vectpt*=pt_mult;
		evtemp->pjet_scalpt*=pt_mult;

		event_dest.push_back(*evtemp);
	}
}


void met_resp( vector<event>& event_source ) {

	TH2D* presp_qt= make_resp_plot(event_source);
	TF1* func=new TF1("func","[0]+[1]*exp([2]*x)");
	presp_qt->Fit(func);

	vector<event> event_corr;
	correct_pt(event_source, event_corr, func);

	TH2D* presp_qt_corr=make_resp_plot(event_corr);

}



int main(int argc, char* argv[]){

   // setup fit results tree
   double jetbinpt=0, psig_nvert_corr=0, psig_qt_corr=0,
          pchi2slope_left=0, pchi2slope_right=0,
          a1=0, a2=0, a3=0, a4=0, a5=0, N1=0, S1=0, P1=0;
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
   tree->Branch("P1", &P1);

   // declarations
   Fitter fitter;
   vector<event> eventvec_MC;
   vector<event> eventvec_data;

   // option flags
   char c;
   int numevents = -1;
   while( (c = getopt(argc, argv, "n:j:h")) != -1 ) {
      switch(c)
      {
         case 'n' :
            numevents = atoi(optarg);
            break;

         case 'j' :
            fitter.jetbinpt = atoi(optarg);
            break;

         case 'h' :
            cout << "Usage: ./DoFit <flags>\n";
            cout << "Flags: \n";
            cout << "\t-n number\t Number of events to fit.  Default at -1.";
            cout << "\t-j number\t Jet bin pt threshold.  Default at 30 GeV.";
            cout << "\t-h\t Display this menu.";
            return -1;
            break;

         default :
            continue;
      }
   }

   //
   // ######################### BEGIN FIT #########################
   //

   bool use_data=false;


   // fill eventvecs
   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Zmumu/"
         "Zmumu_MC_DYJettoLL_TuneZ2_M-50_7TeV_madgraph_tauola_20121221.root",
         eventvec_MC, numevents, true);
   fitter.MatchMCjets( eventvec_MC );

   //don't use data
   /*fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Zmumu/"
         "Zmumu_data_DoubleMu_Run2011A_08Nov2011_v1_20121221.root",
         eventvec_data, numevents/2, false);
   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Zmumu/"
         "Zmumu_data_DoubleMu_Run2011B_19Nov2011_v1_20121221.root",
         eventvec_data, numevents/2, false);*/
   
   cout << "\n  MC EVENTS: " << eventvec_MC.size() << endl;
   /*
    * cout << "DATA EVENTS: " << eventvec_data.size() << endl;
    */

   // minimize
   cout << "\n ############################ " << endl;
   cout << " ###########  MC  ########### " << endl;
   cout << " ############################ \n" << endl;
   fitter.RunMinimizer( eventvec_MC );

   /*cout << "\n ############################ " << endl;
   cout << " ########### Data ########### " << endl;
   cout << " ############################ \n" << endl;
   fitter.RunMinimizer( eventvec_data );

   fitter.PlotsDataMC( eventvec_data, eventvec_MC, "results/plotsDataMC.root" );*/

   //
   // ######################### END FIT #########################
   //


   //
   // ######################### MET RESPONSE CORRECTIONS #########################
   //

   TH2D* presp_qt= make_resp_plot(eventvec_MC);
   TF1* func=new TF1("func","[0]+[1]*exp([2]*x)");
   presp_qt->Fit(func);

   vector<event> eventvec_corr;
   correct_pt(eventvec_MC, eventvec_corr, func);

   TH2D* presp_qt_corr=make_resp_plot(eventvec_corr);

   fitter.RunMinimizer( eventvec_corr );

   //
   // ######################### END MET RESPONSE CORRECTIONS #########################
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
   P1 = par[7];
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
