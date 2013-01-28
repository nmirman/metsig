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


//code to perform met response corrections


TProfile* make_presp_qt( vector<event>* eventref, string title="presp_qt" ) {
//assumes qt, ut_par values are correctly filled

  TProfile* presp_qt = new TProfile (title.c_str(),
		    "Response = |<u_{#parallel}>|/q_{T} vs. q_{T};q_{T} (GeV);Response", 25, 0, 100);


	for( vector<event>::iterator ev = eventref->begin(); ev < eventref->end(); ev++ ){
		presp_qt->Fill( ev->qt, -(ev->ut_par)/(ev->qt), ev->weight );
	}

	return presp_qt;
}

TProfile* make_psig_nvert( vector<event>* eventref, string title="psig_nvert" ) {
	TProfile* p=new TProfile(title.c_str(),
	         "Significance vs. N Vertices;N Vertices;<S_{E}>", 30, 0, 30);

	for( vector<event>::iterator ev = eventref->begin(); ev < eventref->end(); ev++ ){
		p->Fill(ev->nvertices,ev->sig,ev->weight);
	}

	return p;
}

TProfile* make_psig_qt( vector<event>* eventref, string title="psig_qt" ) {
	TProfile* p=new TProfile(title.c_str(),
	         "Significance vs. q_{T};q_{T} (GeV);<S_{E}>", 15, 0, 100);

	for( vector<event>::iterator ev = eventref->begin(); ev < eventref->end(); ev++ ){
		p->Fill(ev->qt,ev->sig,ev->weight);
	}

	return p;
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

   //numevents=100000;

   bool use_data=true;
   bool use_mc=!use_data;


   // fill eventvecs
   if(use_mc) {
	   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Zmumu/"
			 "Zmumu_MC_DYJettoLL_TuneZ2_M-50_7TeV_madgraph_tauola_20121221.root",
			 eventvec_MC, numevents, true);
	   fitter.MatchMCjets( eventvec_MC );
   }

   if(use_data) {
	   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Zmumu/"
			 "Zmumu_data_DoubleMu_Run2011A_08Nov2011_v1_20121221.root",
			 eventvec_data, numevents/2, false);
	   fitter.ReadNtuple( "/eos/uscms/store/user/nmirman/Zmumu/"
			 "Zmumu_data_DoubleMu_Run2011B_19Nov2011_v1_20121221.root",
			 eventvec_data, numevents/2, false);
   }
   if(use_mc)	cout << "\n  MC EVENTS: " << eventvec_MC.size() << endl;
   if(use_data)	cout << "DATA EVENTS: " << eventvec_data.size() << endl;

   // minimize
   if(use_mc) {
	   cout << "\n ############################ " << endl;
	   cout << " ###########  MC  ########### " << endl;
	   cout << " ############################ \n" << endl;
	   fitter.RunMinimizer( eventvec_MC );
   }
   if(use_data) {
	   cout << "\n ############################ " << endl;
	   cout << " ########### Data ########### " << endl;
	   cout << " ############################ \n" << endl;
	   fitter.RunMinimizer( eventvec_data );
   }

   //fitter.PlotsDataMC( eventvec_data, eventvec_MC, "results/plotsDataMC.root" );

   //
   // ######################### END FIT #########################
   //


   //
   // ######################### MET RESPONSE CORRECTIONS #########################
   //

   std::cout << "\n######################### MET RESPONSE CORRECTIONS #########################\n" << std::endl;

   vector<event>* eventvec;
   if(use_data) {
	   eventvec=&eventvec_data;
   } else if(use_mc){
	   eventvec=&eventvec_MC;
   } else {
	   cout << "PROBLEM:  No events to use!" << endl;
	   return 0;
   }

   TProfile* presp_qt= make_presp_qt(eventvec, "presp_qt (original)");
   TProfile* psig_nvert= make_psig_nvert(eventvec, "psig_nvert (original)");
   TProfile* psig_qt= make_psig_qt(eventvec, "psig_qt (original)");

   TF1* func=new TF1("func","[0]+[1]*exp([2]*x)");
   func->SetParName(0,"Offset");
   func->SetParName(1,"Scale");
   func->SetParName(2,"Power");

   func->SetParameter(0,1);
   func->SetParameter(1,-0.4);
   func->SetParameter(2,-0.05);

   presp_qt->Fit(func);


   //vector<event>* eventvec_corr=new vector<event>();
   vector<event> eventvec_corr;

   //TH1D* hpt_mult=correct_pt(eventvec_MC, eventvec_corr, func);

   //actual correction

   TH1D* hpt_mult=new TH1D("hpt_mult","p_{T} Multiplier = 1/|f(q_{T})|;Multiplier",500,0,5);
   TH2D* hf_vs_ut=new TH2D("hf_vs_ut","u_{#parallel}/f(q_{T}) vs. q_{T};q_{T};u_{#parallel}/f(q_{T})", 25, 0, 100, 25, 0, 100 );

   for( vector<event>::iterator ev = eventvec->begin(); ev < eventvec->end(); ev++ ){

	   event* evtemp=new event(*ev);
	   //double pt_mult=func->Eval(ev->qt)/ev->ut_par;	//old pt_mult
	   //double pt_mult= abs( ev->qt/func->Eval(ev->qt) );
	   //double pt_mult=-ev->qt/ev->ut_par;
	   double pt_mult=abs( 1/func->Eval(ev->qt) );

	   hpt_mult->Fill(pt_mult);
	   hf_vs_ut->Fill(ev->qt,ev->ut_par/func->Eval(ev->qt));

	   // high pt jets
	   for(unsigned int i=0; i<evtemp->jet_ptUncor.size(); i++) {
		   evtemp->jet_ptUncor[i] *= pt_mult;
		   evtemp->jet_ptL123[i] *= pt_mult;
		   evtemp->jet_ptT1[i] *= pt_mult;
	   }

	   // pseudojet
	   evtemp->pjet_vectpt*=pt_mult;
	   evtemp->pjet_scalpt*=pt_mult;

	   eventvec_corr.push_back(*evtemp);
   }


   //fitter.FindSignificance(fitter.gMinuit->X(), eventvec_corr);

   fitter.RunMinimizer( eventvec_corr );

   TProfile* presp_qt_corr=make_presp_qt(&eventvec_corr, "presp_qt (corrected)");
   TProfile* psig_nvert_cor= make_psig_nvert(&eventvec_corr, "psig_nvert (corrected)");
   TProfile* psig_qt_cor= make_psig_qt(&eventvec_corr, "psig_qt (corrected)");

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
   
   presp_qt->Write();
   psig_nvert->Write();
   psig_qt->Write();
   presp_qt_corr->Write();
   psig_nvert_cor->Write();
   psig_qt_cor->Write();
   hpt_mult->Write();
   hf_vs_ut->Write();
   
   
   file->Write();
   file->Close();

   return 0;
}
