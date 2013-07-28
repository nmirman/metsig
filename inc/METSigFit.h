#ifndef METSIG_FIT_H
#define METSIG_FIT_H

#include <vector>
#include <cmath>
#include <map>
#include <TH1.h>
#include <TH2.h>

#include "TMath.h"
#include "TLorentzVector.h"
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"
using namespace std;

struct event {

   string channel;

   int nvertices;
   double weight;

   double met;
   double sig;
   double sig_init;
   double det;
   double cov_xx;
   double cov_xy;
   double cov_yy;

   double cov_xx_highpt;
   double cov_xx_pjet;
   double cov_dtt;
   double cov_dff;

   double qt;
   double ut;
   double ut_par;
   double ut_perp;

   double resp_correction;

   // pseudojet
   double pjet_phiL123;
   double pjet_vectptL123;
   double pjet_scalptL123;
   double pjet_phiT1;
   double pjet_vectptT1;
   double pjet_scalptT1;
   double pjet_size;

   // MC MET smearing
   double met_varx;
   double met_vary;
   double met_rho;

   // variables for ROC
   double metsig2011;

   // leptons
   //vector<TLorentzVector> muon_4vect;
   //vector<TLorentzVector> electron_4vect;
   vector<double> muon_pt;
   vector<double> muon_phi;
   vector<double> electron_pt;
   vector<double> electron_phi;

   // high pt jets
   vector<double> jet_phi;
   vector<double> jet_eta;
   vector<double> jet_ptUncor;
   vector<double> jet_ptL123;
   vector<double> jet_ptT1;
}; 

class Fitter{
   public:

      Fitter();
      ~Fitter();

      void ReadNtuple(const char[], vector<event>&, const double, const bool,
            string, const bool=false, const int=0, const int=0);
      //void MatchMCjets(vector<event>&);
      void RunMinimizer(vector<event>&);
      void FindSignificance(const double*, vector<event>&);
      void FillHists(vector<event>&, string);
      void PrintHists( const char*, string );
      void PJetReweight(vector<event>&, vector<event>&);
      void ResponseCorrection(vector<event>&);

      bool significance_cut;
      double jetbinpt, jetcorrpt;
      double psig_nvert_corr, psig_qt_corr, pchi2slope_left, pchi2slope_right;

      ROOT::Minuit2::Minuit2Minimizer* gMinuit;

      double dpt_(double x, double _eta){
         double feta = fabs(_eta);
         int ieta =feta<4.5? feta/0.5 : 9;
         double p0 = sigmaPt[ieta][0];
         double p1 = sigmaPt[ieta][1];
         double p2 = sigmaPt[ieta][2];
         double p3 = sigmaPt[ieta][3];
         return sqrt(((TMath::Sign(1.,p0)*p0*p0/x/x)+(p1*p1*pow(x,p3-1)))+p2*p2);
      }

      double dph_(double x, double _eta){
         double feta = fabs(_eta);
         int ieta =feta<4.5? feta/0.5 : 9;
         double p0 = sigmaPhi[ieta][0];
         double p1 = sigmaPhi[ieta][1];
         double p2 = sigmaPhi[ieta][2];
         double p3 = sigmaPhi[ieta][3];
         double p4 = sigmaPhi[ieta][4];
         return (sqrt((p0*p0/x/x+(p1*p1/x))+p2*p2)+(p3/x))   + ((p4/x)/sqrt(x));
      }

      // histograms
      map<string, TH1*> histsData_;
      map<string, TH1*> histsMC_;
      map<string, TH2*> profsData_;
      map<string, TH2*> profsMC_;

      map<string, TH1*> histsMC_signal_;
      map<string, TH1*> histsMC_top_;
      map<string, TH1*> histsMC_EWK_;
      map<string, TH1*> histsMC_QCD_;
      map<string, TH1*> histsMC_gamma_;
      map<string, TH1*> histsMC_DY_;

   private:
      double Min2LL(const double*);

      ROOT::Math::IMultiGenFunction* fFunc;

      vector<event>* eventvecPnt;

      //jet resolutions from 2010
      static const double sigmaPt[10][4];
      static const double sigmaPhi[10][5];
      static const double puWeights[50];


}; 
#endif
