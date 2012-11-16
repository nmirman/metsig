#ifndef METSIG_FIT_H
#define METSIG_FIT_H

#include <vector>
#include <cmath>

#include "TMath.h"
#include "TLorentzVector.h"
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"
using namespace std;

struct event {

   int nvertices;
   double weight;

   double met;
   double sig;
   double det;
   double cov_xx;
   double cov_xy;
   double cov_yy;

   double qt;
   double ut;
   double ut_par;
   double ut_perp;

   vector<double> muon_pt;
   vector<double> muon_phi;
   vector<TLorentzVector> muon_4vect;

   // high pt jets
   vector<double> jet_phi;
   vector<double> jet_eta;
   vector<double> jet_ptUncor;
   vector<double> jet_ptL123;
   vector<double> jet_ptT1;
   vector<int> jet_matchIndex;
   vector<double> jet_matchdR;
   vector<bool> jet_id;

   // pseudojet
   double pjet_phi;
   double pjet_vectpt;
   double pjet_scalpt;
   double pjet_size;

   vector<double> genjet_pt;
   vector<double> genjet_phi;
   vector<double> genjet_eta;

};

class Fitter{
   public:

      Fitter();
      ~Fitter();

      void ReadNtuple(const char[], vector<event>&, const int, const bool);
      void MatchMCjets(vector<event>&);
      void GetPUWeights(vector<event>&, vector<event>&);
      void RunMinimizer(vector<event>&);
      void FindSignificance(const double*, vector<event>&);
      void PlotsDataMC(vector<event>&, vector<event>&, const char[]);

      double jetfitLOW, jetfitHIGH, jetcorrMIN;

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
