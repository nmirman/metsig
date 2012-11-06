#ifndef METSIG_FIT_H
#define METSIG_FIT_H

#include <vector>
#include <cmath>

#include "TMath.h"
#include "TLorentzVector.h"
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"

struct event {

   int nvertices;

   double sig;
   double det;
   double cov_xx;
   double cov_xy;
   double cov_yy;

   double qt;
   double ut;
   double ut_par;

   std::vector<double> muon_pt;
   std::vector<double> muon_phi;
   std::vector<TLorentzVector> muon_4vect;

   // high pt jets
   std::vector<double> jet_pt;
   std::vector<double> jet_phi;
   std::vector<double> jet_eta;
   std::vector<double> jet_ptL123;
   std::vector<double> jet_ptT1;
   std::vector<int> jet_matchIndex;
   std::vector<double> jet_matchdR;

   // pseudojet
   double pjet_phi;
   double pjet_vectpt;
   double pjet_scalpt;
   double pjet_size;

   std::vector<double> genjet_pt;
   std::vector<double> genjet_phi;
   std::vector<double> genjet_eta;

};

class Fitter{
   public:

      Fitter();
      ~Fitter();

      void ReadNtuple(const char[], std::vector<event>&, const bool);
      void MatchMCjets(std::vector<event>&);
      void RunMinimizer(std::vector<event>&);
      void FindSignificance(const double*, std::vector<event>&);

      double jetfitLOW, jetfitHIGH, jetcorrMIN;

      const double *parameter;

   private:
      double Min2LL(const double*);

      ROOT::Minuit2::Minuit2Minimizer* gMinuit;
      ROOT::Math::IMultiGenFunction* fFunc;

      std::vector<event>* eventvecPnt;

      //jet resolutions from 2010
      static const double sigmaPt[10][4];
      static const double sigmaPhi[10][5];

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

}; 
#endif
