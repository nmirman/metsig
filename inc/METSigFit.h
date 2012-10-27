#ifndef METSIG_FIT_H
#define METSIG_FIT_H

#include <vector>
#include <cmath>

#include "TMath.h"

class Fitter{
   public:
      Fitter();
      ~Fitter();

      void ReadNtuple(const char[], bool);

   private:
      double Min2LL(const double*);
      void MatchMCjets();
/*
      //jet resolutions from 2010
      double sigmaPt[10][4]={{-0.349206, 0.297831, 0, 0.471121},
         {-0.499735, 0.336391, 0, 0.430689},
         {-0.561649, 0.420293, 0, 0.392398},
         {-1.12329, 0.657891, 0, 0.139595},
         {1.04792, 0.466763, 0, 0.193137},
         {2.56933, 0.305802, 0, 0.398929},
         {2.81978, 0.272373, 0, 0.579396},
         {1.65966, 0.223683, 0, 0.60873},
         {1.41584, 0.209477, 0, 0.588872},
         {1.41584, 0.209477, 0, 0.588872}};

      double sigmaPhi[10][5] ={   {    926.978,      2.52747,    0.0304001,  -926.224,     -1.94117},
         {3.32512e-06,     0.063941,  -0.00387593,  0.301932,    -0.825352},
         {    0.38469,    0.0755727,   -0.0044353,  0.453887,      -1.8947},
         {2.92001e-07,    0.0718389,  -0.00385579,  0.403668,     -0.62698},
         { 0.00336639,    0.0880209,   -0.0023084,  0.214304,    -0.416353},
         {    11.1957,     0.643236,   0.00711422,  -10.7613,     0.280927},
         {     1.9027, -4.56986e-06,    0.0304793,  -1.09124,    -0.136402},
         {    2.11446,     0.203329,   -0.0175832,  -1.67946,  -0.00853474},
         {   0.765787, -3.90638e-06, -4.70224e-08,   0.11831,      -1.4675},
         {    259.189,   0.00132792,    -0.311411,  -258.647,            0}};

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
*/
      struct event {

         int nvertices;

         double sig;
         double qt;
         double ut;
         double ut_par;
         
         std::vector<double> muon_pt;
         std::vector<double> muon_phi;

         // high pt jets
         std::vector<double> jet_pt;
         std::vector<double> jet_phi;
         std::vector<double> jet_eta;
         std::vector<double> jet_corrL123;
         std::vector<double> jet_corrL1;
         std::vector<int> jet_matchIndex;
         std::vector<double> jet_matchdR;

         // pseudojet
         double pjet_pt; // vector sum of unclustered energy
         double pjet_phi;
         double pjet_ht; // scalar sum of unclustered energy

         std::vector<double> genjet_pt;
         std::vector<double> genjet_phi;
         std::vector<double> genjet_eta;

      };

      std::vector<event> eventvec_MC;
      std::vector<event> eventvec_data;
     
};
#endif
