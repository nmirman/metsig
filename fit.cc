#include "TRandom3.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TROOT.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

const float pi = TMath::Pi();


float met;
int nv;
int nj;
float jpt[1000];
float jph[1000];
float jta[1000];
int nm;
float mpt[100];
float mph[100];
float pfj_l1[1000];
float pfj_l1l2l3[1000];
int gennj;
float genpt[1000];
float genph[1000];
float genta[1000];

TTree *t;
TTree *ttt;
int n;
float ppt[1000];
float phi[1000];
float eta[1000];
float jetcorrL123[1000];
float jetcorrL1[1000];
float genppt[1000];
float genphi[1000];
float geneta[1000];

float lpt;
float lph;
float lst;


TH1D *hpc2 = new TH1D("hpc2", "", 100, 0, 1);
TProfile *psig = new TProfile("psig", "psig", 15, 0, 15);
TH1D *hqt = new TH1D("hqt", ";q_{T} Z (GeV);Events/5 GeV", 30, 0, 150);
TProfile *pqt = new TProfile("pqt", "pqt", 15, 0, 100);
TH1D *hut_par = new TH1D("hut_par", ";u_{#parallel} (GeV);Events/10 GeV", 17, -150, 20);
TH1D *hut_perp = new TH1D("hut_perp", ";u_{#perp} (GeV);Events/4 GeV", 25, -50, 50);
TProfile *put = new TProfile("put", ";q_{T} (GeV);|<u_{#parallel}>|/q_{T}", 25, 0, 100);
bool FillTree=false;
bool FillHist=false;
bool FirstFit=true;

double Min2LL(const double *x);

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



void fit(){

   gROOT->Reset();

   TFile *f = new TFile("Zmumu_ntuple_20120919.root");
   TTree *tt = (TTree*)f->Get("events");
   tt->SetBranchAddress("v_size", &nv);
   tt->SetBranchAddress("met_et",&met);

   tt->SetBranchAddress("pfj_size", &nj);
   tt->SetBranchAddress("pfj_pt", jpt);
   tt->SetBranchAddress("pfj_phi", jph);
   tt->SetBranchAddress("pfj_eta", jta);

   tt->SetBranchAddress("mu_size", &nm);
   tt->SetBranchAddress("mu_pt", mpt);
   tt->SetBranchAddress("mu_phi", mph);

   tt->SetBranchAddress("pfj_l1", &pfj_l1);
   tt->SetBranchAddress("pfj_l1l2l3", &pfj_l1l2l3);

   tt->SetBranchAddress("genj_size", &gennj);
   tt->SetBranchAddress("genj_pt", genpt);
   tt->SetBranchAddress("genj_phi", genph);
   tt->SetBranchAddress("genj_eta", genta);

   t = new TTree("fit", "fit");
   t->SetDirectory(0);
   t->Branch("nv", &nv, "nv/I");
   t->Branch("nmu", &nm, "nmu/I");
   t->Branch("mpt", mpt, "mpt[nmu]/F");
   t->Branch("mph", mph, "mph[nmu]/F");
   t->Branch("n", &n, "n/I");
   t->Branch("ppt", ppt, "ppt[n]/F");
   t->Branch("phi", phi, "phi[n]/F");
   t->Branch("eta", eta, "eta[n]/F");
   t->Branch("jetcorrL123", jetcorrL123, "jetcorrL123[n]/F");
   t->Branch("jetcorrL1", jetcorrL1, "jetcorrL1[n]/F");
   t->Branch("lpt", &lpt, "lpt/F");
   t->Branch("lph", &lph, "lph/F");
   t->Branch("lst", &lst, "lst/F");
   t->Branch("genppt", genppt, "genppt[n]/F");
   t->Branch("genphi", genphi, "genphi[n]/F");
   t->Branch("geneta", geneta, "geneta[n]/F");

   ttt=t->CloneTree(0);
   ttt->SetDirectory(0);

   for(Long64_t ev=0; ev<tt->GetEntries(); ++ev){
      tt->GetEntry(ev);

      if(nm!=2 || mpt[0]<25  || mpt[1]<20) continue;
      n=0;
      float lpx=0;
      float lpy=0;
      lst=0;
      for(int i=0; i<nj; ++i){

         jetcorrL123[i] = pfj_l1l2l3[i];
         jetcorrL1[i] = pfj_l1[i];

         if(jpt[i]>3.0){
            ppt[n]=jpt[i];
            phi[n]=jph[i];
            eta[n]=jta[i];
            ++n;
         }
         else{
            lpx+=jpt[i]*cos(jph[i]);
            lpy+=jpt[i]*sin(jph[i]);
            lst+=jpt[i];
         }

      }
      lpt = sqrt(lpx*lpx+lpy*lpy);
      lph = atan2(lpy, lpx);
      t->Fill();
   }

   ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
   min->SetTolerance(0.001);
   min->SetStrategy(0);
   min->SetPrintLevel(3);
   ROOT::Math::Functor ff(&Min2LL, 12); 
   min->SetFunction(ff);
   min->SetVariable(0, "a1", 1.5, 0.01);
   min->SetVariable(1, "a2", 1.5, 0.01);
   min->SetVariable(2, "a3", 1.5, 0.01);
   min->SetVariable(3, "a4", 1.5, 0.01);
   min->SetVariable(4, "a5", 1.5, 0.01);
   min->SetVariable(5, "k0", 1, 0.01);
   min->SetVariable(6, "k1", 1, 0.01);
   min->SetVariable(7, "k2", 1, 0.01);
   min->SetVariable(8, "N1", 4, 0.01);
   min->SetVariable(9, "S1", 0.5, 0.01);
   min->SetVariable(10,"N2", 4, 0.01);
   min->SetVariable(11,"S2", 0.5, 0.01);
   //fit, first round
   min->Minimize();				
   FillTree=true;
   //fill the tree, remove the peak, re-fit only gaussian part
   Min2LL(min->X());
   FillTree=false;
   FirstFit=false;
   min->SetStrategy(1);
   min->Minimize();
   min->Hesse();
   //Fill the histogram (tails are included now).
   FillHist=true;
   Min2LL(min->X());

   TCanvas *canvas = new TCanvas("canvas","canvas",900,900);
   canvas->Divide(2,2);

   canvas->cd(1);
   hqt->Draw();

   canvas->cd(2);
   hut_par->Draw();

   canvas->cd(3);
   hut_perp->Draw();

   canvas->cd(4);
   put->Draw();

   TCanvas *canvas2 = new TCanvas();
   hpc2->SetMinimum(0);
   hpc2->Draw("E");
   hpc2->GetXaxis()->SetTitle("p(#chi^{2})");

   TF1 *f1 = new TF1("f1","pol0",0,1);
   hpc2->Fit("f1");
   hpc2->Draw("E");

   std::cout<< "CHI2 = " << f1->GetChisquare() << "/" << f1->GetNDF() << 
      ", " << f1->GetProb() << std::endl;

   TCanvas *canvas3 = new TCanvas();
   psig->SetMinimum(0);
   psig->SetMaximum(4);
   psig->Draw();
   psig->GetXaxis()->SetTitle("n vertices");
   psig->GetYaxis()->SetTitle("average S_{E}");

   TF1 *f2 = new TF1("f2","pol0",0,15);
   psig->Fit("f2");
   psig->Draw("E");

   TCanvas *canvas4 = new TCanvas();
   put->SetLineWidth(2);
   put->Draw();

   TF1 *f3 = new TF1("f3","pol0",0,150);
   f3->SetParameter(0,1.0);
   f3->SetLineColor(2);
   f3->Draw("same");

   std::cout<< "CHI2 = " << f2->GetChisquare() << "/" << f2->GetNDF() << 
      ", " << f2->GetProb() << std::endl;

}

double Min2LL(const double *x){
   double m2ll=0;
   int nevets = FirstFit||FillHist ? t->GetEntries() : ttt->GetEntries();

   for(Long64_t ev=0; ev<nevets; ++ev){

      if(FirstFit||FillHist) t->GetEntry(ev);
      else ttt->GetEntry(ev);

      double mex=0;
      double mey=0;
      double cxx=0;
      double cxy=0;
      double cyy=0;

      //implement jet energy correction for Type-I MET
      std::vector<float> pptcorr;
      for(int i=0; i<n; i++){
         pptcorr.push_back( ppt[i]*(jetcorrL123[i] + 1 - jetcorrL1[i]) );
      }

      //jets 6 GeV and above, scale 2010 resolutions
      for(int i=0; i<n; ++i){
         float feta = fabs(eta[i]);
         double c = cos(phi[i]);
         double s = sin(phi[i]);

         mex-=pptcorr[i]*c;
         mey-=pptcorr[i]*s;

         double dpt=0;
         double dph=0;

         if(ppt[i]>10){
            int index=-1;
            if(feta<0.5) index=0;
            else if(feta<1.1) index=1;
            else if(feta<1.7) index=2;
            else if(feta<2.3) index=3;
            else index=4;

            dpt = x[index]*ppt[i]*dpt_(ppt[i], eta[i]);
            dph = ppt[i]*dph_(ppt[i], eta[i]);
         }
         else{
            int index=-1;
            if(feta<2.4) index=0;
            else if(feta<3) index=1;
            else index=2;
            dpt = x[5+index]*sqrt(ppt[i]);
            dph = 0;
         }

         double dtt = dpt*dpt;
         double dff = dph*dph;
         cxx += dtt*c*c+ dff*s*s;
         cxy += c*s*(dtt-dff);
         cyy += dff*c*c+ dtt*s*s;
      }

      //muons with 0 resolutions
      mex-=mpt[0]*cos(mph[0]);
      mey-=mpt[0]*sin(mph[0]);
      mex-=mpt[1]*cos(mph[1]);
      mey-=mpt[1]*sin(mph[1]);

      //unclustered energy, parametrize by e.g.sumEt
      mex-=lpt*cos(lph);
      mey-=lpt*sin(lph);

      double ctt = x[8]*x[8]+x[9]*x[9]*lst;
      double cff = x[10]*x[10]+x[11]*x[11]*lst;

      cxx+=ctt*cos(lph)*cos(lph)+cff*sin(lph)*sin(lph);
      cxy+=cos(lph)*sin(lph)*(ctt-cff);
      cyy+=cff*cos(lph)*cos(lph)+ctt*sin(lph)*sin(lph);

      double det = cxx*cyy-cxy*cxy;
      double nxx = cyy/det;
      double nxy =-cxy/det;
      double nyy = cxx/det;

      double sig=mex*mex*nxx + 2*mex*mey*nxy +mey*mey*nyy;

      double qtx = mpt[0]*cos(mph[0])+mpt[1]*cos(mph[1]);
      double qty = mpt[0]*sin(mph[0])+mpt[1]*sin(mph[1]);
      double qt = sqrt( qtx*qtx + qty*qty );

      double utx = -mex - qtx;
      double uty = -mey - qty;
      double ut = sqrt( utx*utx + uty*uty );

      double ut_par = (utx*qtx + uty*qty)/qt;
      double ut_perp = (uty*qtx - qty*utx)/qt;

      if(FillTree && sig<9.0) ttt->Fill();
      if(FillHist){
         hpc2->Fill(TMath::Prob(sig, 2));
         psig->Fill(nv-1, sig);

         hqt->Fill(qt);
         pqt->Fill(qt, sig);

         hut_par->Fill(ut_par);
         hut_perp->Fill(ut_perp);
         put->Fill(qt, -ut_par/qt);
      }

      m2ll += sig+log(det);
   }

   return m2ll;
}
