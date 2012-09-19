#include "TRandom3.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TTree.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"


double Strue=0.1;
double Ctrue=0.01;

const int NEVENTS=10000;

const float pi = TMath::Pi();

int nj;
double jpt[1000];
double jph[1000];

TTree *t = new TTree("events", "events");

//TH1D *hpc2 = new TH1D("hpc2", "", 100, 0, 1);


bool FILL=false;

double Min2LL(const double *x);


void toy(){
    TRandom3 r(0);
    delete gRandom;
    gRandom = new TRandom3(0);

    t->Branch("nj",  &nj, "nj/I");
    t->Branch("jph", jph, "jph[nj]/D");
    t->Branch("jpt", jpt, "jpt[nj]/D");

    // 'true' jet resolution
    TF1 *freso = new TF1("fresol", "sqrt([0]*[0]*x + [1]*[1]*x*x)", 10, 200);
    freso->SetParameter(0, Strue);
    freso->SetParameter(1, Ctrue);

    // create jets uniform in phi,
    // smeared with freso
    for(Long64_t ev=0; ev<NEVENTS; ++ev){

	int njets   = r.Poisson(20);         if(njets<2) continue;
	int phi0    = r.Uniform(-pi, pi);  
	double pt0  = r.Exp(20);

	nj=njets;
	bool save=true;
	for(int i=0; i<njets; ++i){
	    jph[i] = phi0+2*pi/njets*i;
	    jpt[i] = r.Gaus(pt0, freso->Eval(pt0));
	    if(jpt[i]<0) save = false;
	}

	if(save) t->Fill();
    }

    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    min->SetTolerance(0.001);
    min->SetStrategy(2);
    min->SetPrintLevel(2);
    ROOT::Math::Functor ff(&Min2LL, 2); 
    min->SetFunction(ff);
    min->SetVariable(0, "S", Strue*1.5, 0.01*Strue);
    min->SetVariable(1, "C", Ctrue*0.5, 0.01*Ctrue);
    min->Minimize();				
    FILL=true;
    Min2LL(min->X());
    //hpc2->SetMinimum(0);
    //hpc2->Draw("E");

    cout << "Strue = " << Strue << endl;
    cout << "Ctrue = " << Ctrue << endl;
}


double Min2LL(const double *x){
    double m2ll=0;

    for(Long64_t ev=0; ev<t->GetEntries(); ++ev){

	t->GetEntry(ev);
	
	double mex=0;
	double mey=0;
	double cxx=0;
	double cxy=0;
	double cyy=0;

	for(int i=0; i<nj; ++i){

	    double pt  = jpt[i];
	    double phi = jph[i];
	    double c = cos(phi);
	    double s = sin(phi);

	    mex-=pt*c;
	    mey-=pt*s;

       // parameterize jet resolution
	    double dpt=sqrt(x[0]*x[0]*pt+x[1]*x[1]*pt*pt);
	    double dtt = dpt*dpt;

       // covariance matrix
       // taking zero phi resolution
	    cxx += dtt*c*c+ 0.0*s*s;
	    cxy += c*s*(dtt-0);
	    cyy += 0.0*c*c+ dtt*s*s;
	}

	double det = cxx*cyy-cxy*cxy;
	double nxx = cyy/det;
	double nxy =-cxy/det;
	double nyy = cxx/det;

	double sig=mex*mex*nxx + 2*mex*mey*nxy +mey*mey*nyy;

	//if(FILL) hpc2->Fill(TMath::Prob(sig, 2));

	m2ll += sig+log(det);
    }

    return m2ll;
}
