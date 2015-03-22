#define gepep_fastpipill_cxx
#include "gepep_fastpipill.h"
#include <TH1D.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include "TPaveText.h"
#include "function.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include <fstream>
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooBreitWigner.h"
#include "RooPolynomial.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooMsgService.h"
#include <vector>
extern std::string outputdir;
using namespace RooFit;
namespace PSIP{
  void FitSpectrum(TTree *&dataraw, char* namesfx, bool out=false);
}

bool gepep_fastpipill::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L gepep_fastpipill.C
//      Root > gepep_fastpipill t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
  //ofstream logfile("log.txt");
  //std::cout.rdbuf(logfile.rdbuf());
  //std::clog.rdbuf(logfile.rdbuf());
  // std::cerr.rdbuf(logfile.rdbuf());
  std::cout<<"Use algorithm: fastpipill_check"<<std::endl;
  if (fChain == 0) return false;
  Long64_t nentries = fChain->GetEntriesFast();

  std::cout<<"Total entry is "<<nentries<<std::endl;
  //int nBins=12;
  double me=0.000511;
  double mmu=0.105658;
  double mpi=0.13957;
  double psiplow=3.665;
  double psipup=3.705;

  // for factor fit
  char fname[1000];
  std::cout<<"0 mass window ("<<psiplow<<", "<<psipup<<")"<<std::endl;
  sprintf(fname,"%s/plot_psip.root",outputdir.c_str());
  std::cout<<"0 mass window ("<<psiplow<<", "<<psipup<<")"<<std::endl;
  TFile *f=new TFile(fname,"RECREATE");

  TCanvas *c1=new TCanvas("c1","",800,600);
  TTree *vars = new TTree("vars","vars");
  double mass;
  double phi,phi1,phi2;
  double costheta,costheta1,costheta2;
  double p1,p2;
  vars->Branch("phi",&phi,"phi/D");
  vars->Branch("phi1",&phi1,"phi1/D");
  vars->Branch("phi2",&phi2,"phi2/D");
  vars->Branch("costheta",&costheta,"costheta/D");
  vars->Branch("costheta1",&costheta1,"costheta1/D");
  vars->Branch("costheta2",&costheta2,"costheta2/D");
  vars->Branch("p1",&p1,"p1/D");
  vars->Branch("p2",&p2,"p2/D");
  vars->Branch("mass",&mass,"mass/D");

  // for initial spectrum
  Long64_t nbytes = 0, nb = 0;

  const int Npart=1;
  // cut phi
  //double phicut[Npart+1];//={0.0,0.5,1.0,1.5,2.0};
  //double start=0.0;
  //double stop =2*TMath::Pi();
  //for(int i=0;i<Npart+1;i++){
  //  phicut[i] = (stop-start)/Npart*i+start;
  //}
  //
  // cut cos theta
  //double costhecut[Npart+1];//={0.0,0.5,1.0,1.5,2.0};
  //double start=-1;
  //double stop =1;
  //for(int i=0;i<Npart+1;i++){
  //  costhecut[i] = (stop-start)/Npart*i+start;
  //}
  // cut p
  double pcut[Npart+1];//={0.0,0.5,1.0,1.5,2.0};
  double start=0.1;
  double stop =0.4;
  for(int i=0;i<Npart+1;i++){
    pcut[i] = (stop-start)/Npart*i+start;
  }

  char name[100];
  // ~~~~~~~~ draw nxn histogram, m distribution in different range
 
  Psip evt;
  std::vector<Psip> evts[Npart][Npart];

  // loop the data
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;

     evt.Setval(pipx4[0],pipy4[0],pipz4[0],
	            pipx4[1],pipy4[1],pipz4[1],
	            lepx4[0],lepy4[0],lepz4[0],
	            lepx4[1],lepy4[1],lepz4[1]
	           );

     // total invariant mass
     if(cos(angle4)>0.90) continue; // cut bhabha 
     if (decay_ee == 0) evt.SetLeptonM(mmu);
     else if (decay_ee==1) evt.SetLeptonM(me);
     else {
       std::cout<<"can not identify lepton. "<<std::endl;
       continue;;
     } 
     // if (Cut(ientry) < 0) continue;
	 //int parti, partj;
	 mass = evt.InvMass();
     //std::cout<<"jentry "<< jentry <<", mass "<<mass<<std::endl;
     //std::cout<<evt.pipx1<<" "<<evt.pipy1<<" "<<evt.pipz1<<" "<<evt.ml<<std::endl;
     p1 = evt.GetP1();
     p2 = evt.GetP2();
	 if (p1<0.10 || p1>0.4) continue;
	 if (p2<0.10 || p2>0.4) continue;
	 //if (p1+p2<0.4 || p1+p2>0.6) continue;
     costheta1 = evt.GetCostheta1();
     costheta2 = evt.GetCostheta2();
     phi1 = evt.GetPhi1();
     phi2 = evt.GetPhi2();
     if ( mass>psiplow && mass<psipup ) {
	   vars->Fill();
	 }
	 int parti=-1,partj=-1;
	 for (int i=0;i<Npart;i++){
       if (p1>=pcut[i]&&p1<pcut[i+1]){
         parti=i;
         break;
       }
     }

     for (int i=0;i<Npart;i++){
       if (p2>=pcut[i]&&p2<pcut[i+1]){
         partj = i;
         break;
       }
     }
     if(parti==-1 || partj==-1) continue;
     if (mass>psiplow-0.002 && mass<psipup+0.002){
	   evts[parti][partj].push_back(evt);
     }

  }
  vars->Write();
   for (int parti=0; parti<Npart;parti++)
   for (int partj=0; partj<Npart;partj++){
     sprintf(fname,"part%02d%02d",parti,partj);
	 if(evts[parti][partj].size()<200) continue;
	 FitSpe(evts[parti][partj],fname);
   }

  //FitSpe(evts,"all");
  return true;
}

void gepep_fastpipill::FitSpe(std::vector<Psip> &evts,const char *namesfx)
{
  double psiplow=3.665;
  double psipup=3.705;
  // for factor fit
 
  TTree *datarawo = new TTree("datarawo","dataraw");
  TTree *dataraw = new TTree("dataraw","dataraw");
  TTree *datarawl = new TTree("datarawl","dataraw");
  TTree *datarawu = new TTree("datarawu","dataraw");
  double mass;
  datarawo->Branch("x",&mass,"x/D");
  dataraw->Branch("x",&mass,"x/D");
  datarawl->Branch("x",&mass,"x/D");
  datarawu->Branch("x",&mass,"x/D");
 
  // try to correct the spectrum
  // iniialize the fit function
  //double factor4,factor4err;// for pi
 
  int Npart=30;
  //int Ncos=10;
  double pcut[Npart+1];//={0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,
		  // 0.60,0.70,0.80,0.90,1.00,1.20,1.40,1.60,1.80,2.00};//={0.0,0.5,1.0,1.5,2.0};
   //pcut[0]=0.1;
   //pcut[1]=0.4;
   //pcut[2]=0.9;
   //double coscut[Ncos+1];
   double  facmap[Npart];
   double facemap[Npart];

//  set normal factor in (0.2, 0.3) to 1.00061, get factor in different range
/*
   pcut[0] =0.0 ;  facmap[0] =1.0;       facemap[0] =1.0;     
   pcut[1] =0.05;  facmap[1] =1.0;       facemap[1] =1.0;     
   pcut[2] =0.10;  facmap[2] =1.0;       facemap[2] =1.0;
   pcut[3] =0.15;  facmap[3] =1.00196 ;  facemap[3] =6.16454e-05;
   pcut[4] =0.20;  facmap[4] =1.00103 ;  facemap[4] =6.03674e-05;
   pcut[5] =0.25;  facmap[5] =1.00037 ;  facemap[5] =7.83072e-05;
   pcut[6] =0.30;  facmap[6] =1.00029 ;  facemap[6] =0.000100928;
   pcut[7] =0.35;  facmap[7] =1.00006 ;  facemap[7] =0.000131676;
   pcut[8] =0.40;  facmap[8] =0.999687;  facemap[8] =0.000174975;
   pcut[9] =0.45;  facmap[9] =0.999197;  facemap[9] =0.000210539;
   pcut[10]=0.50;  facmap[10]=0.999238;  facemap[10]=0.000194209;
   pcut[11]=0.60;  facmap[11]=0.999154;  facemap[11]=0.00028029 ;
   pcut[12]=0.70;  facmap[12]=0.999249;  facemap[12]=0.00044126 ;
   pcut[13]=0.80;  facmap[13]=0.999974;  facemap[13]=0.000672389;
   pcut[14]=0.90;  facmap[14]=1.00274 ;  facemap[14]=0.00108137 ;
   pcut[15]=1.00;  facmap[15]=0.999595;  facemap[15]=0.00119913 ;
   pcut[16]=1.20;  facmap[16]=1.00089 ;  facemap[16]=0.00206829 ;
   pcut[17]=1.40;  facmap[17]=0.988005;  facemap[17]=0.00665777 ;
   pcut[18]=1.60;  facmap[18]=1.0;       facemap[18]=1.0;     
   pcut[19]=1.80;  facmap[19]=1.0;       facemap[19]=1.0;     
   pcut[20]=2.00;  //facmap[20]=2.00;
*/
//  set normal factor in (0.1, 0.4) to 1.000815, get factor in different range
//  only for r value data, limit pi+ split pi-
/*
   pcut[0] =0.0 ;  facmap[0] =1.0;       facemap[0] =1.0;     
   pcut[1] =0.05;  facmap[1] =1.0;       facemap[1] =1.0;     
   pcut[2] =0.10;  facmap[2] =1.00283 ;  facemap[2] =0.000250076;
   pcut[3] =0.15;  facmap[3] =1.00132 ;  facemap[3] =0.000131411;
   pcut[4] =0.20;  facmap[4] =1.00073 ;  facemap[4] =0.000129736;
   pcut[5] =0.25;  facmap[5] =1.00052 ;  facemap[5] =0.000116304;
   pcut[6] =0.30;  facmap[6] =1.00051 ;  facemap[6] =0.00012692 ;
   pcut[7] =0.35;  facmap[7] =1.00015 ;  facemap[7] =0.000155205;
   pcut[8] =0.40;  facmap[8] =0.999884;  facemap[8] =0.000200586;
   pcut[9] =0.45;  facmap[9] =0.999119;  facemap[9] =0.000245204;
   pcut[10]=0.50;  facmap[10]=0.99943 ;  facemap[10]=0.000135434;
   pcut[11]=0.60;  facmap[11]=0.999845;  facemap[11]=0.000352511;
   pcut[12]=0.70;  facmap[12]=0.999938;  facemap[12]=0.000551319;
   pcut[13]=0.80;  facmap[13]=0.999888;  facemap[13]=0.000923838;
   pcut[14]=0.90;  facmap[14]=1.00093 ;  facemap[14]=0.00142603 ;
   pcut[15]=1.00;  facmap[15]=0.999717;  facemap[15]=0.00181855 ;
   pcut[16]=1.20;  facmap[16]=0.993496;  facemap[16]=0.00501287 ;
   pcut[17]=1.40;  facmap[17]=0.981287;  facemap[17]=0.0131427  ;
   pcut[18]=1.60;  facmap[18]=1.0;       facemap[18]=1.0;     
   pcut[19]=1.80;  facmap[19]=1.0;       facemap[19]=1.0;     
   pcut[20]=2.00;  //facmap[20]=2.00;
*/
  //  set normal factor in (0.1, 0.4) to 1.000815, get factor in different range
//  only for r value data, limit pi- split pi+
/*
   pcut[0] =0.0 ;  facmap[0] =1.0;       facemap[0] =1.0;     
   pcut[1] =0.05;  facmap[1] =1.0;       facemap[1] =1.0;     
   pcut[2] =0.10;  facmap[2] =1.00172 ;  facemap[2] =0.000197672;
   pcut[3] =0.15;  facmap[3] =1.00101 ;  facemap[3] =0.000120255;
   pcut[4] =0.20;  facmap[4] =1.00081 ;  facemap[4] =0.000122614;
   pcut[5] =0.25;  facmap[5] =1.00048 ;  facemap[5] =0.000122483;
   pcut[6] =0.30;  facmap[6] =1.00074 ;  facemap[6] =0.000145777;
   pcut[7] =0.35;  facmap[7] =1.0003  ;  facemap[7] =0.000182639;
   pcut[8] =0.40;  facmap[8] =0.999589;  facemap[8] =0.000226131;
   pcut[9] =0.45;  facmap[9] =0.999426;  facemap[9] =0.000266261;
   pcut[10]=0.50;  facmap[10]=0.998955;  facemap[10]=0.000260378;
   pcut[11]=0.60;  facmap[11]=0.999262;  facemap[11]=0.000221771;
   pcut[12]=0.70;  facmap[12]=0.999362;  facemap[12]=0.000590433;
   pcut[13]=0.80;  facmap[13]=0.999897;  facemap[13]=0.000976447;
   pcut[14]=0.90;  facmap[14]=1.00053 ;  facemap[14]=0.00145681 ;
   pcut[15]=1.00;  facmap[15]=0.998737;  facemap[15]=0.00191489 ;
   pcut[16]=1.20;  facmap[16]=1.01641 ;  facemap[16]=0.0126093  ;
   pcut[17]=1.40;  facmap[17]=1.00226 ;  facemap[17]=0.00331346 ;
   pcut[18]=1.60;  facmap[18]=1.0;       facemap[18]=1.0;     
   pcut[19]=1.80;  facmap[19]=1.0;       facemap[19]=1.0;     
   pcut[20]=2.00;  //facmap[20]=2.00;
*/
   //  set normal factor in (0.1, 0.4) to 1.000815, get factor in different range
//  only for r value data, combine both part
/*
   pcut[0] =0.0 ;  facmap[0] =1.0;       facemap[0] =1.0;     
   pcut[1] =0.05;  facmap[1] =1.0;       facemap[1] =1.0;     
   pcut[2] =0.10;  facmap[2] =1.00215 ;  facemap[2] =0.000155076;
   pcut[3] =0.15;  facmap[3] =1.00115 ;  facemap[3] =8.87154e-05;
   pcut[4] =0.20;  facmap[4] =1.00077 ;  facemap[4] =8.91127e-05;
   pcut[5] =0.25;  facmap[5] =1.0005  ;  facemap[5] =8.43392e-05;
   pcut[6] =0.30;  facmap[6] =1.00061 ;  facemap[6] =9.57233e-05;
   pcut[7] =0.35;  facmap[7] =1.00021 ;  facemap[7] =0.000118269;
   pcut[8] =0.40;  facmap[8] =0.999754;  facemap[8] =0.000150058;
   pcut[9] =0.45;  facmap[9] =0.99926 ;  facemap[9] =0.000180371;
   pcut[10]=0.50;  facmap[10]=0.999329;  facemap[10]=0.000120152;
   pcut[11]=0.60;  facmap[11]=0.999427;  facemap[11]=0.000187713;
   pcut[12]=0.70;  facmap[12]=0.99967 ;  facemap[12]=0.00040296 ;
   pcut[13]=0.80;  facmap[13]=0.999892;  facemap[13]=0.00067108 ;
   pcut[14]=0.90;  facmap[14]=1.00073 ;  facemap[14]=0.00101906 ;
   pcut[15]=1.00;  facmap[15]=0.999252;  facemap[15]=0.00131865 ;
   pcut[16]=1.20;  facmap[16]=0.996623;  facemap[16]=0.00465825 ;
   pcut[17]=1.40;  facmap[17]=1.00101 ;  facemap[17]=0.00321292 ;
   pcut[18]=1.60;  facmap[18]=1.0;       facemap[18]=1.0;     
   pcut[19]=1.80;  facmap[19]=1.0;       facemap[19]=1.0;     
   pcut[20]=2.00;  //facmap[20]=2.00;
*/
   //  set normal factor in (0.2, 0.3) to 1.00061, get factor in different range
//  only for r value data, combine both part
/*
   pcut[0] =0.0 ;  facmap[0] =1.0;       facemap[0] =1.0;     
   pcut[1] =0.05;  facmap[1] =1.0;       facemap[1] =1.0;     
   pcut[2] =0.10;  facmap[2] =1.0;       facemap[2] =1.0;
   pcut[3] =0.15;  facmap[3] =1.00193 ;  facemap[3] =0.00014646 ;
   pcut[4] =0.20;  facmap[4] =1.00092 ;  facemap[4] =0.00014185 ;
   pcut[5] =0.25;  facmap[5] =1.00027 ;  facemap[5] =0.00017408 ;
   pcut[6] =0.30;  facmap[6] =1.00028 ;  facemap[6] =0.00025041 ;
   pcut[7] =0.35;  facmap[7] =0.999665;  facemap[7] =0.00028309 ;
   pcut[8] =0.40;  facmap[8] =0.999421;  facemap[8] =0.00040939 ;
   pcut[9] =0.45;  facmap[9] =0.998626;  facemap[9] =0.00048584 ;
   pcut[10]=0.50;  facmap[10]=0.999499;  facemap[10]=0.00043914 ;
   pcut[11]=0.60;  facmap[11]=0.998844;  facemap[11]=0.00070548 ;
   pcut[12]=0.70;  facmap[12]=0.998429;  facemap[12]=0.00139006 ;
   pcut[13]=0.80;  facmap[13]=1.00149 ;  facemap[13]=0.00189529 ;
   pcut[14]=0.90;  facmap[14]=1.00292 ;  facemap[14]=0.00278649 ;
   pcut[15]=1.00;  facmap[15]=1.00281 ;  facemap[15]=0.00383812 ;
   pcut[16]=1.20;  facmap[16]=1.02283 ;  facemap[16]=0.0085055  ;
   pcut[17]=1.40;  facmap[17]=1.0;       facemap[17]=1.0;
   pcut[18]=1.60;  facmap[18]=1.0;       facemap[18]=1.0;     
   pcut[19]=1.80;  facmap[19]=1.0;       facemap[19]=1.0;     
   pcut[20]=2.00;  //facmap[20]=2.00;
*/
   //  set factors in (0.15, 0.6) to up values, get factor in different range
//  only for r value data, combine both part
/*
   pcut[0] =0.0 ;  facmap[0] =1.0;       facemap[0] =1.0;     
   pcut[1] =0.05;  facmap[1] =1.01015 ;  facemap[1] =0.00242001;     
   pcut[2] =0.10;  facmap[2] =1.00338 ;  facemap[2] =0.00025419;
   pcut[3] =0.15;  facmap[3] =1.00186 ;  facemap[3] =0.00013994;
   pcut[4] =0.20;  facmap[4] =1.00092 ;  facemap[4] =0.00013154;
   pcut[5] =0.25;  facmap[5] =1.00019 ;  facemap[5] =0.00012394;
   pcut[6] =0.30;  facmap[6] =0.999856;  facemap[6] =0.00015996;
   pcut[7] =0.35;  facmap[7] =0.999711;  facemap[7] =0.00016778;
   pcut[8] =0.40;  facmap[8] =0.999856;  facemap[8] =0.00023577;
   pcut[9] =0.45;  facmap[9] =0.998991;  facemap[9] =0.00027663;
   pcut[10]=0.50;  facmap[10]=0.999686;  facemap[10]=0.00026146;
   pcut[11]=0.60;  facmap[11]=1.00061 ;  facemap[11]=0.00036703;
   pcut[12]=0.70;  facmap[12]=1.00096 ;  facemap[12]=0.00057256;
   pcut[13]=0.80;  facmap[13]=1.00081 ;  facemap[13]=0.00093882;
   pcut[14]=0.90;  facmap[14]=1.00299 ;  facemap[14]=0.00142871;
   pcut[15]=1.00;  facmap[15]=1.00248 ;  facemap[15]=0.00182948;
   pcut[16]=1.20;  facmap[16]=0.994156;  facemap[16]=0.00528607;
   pcut[17]=1.40;  facmap[17]=0.989396;  facemap[17]=0.00765296;
   pcut[18]=1.60;  facmap[18]=1.0;       facemap[18]=1.0;     
   pcut[19]=1.80;  facmap[19]=1.0;       facemap[19]=1.0;     
   pcut[20]=2.00;  //facmap[20]=2.00;
*/
   //  set factors in (0.15, 0.6) to up values, get factor in different range
//  only for r value data, combine both part

   pcut[0] = 0.0  ;  facmap[0] = 1.0     ;  facemap[0] = 1.0     ;
   pcut[1] = 0.025;  facmap[1] = 1.0     ;  facemap[1] = 1.0     ;
   pcut[2] = 0.05 ;  facmap[2] = 1.0     ;  facemap[2] = 1.0     ;
   pcut[3] = 0.075;  facmap[3] = 1.00919 ;  facemap[3] = 0.00191258 ;
   pcut[4] = 0.10 ;  facmap[4] = 1.00559 ;  facemap[4] = 0.000552638;
   pcut[5] = 0.125;  facmap[5] = 1.00256 ;  facemap[5] = 0.000273443;
   pcut[6] = 0.15 ;  facmap[6] = 1.00197 ;  facemap[6] = 0.000202182;
   pcut[7] = 0.175;  facmap[7] = 1.00177 ;  facemap[7] = 0.000193008;
   pcut[8] = 0.20 ;  facmap[8] = 1.00127 ;  facemap[8] = 0.000189631;
   pcut[9] = 0.225;  facmap[9] = 1.00058 ;  facemap[9] = 0.00017915 ;
   pcut[10]= 0.25 ;  facmap[10]= 1.00048 ;  facemap[10]= 0.000172806;
   pcut[11]= 0.275;  facmap[11]= 0.999882;  facemap[11]= 0.000178129;
   pcut[12]= 0.30 ;  facmap[12]= 0.999897;  facemap[12]= 0.000200975;
   pcut[13]= 0.325;  facmap[13]= 0.999849;  facemap[13]= 0.000233133;
   pcut[14]= 0.35 ;  facmap[14]= 0.999335;  facemap[14]= 0.000277518;
   pcut[15]= 0.375;  facmap[15]= 0.99998 ;  facemap[15]= 0.000284219;
   pcut[16]= 0.40 ;  facmap[16]= 0.999656;  facemap[16]= 0.000311098;
   pcut[17]= 0.425;  facmap[17]= 1.00012 ;  facemap[17]= 0.00035162 ;
   pcut[18]= 0.45 ;  facmap[18]= 0.998626;  facemap[18]= 0.000378786;
   pcut[19]= 0.475;  facmap[19]= 0.999348;  facemap[19]= 0.000407338;
   pcut[20]= 0.50 ;  facmap[20]= 0.999698;  facemap[20]= 0.000259805;
   pcut[21]= 0.60 ;  facmap[21]= 1.00061 ;  facemap[21]= 0.000366871;
   pcut[22]= 0.70 ;  facmap[22]= 1.00096 ;  facemap[22]= 0.000572888;
   pcut[23]= 0.80 ;  facmap[23]= 1.00081 ;  facemap[23]= 0.000938436;
   pcut[24]= 0.90 ;  facmap[24]= 1.00299 ;  facemap[24]= 0.00143016 ;
   pcut[25]= 1.00 ;  facmap[25]= 1.00247 ;  facemap[25]= 0.00182197 ;
   pcut[26]= 1.20 ;  facmap[26]= 0.994152;  facemap[26]= 0.00529237 ;
   pcut[27]= 1.40 ;  facmap[27]= 0.989408;  facemap[27]= 0.00764222 ;
   pcut[28]= 1.60 ;  facmap[28]= 1.0     ;  facemap[28]= 1.0     ;
   pcut[29]= 1.80 ;  facmap[29]= 1.0     ;  facemap[29]= 1.0     ;
   pcut[30]= 2.00 ;

 //pcut[0] = 0.0  ; //coscut[0] = -1.0;
 //pcut[1] = 0.025; //coscut[1] = -0.8;   
 //pcut[2] = 0.05 ; //coscut[2] = -0.6;
 //pcut[3] = 0.075; //coscut[3] = -0.4;
 //pcut[4] = 0.10 ; //coscut[4] = -0.2;   
 //pcut[5] = 0.125; //coscut[5] =  0.0;   
 //pcut[6] = 0.15 ; //coscut[6] =  0.2;   
 //pcut[7] = 0.175; //coscut[7] =  0.4;   
 //pcut[8] = 0.20 ; //coscut[8] =  0.6;   
 //pcut[9] = 0.225; //coscut[9] =  0.8;   
 //pcut[10]= 0.25 ; //coscut[10]=  1.0;   
 //pcut[11]= 0.275;   
 //pcut[12]= 0.30 ;   
 //pcut[13]= 0.325;   
 //pcut[14]= 0.35 ;   
 //pcut[15]= 0.375;   
 //pcut[16]= 0.40 ;   
 //pcut[17]= 0.425;   
 //pcut[18]= 0.45 ;   
 //pcut[19]= 0.475;   
 //pcut[20]= 0.50 ;   
 //pcut[21]= 0.60 ;   
 //pcut[22]= 0.70 ;   
 //pcut[23]= 0.80 ;   
 //pcut[24]= 0.90 ;   
 //pcut[25]= 1.00 ;   
 //pcut[26]= 1.20 ;   
 //pcut[27]= 1.40 ;   
 //pcut[28]= 1.60 ;   
 //pcut[29]= 1.80 ;   
 //pcut[30]= 2.00 ;   
/*
  facmap[0][0] =  1. ;  facemap[0][0] =  1.;
  facmap[0][1] =  1. ;  facemap[0][1] =  1.;
  facmap[0][2] =  1. ;  facemap[0][2] =  1.;
  facmap[0][3] =  1. ;  facemap[0][3] =  1.;
  facmap[0][4] =  1. ;  facemap[0][4] =  1.;
  facmap[0][5] =  1. ;  facemap[0][5] =  1.;
  facmap[0][6] =  1. ;  facemap[0][6] =  1.;
  facmap[0][7] =  1. ;  facemap[0][7] =  1.;
  facmap[0][8] =  1. ;  facemap[0][8] =  1.;
  facmap[0][9] =  1. ;  facemap[0][9] =  1.;
  facmap[1][0] =  1. ;  facemap[1][0] =  1.;
  facmap[1][1] =  1. ;  facemap[1][1] =  1.;
  facmap[1][2] =  1. ;  facemap[1][2] =  1.;
  facmap[1][3] =  1. ;  facemap[1][3] =  1.;
  facmap[1][4] =  1. ;  facemap[1][4] =  1.;
  facmap[1][5] =  1. ;  facemap[1][5] =  1.;
  facmap[1][6] =  1. ;  facemap[1][6] =  1.;
  facmap[1][7] =  1. ;  facemap[1][7] =  1.;
  facmap[1][8] =  1. ;  facemap[1][8] =  1.;
  facmap[1][9] =  1. ;  facemap[1][9] =  1.;
  facmap[2][0] =  1. ;  facemap[2][0] =  1.;
  facmap[2][1] =  1. ;  facemap[2][1] =  1.;
  facmap[2][2] =  1. ;  facemap[2][2] =  1.;
  facmap[2][3] =  1. ;  facemap[2][3] =  1.;
  facmap[2][4] =  1. ;  facemap[2][4] =  1.;
  facmap[2][5] =  1. ;  facemap[2][5] =  1.;
  facmap[2][6] =  1. ;  facemap[2][6] =  1.;
  facmap[2][7] =  1. ;  facemap[2][7] =  1.;
  facmap[2][8] =  1. ;  facemap[2][8] =  1.;
  facmap[2][9] =  1. ;  facemap[2][9] =  1.;
  facmap[3][0] =  1.02549  ;  facemap[3][0] = 0.057221   ;
  facmap[3][1] =  1.00705  ;  facemap[3][1] = 0.00655242 ;
  facmap[3][2] =  0.999087 ;  facemap[3][2] = 0.00674187 ;
  facmap[3][3] =  1.00403  ;  facemap[3][3] = 0.00280597 ;
  facmap[3][4] =  1.00985  ;  facemap[3][4] = 0.00348275 ;
  facmap[3][5] =  1.00937  ;  facemap[3][5] = 0.00430175 ;
  facmap[3][6] =  1.00833  ;  facemap[3][6] = 0.00332264 ;
  facmap[3][7] =  1.0148   ;  facemap[3][7] = 0.00705234 ;
  facmap[3][8] =  1.01301  ;  facemap[3][8] = 0.0170698  ;
  facmap[3][9] =  1.0152   ;  facemap[3][9] = 0.0139321  ;
  facmap[4][0] =  0.99162  ;  facemap[4][0] = 0.0189203  ;
  facmap[4][1] =  1.00686  ;  facemap[4][1] = 0.00243901 ;
  facmap[4][2] =  1.00176  ;  facemap[4][2] = 0.00135396 ;
  facmap[4][3] =  1.00301  ;  facemap[4][3] = 0.00114467 ;
  facmap[4][4] =  1.00485  ;  facemap[4][4] = 0.00033575 ;
  facmap[4][5] =  1.00536  ;  facemap[4][5] = 0.00143191 ;
  facmap[4][6] =  1.00459  ;  facemap[4][6] = 0.00127563 ;
  facmap[4][7] =  1.00835  ;  facemap[4][7] = 0.00186224 ;
  facmap[4][8] =  1.01304  ;  facemap[4][8] = 0.00318143 ;
  facmap[4][9] =  1.0259   ;  facemap[4][9] = 0.00711918 ;
  facmap[5][0] =  0.995606 ;  facemap[5][0] = 0.00516567 ;
  facmap[5][1] =  0.999095 ;  facemap[5][1] = 0.00124138 ;
  facmap[5][2] =  0.999311 ;  facemap[5][2] = 0.00079172 ;
  facmap[5][3] =  1.00009  ;  facemap[5][3] = 0.00062837 ;
  facmap[5][4] =  1.0022   ;  facemap[5][4] = 0.00059747 ;
  facmap[5][5] =  1.00334  ;  facemap[5][5] = 0.00068348 ;
  facmap[5][6] =  1.00425  ;  facemap[5][6] = 0.00068889 ;
  facmap[5][7] =  1.00505  ;  facemap[5][7] = 0.00084526 ;
  facmap[5][8] =  1.00804  ;  facemap[5][8] = 0.00167035 ;
  facmap[5][9] =  1.00503  ;  facemap[5][9] = 0.00537048 ;
  facmap[6][0] =  1.00025  ;  facemap[6][0] = 0.00243444  ;
  facmap[6][1] =  0.997159 ;  facemap[6][1] = 0.000714729 ;
  facmap[6][2] =  0.999243 ;  facemap[6][2] = 0.000526227 ;
  facmap[6][3] =  1.00145  ;  facemap[6][3] = 0.00048751  ;
  facmap[6][4] =  1.00162  ;  facemap[6][4] = 0.000514444 ;
  facmap[6][5] =  1.00174  ;  facemap[6][5] = 0.00050462  ;
  facmap[6][6] =  1.00364  ;  facemap[6][6] = 0.000546889 ;
  facmap[6][7] =  1.00366  ;  facemap[6][7] = 0.000575031 ;
  facmap[6][8] =  1.00529  ;  facemap[6][8] = 0.000810584 ;
  facmap[6][9] =  1.00897  ;  facemap[6][9] = 0.00376518  ;
  facmap[7][0] =  0.998596 ;  facemap[7][0] = 0.00181152  ;
  facmap[7][1] =  0.998134 ;  facemap[7][1] = 0.000651994 ;
  facmap[7][2] =  0.998765 ;  facemap[7][2] = 0.000506631 ;
  facmap[7][3] =  0.999931 ;  facemap[7][3] = 0.000482376 ;
  facmap[7][4] =  1.00084  ;  facemap[7][4] = 0.000490128 ;
  facmap[7][5] =  1.00127  ;  facemap[7][5] = 0.000515299 ;
  facmap[7][6] =  1.00374  ;  facemap[7][6] = 0.000556172 ;
  facmap[7][7] =  1.00418  ;  facemap[7][7] = 0.000560098 ;
  facmap[7][8] =  1.00615  ;  facemap[7][8] = 0.000720018 ;
  facmap[7][9] =  1.007    ;  facemap[7][9] = 0.00225047  ;
  facmap[8][0] =  0.995945 ;  facemap[8][0] = 0.00191146 ;
  facmap[8][1] =  0.996704 ;  facemap[8][1] = 0.00057081 ;
  facmap[8][2] =  0.999492 ;  facemap[8][2] = 0.00052991 ;
  facmap[8][3] =  0.999143 ;  facemap[8][3] = 0.00053043 ;
  facmap[8][4] =  1.00014  ;  facemap[8][4] = 0.00048051 ;
  facmap[8][5] =  1.00161  ;  facemap[8][5] = 0.00051432 ;
  facmap[8][6] =  1.00336  ;  facemap[8][6] = 0.00054595 ;
  facmap[8][7] =  1.00331  ;  facemap[8][7] = 0.00054867 ;
  facmap[8][8] =  1.00581  ;  facemap[8][8] = 0.00066739 ;
  facmap[8][9] =  1.00888  ;  facemap[8][9] = 0.00216271 ;
  facmap[9][0] =  0.998014 ;  facemap[9][0] = 0.00147886 ;
  facmap[9][1] =  0.996598 ;  facemap[9][1] = 0.00051613 ;
  facmap[9][2] =  0.997804 ;  facemap[9][2] = 0.00053125 ;
  facmap[9][3] =  0.999737 ;  facemap[9][3] = 0.00054737 ;
  facmap[9][4] =  1.00042  ;  facemap[9][4] = 0.00055447 ;
  facmap[9][5] =  0.999856 ;  facemap[9][5] = 0.00049475 ;
  facmap[9][6] =  1.00264  ;  facemap[9][6] = 0.00055467 ;
  facmap[9][7] =  1.00265  ;  facemap[9][7] = 0.00046703 ;
  facmap[9][8] =  1.00379  ;  facemap[9][8] = 0.00056445 ;
  facmap[9][9] =  1.0065   ;  facemap[9][9] = 0.0017228  ;
  facmap[10][0] = 0.9957    ;  facemap[10][0] = 0.00125198 ;
  facmap[10][1] = 0.99751   ;  facemap[10][1] = 0.00050237 ;
  facmap[10][2] = 0.997254  ;  facemap[10][2] = 0.00052201 ;
  facmap[10][3] = 0.998913  ;  facemap[10][3] = 0.00046551 ;
  facmap[10][4] = 1.00028   ;  facemap[10][4] = 0.00041825 ;
  facmap[10][5] = 1.0004    ;  facemap[10][5] = 0.00043782 ;
  facmap[10][6] = 1.00212   ;  facemap[10][6] = 0.00050719 ;
  facmap[10][7] = 1.00272   ;  facemap[10][7] = 0.00051503 ;
  facmap[10][8] = 1.00422   ;  facemap[10][8] = 0.00061772 ;
  facmap[10][9] = 1.00621   ;  facemap[10][9] = 0.0016175  ;
  facmap[11][0] = 0.99672   ;  facemap[11][0] = 0.00135414 ;
  facmap[11][1] = 0.996279  ;  facemap[11][1] = 0.00059237 ;
  facmap[11][2] = 0.996757  ;  facemap[11][2] = 0.00053281 ;
  facmap[11][3] = 0.998351  ;  facemap[11][3] = 0.00046704 ;
  facmap[11][4] = 0.999418  ;  facemap[11][4] = 0.00048226 ;
  facmap[11][5] = 1.00049   ;  facemap[11][5] = 0.00048814 ;
  facmap[11][6] = 1.00164   ;  facemap[11][6] = 0.00048534 ;
  facmap[11][7] = 1.00307   ;  facemap[11][7] = 0.00058806 ;
  facmap[11][8] = 1.00358   ;  facemap[11][8] = 0.00062555 ;
  facmap[11][9] = 1.00156   ;  facemap[11][9] = 0.00134879 ;
  facmap[12][0] = 0.995771  ;  facemap[12][0] = 0.00122583 ;
  facmap[12][1] = 0.996813  ;  facemap[12][1] = 0.00063765 ;
  facmap[12][2] = 0.997734  ;  facemap[12][2] = 0.00052545 ;
  facmap[12][3] = 0.998646  ;  facemap[12][3] = 0.00056599 ;
  facmap[12][4] = 0.999109  ;  facemap[12][4] = 0.00059419 ;
  facmap[12][5] = 1.00022   ;  facemap[12][5] = 0.00059641 ;
  facmap[12][6] = 1.00105   ;  facemap[12][6] = 0.00055426 ;
  facmap[12][7] = 1.00281   ;  facemap[12][7] = 0.00058796 ;
  facmap[12][8] = 1.00396   ;  facemap[12][8] = 0.00065522 ;
  facmap[12][9] = 1.00155   ;  facemap[12][9] = 0.00142705 ;
  facmap[13][0] = 0.996244  ;  facemap[13][0] = 0.00153883 ;
  facmap[13][1] = 0.994997  ;  facemap[13][1] = 0.00069699 ;
  facmap[13][2] = 0.996522  ;  facemap[13][2] = 0.00065030 ;
  facmap[13][3] = 0.997438  ;  facemap[13][3] = 0.00067146 ;
  facmap[13][4] = 1.00015   ;  facemap[13][4] = 0.00070481 ;
  facmap[13][5] = 1.00098   ;  facemap[13][5] = 0.00077353 ;
  facmap[13][6] = 1.00262   ;  facemap[13][6] = 0.00074108 ;
  facmap[13][7] = 1.00236   ;  facemap[13][7] = 0.00067589 ;
  facmap[13][8] = 1.00457   ;  facemap[13][8] = 0.00074389 ;
  facmap[13][9] = 1.0043    ;  facemap[13][9] = 0.00136408 ;
  facmap[14][0] = 0.998371  ;  facemap[14][0] = 0.00126655 ;
  facmap[14][1] = 0.996886  ;  facemap[14][1] = 0.00069161 ;
  facmap[14][2] = 0.997737  ;  facemap[14][2] = 0.00068289 ;
  facmap[14][3] = 0.997049  ;  facemap[14][3] = 0.00084988 ;
  facmap[14][4] = 0.998171  ;  facemap[14][4] = 0.0008696  ;
  facmap[14][5] = 0.999393  ;  facemap[14][5] = 0.00089048 ;
  facmap[14][6] = 1.00174   ;  facemap[14][6] = 0.00086111 ;
  facmap[14][7] = 1.00197   ;  facemap[14][7] = 0.00075576 ;
  facmap[14][8] = 1.00333   ;  facemap[14][8] = 0.00069103 ;
  facmap[14][9] = 1.00309   ;  facemap[14][9] = 0.00147684 ;
  facmap[15][0] = 0.998735  ;  facemap[15][0] = 0.00126323 ;
  facmap[15][1] = 0.996166  ;  facemap[15][1] = 0.00069323 ;
  facmap[15][2] = 0.997322  ;  facemap[15][2] = 0.00081431 ;
  facmap[15][3] = 0.998059  ;  facemap[15][3] = 0.00089935 ;
  facmap[15][4] = 0.998851  ;  facemap[15][4] = 0.00099198 ;
  facmap[15][5] = 1.00182   ;  facemap[15][5] = 0.00113327 ;
  facmap[15][6] = 1.00156   ;  facemap[15][6] = 0.00098787 ;
  facmap[15][7] = 1.00381   ;  facemap[15][7] = 0.00089346 ;
  facmap[15][8] = 1.00303   ;  facemap[15][8] = 0.00066777 ;
  facmap[15][9] = 1.00299   ;  facemap[15][9] = 0.00189277 ;
  facmap[16][0] = 0.995642  ;  facemap[16][0] = 0.00144649 ;
  facmap[16][1] = 0.994434  ;  facemap[16][1] = 0.00080775 ;
  facmap[16][2] = 0.996004  ;  facemap[16][2] = 0.00096512 ;
  facmap[16][3] = 0.997221  ;  facemap[16][3] = 0.00115292 ;
  facmap[16][4] = 0.999209  ;  facemap[16][4] = 0.00114151 ;
  facmap[16][5] = 1.00143   ;  facemap[16][5] = 0.00113603 ;
  facmap[16][6] = 1.00458   ;  facemap[16][6] = 0.00123215 ;
  facmap[16][7] = 1.00339   ;  facemap[16][7] = 0.00089603 ;
  facmap[16][8] = 1.00332   ;  facemap[16][8] = 0.00075895 ;
  facmap[16][9] = 1.00277   ;  facemap[16][9] = 0.00129257 ;
  facmap[17][0] = 0.994591  ;  facemap[17][0] = 0.00135669 ;
  facmap[17][1] = 0.995838  ;  facemap[17][1] = 0.00077980 ;
  facmap[17][2] = 0.995069  ;  facemap[17][2] = 0.00114158 ;
  facmap[17][3] = 0.998819  ;  facemap[17][3] = 0.00136935 ;
  facmap[17][4] = 1.00038   ;  facemap[17][4] = 0.00128658 ;
  facmap[17][5] = 1.00086   ;  facemap[17][5] = 0.00096132 ;
  facmap[17][6] = 1.00351   ;  facemap[17][6] = 0.00133465 ;
  facmap[17][7] = 1.00509   ;  facemap[17][7] = 0.00115317 ;
  facmap[17][8] = 1.0026    ;  facemap[17][8] = 0.00090339 ;
  facmap[17][9] = 1.0049    ;  facemap[17][9] = 0.00127241 ;
  facmap[18][0] = 0.996623  ;  facemap[18][0] = 0.00115259 ;
  facmap[18][1] = 0.993112  ;  facemap[18][1] = 0.001056   ;
  facmap[18][2] = 0.994203  ;  facemap[18][2] = 0.00171305 ;
  facmap[18][3] = 0.994541  ;  facemap[18][3] = 0.00147581 ;
  facmap[18][4] = 0.997952  ;  facemap[18][4] = 0.00147878 ;
  facmap[18][5] = 1.0015    ;  facemap[18][5] = 0.00155263 ;
  facmap[18][6] = 1.00419   ;  facemap[18][6] = 0.00153495 ;
  facmap[18][7] = 1.00009   ;  facemap[18][7] = 0.00130044 ;
  facmap[18][8] = 1.00292   ;  facemap[18][8] = 0.00096275 ;
  facmap[18][9] = 1.00347   ;  facemap[18][9] = 0.00130109 ;
  facmap[19][0] = 0.995813  ;  facemap[19][0] = 0.0011995  ;
  facmap[19][1] = 0.993016  ;  facemap[19][1] = 0.00115693 ;
  facmap[19][2] = 0.992949  ;  facemap[19][2] = 0.00155091 ;
  facmap[19][3] = 0.996068  ;  facemap[19][3] = 0.00146682 ;
  facmap[19][4] = 0.997903  ;  facemap[19][4] = 0.00158468 ;
  facmap[19][5] = 1.00227   ;  facemap[19][5] = 0.0017465  ;
  facmap[19][6] = 1.00623   ;  facemap[19][6] = 0.00180336 ;
  facmap[19][7] = 1.0034    ;  facemap[19][7] = 0.001334   ;
  facmap[19][8] = 1.00472   ;  facemap[19][8] = 0.00107403 ;
  facmap[19][9] = 1.00264   ;  facemap[19][9] = 0.00113923 ;
  facmap[20][0] = 0.995317  ;  facemap[20][0] = 0.00078434 ;
  facmap[20][1] = 0.993119  ;  facemap[20][1] = 0.00073807 ;
  facmap[20][2] = 0.995303  ;  facemap[20][2] = 0.00085557 ;
  facmap[20][3] = 0.997651  ;  facemap[20][3] = 0.00095582 ;
  facmap[20][4] = 0.999011  ;  facemap[20][4] = 0.00103992 ;
  facmap[20][5] = 1.00088   ;  facemap[20][5] = 0.00100889 ;
  facmap[20][6] = 1.00369   ;  facemap[20][6] = 0.00108691 ;
  facmap[20][7] = 1.0053    ;  facemap[20][7] = 0.00093303 ;
  facmap[20][8] = 1.00585   ;  facemap[20][8] = 0.00074471 ;
  facmap[20][9] = 1.00331   ;  facemap[20][9] = 0.00076694 ;
  facmap[21][0] = 0.997555  ;  facemap[21][0] = 0.00096530 ;
  facmap[21][1] = 0.992248  ;  facemap[21][1] = 0.00116718 ;
  facmap[21][2] = 0.996053  ;  facemap[21][2] = 0.00103887 ;
  facmap[21][3] = 0.999313  ;  facemap[21][3] = 0.00133111 ;
  facmap[21][4] = 1.00155   ;  facemap[21][4] = 0.00143243 ;
  facmap[21][5] = 1.00383   ;  facemap[21][5] = 0.00152881 ;
  facmap[21][6] = 1.00285   ;  facemap[21][6] = 0.00137629 ;
  facmap[21][7] = 1.00635   ;  facemap[21][7] = 0.00134917 ;
  facmap[21][8] = 1.00387   ;  facemap[21][8] = 0.00098909 ;
  facmap[21][9] = 1.00423   ;  facemap[21][9] = 0.00104381 ;
  facmap[22][0] = 0.995169  ;  facemap[22][0] = 0.00165177 ;
  facmap[22][1] = 0.993842  ;  facemap[22][1] = 0.00151265 ;
  facmap[22][2] = 0.997564  ;  facemap[22][2] = 0.0016798  ;
  facmap[22][3] = 1.00104   ;  facemap[22][3] = 0.00213766 ;
  facmap[22][4] = 1.00565   ;  facemap[22][4] = 0.002785   ;
  facmap[22][5] = 1.00195   ;  facemap[22][5] = 0.00200829 ;
  facmap[22][6] = 1.00325   ;  facemap[22][6] = 0.0026152  ;
  facmap[22][7] = 1.00584   ;  facemap[22][7] = 0.00223157 ;
  facmap[22][8] = 1.00246   ;  facemap[22][8] = 0.00136671 ;
  facmap[22][9] = 1.00438   ;  facemap[22][9] = 0.00158491 ;
  facmap[23][0] = 0.99278   ;  facemap[23][0] = 0.00097255 ;
  facmap[23][1] = 0.995206  ;  facemap[23][1] = 0.00248032 ;
  facmap[23][2] = 0.992371  ;  facemap[23][2] = 0.00396394 ;
  facmap[23][3] = 0.995349  ;  facemap[23][3] = 0.00427303 ;
  facmap[23][4] = 1.0003    ;  facemap[23][4] = 0.00294538 ;
  facmap[23][5] = 1.00142   ;  facemap[23][5] = 0.00359154 ;
  facmap[23][6] = 1.00584   ;  facemap[23][6] = 0.00139165 ;
  facmap[23][7] = 1.00421   ;  facemap[23][7] = 0.00383721 ;
  facmap[23][8] = 1.00581   ;  facemap[23][8] = 0.00250709 ;
  facmap[23][9] = 1.00171   ;  facemap[23][9] = 0.00235628 ;
  facmap[24][0] = 0.999002  ;  facemap[24][0] = 0.00326114 ;
  facmap[24][1] = 0.996641  ;  facemap[24][1] = 0.00303282 ;
  facmap[24][2] = 0.995125  ;  facemap[24][2] = 0.00437531 ;
  facmap[24][3] = 1.00918   ;  facemap[24][3] = 0.0079375  ;
  facmap[24][4] = 1.00619   ;  facemap[24][4] = 0.0048033  ;
  facmap[24][5] = 1.0094    ;  facemap[24][5] = 0.00895654 ;
  facmap[24][6] = 1.00337   ;  facemap[24][6] = 0.00745333 ;
  facmap[24][7] = 1.0071    ;  facemap[24][7] = 0.00693747 ;
  facmap[24][8] = 1.00672   ;  facemap[24][8] = 0.00353916 ;
  facmap[24][9] = 1.0069    ;  facemap[24][9] = 0.00417591 ;
*/

  char tmpchr[100];

//~~~~~~~~~~pion part start~~~~~~~~

  for (Long64_t jentry=0; jentry<evts.size();jentry++) {
     
     double p1,p2;
     double costheta1,costheta2;
     // total invariant mass
     p1 = evts.at(jentry).GetP1();
     p2 = evts.at(jentry).GetP2();
	 costheta1 = evts.at(jentry).GetCostheta1();
	 costheta2 = evts.at(jentry).GetCostheta2();
     if (p1<0.10 || p1>1.0) continue;
     if (p2<0.10 || p2>1.0) continue;
	 // without correction
     mass = evts.at(jentry).InvMass();
     if (mass>psiplow-0.001 && mass<psipup+0.001)
       datarawo->Fill();
 
	 double factori=0,factorj=0;
	 // for average correction factor
     for (int i=0;i<Npart;i++){
        if (p1>=pcut[i]&&p1<pcut[i+1]){
           factori = facmap[i];
		   break;
        }
      }
      for (int i=0;i<Npart;i++){
        if (p2>=pcut[i]&&p2<pcut[i+1]){
          factorj = facmap[i];
          break;
        }
      }
	 if (factori==0 || factorj==0){
	   std::cout<<"Waring: factor is 0, id "<<jentry<<", factori "<<factori<<", factorj "<<factorj<<std::endl;
	   std::cout<<"p1 "<<p1<<", p2 "<<p2<<", cos1 "<<costheta1<<std::endl;
	   continue;
	 }
     mass = evts.at(jentry).InvMass(factori,factorj);
     //mass = evts.at(jentry).InvMass(1.001,1.001);
     if (mass>psiplow-0.001 && mass<psipup+0.001)
       dataraw->Fill();
	 // average factor end
	 //
	 // factor at low edge
     for (int i=0;i<Npart;i++){
        if (p1>=pcut[i]&&p1<pcut[i+1]){
              factori = facmap[i]-facemap[i];
			  break;
        }
      }
      for (int i=0;i<Npart;i++){
        if (p2>=pcut[i]&&p2<pcut[i+1]){
              factorj = facmap[i]-facemap[i];
			  break;
        }
      }
     if (factori==0 || factorj==0){
	   std::cout<<"Waring: factor is 0, id "<<jentry<<std::endl;
	   continue;
	 }
     mass = evts.at(jentry).InvMass(factori,factorj);
     //mass = evts.at(jentry).InvMass(1.0005,1.0005);
     if (mass>psiplow-0.001 && mass<psipup+0.001)
       datarawl->Fill();
	 // factor at low edge end
	 //
	 // factor at up edge
     for (int i=0;i<Npart;i++){
        if (p1>=pcut[i]&&p1<pcut[i+1]){
              factori = facmap[i]+facemap[i];
			  break;
        }
      }
      for (int i=0;i<Npart;i++){
        if (p2>=pcut[i]&&p2<pcut[i+1]){
              factorj = facmap[i]+facemap[i];
			  break;
        }
      }
     if (factori==0 || factorj==0){
	   std::cout<<"Waring: factor is 0, id "<<jentry<<std::endl;
	   continue;
	 }
     mass = evts.at(jentry).InvMass(factori,factorj);
     //mass = evts.at(jentry).InvMass(1.002,1.002);
     if (mass>psiplow-0.001 && mass<psipup+0.001)
       datarawu->Fill();
  }
  //dataraw->Write();
  // no correction
  sprintf(tmpchr,"raw_%s",namesfx);
  PSIP::FitSpectrum(datarawo,tmpchr,true);
  std::cout<<"cccccccccca"<<std::endl;

  // factor at average
  sprintf(tmpchr,"nom_%s",namesfx);
  PSIP::FitSpectrum(dataraw,tmpchr,true);
  std::cout<<"dddddddddda"<<std::endl;
 
  // factor at low edge
  sprintf(tmpchr,"low_%s",namesfx);
  PSIP::FitSpectrum(datarawl,tmpchr,true);

  // factor at up edge
  sprintf(tmpchr,"up_%s",namesfx);
  PSIP::FitSpectrum(datarawu,tmpchr,true);

//~~~~~~~~~~pion part end~~~~~~~~
  return;
}

void PSIP::FitSpectrum(TTree *&dataraw, char* namesfx, bool out)
{
  int nBins=100;
  double mparticle = 3.686109;
  double psiplow=3.665;
  double psipup=3.705;
  int Npar;
  char tmpchr[100];
  TCanvas *c1=new TCanvas("c1_1","",800,600);
  // try to use roofit
  RooRealVar x("x","energy",mparticle,psiplow,psipup,"GeV");
  RooRealVar mean("mean","mean of gaussian",mparticle,psiplow,psipup);
  RooRealVar sigma("sigma","width of gaussian",0.0017,0.001,0.0025);
  RooRealVar sigma2("sigma2","width of gaussian",0.0025,0.002,0.004);
  RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
  RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
  
  RooRealVar signal("signal"," ",3800,0,1000000);//event number
  RooRealVar signal2("signal2"," ",3000,0,1000000);//event number
  RooRealVar background("background"," ",2000,0,1000000);
  
  RooRealVar a0("a0","coefficient #0",100,-100000,100000);
  RooRealVar a1("a1","coefficient #1",-1,-100000,100000);
  RooPolynomial ground("ground","ground",x,RooArgList(a0,a1));
 
  mean.setVal(3.68611);
  mean.setError(0.000035);
  sigma.setError(0.000048);
  sigma2.setError(0.000034);
  signal.setError(sqrt(signal.getVal()));
  signal2.setError(sqrt(signal2.getVal()));
  background.setError(sqrt(background.getVal()));
  // Roofit part
  RooAddPdf *sum;
  RooDataSet *dataset;
  RooPlot *xframe;  
  
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
  
  // factor at average
  xframe = x.frame(Title("fit #psi'"));
  sprintf(tmpchr,"data_pi_%s",namesfx);
  dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
  sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,ground),RooArgList(signal,signal2,background));
  Npar=8;
  sum->fitTo(*dataset,Range(psiplow,psipup));
  dataset->plotOn(xframe,Binning(nBins));
  sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
  sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
  sum->plotOn(xframe,Components(ground),LineStyle(3),LineColor(3));
  sum->plotOn(xframe);
  xframe->Draw();
  TPaveText *pt = new TPaveText(0.60,0.5,0.90,0.90,"BRNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(4000);
  pt->SetTextAlign(12);
  pt->SetTextFont(42);
  pt->SetTextSize(0.035);
  sprintf(tmpchr,"#mu_{1} = %1.6f #pm %1.6f",mean.getVal(),mean.getError());
  pt->AddText(tmpchr);
  sprintf(tmpchr,"#sigma_{1} = %1.6f #pm %1.6f",sigma.getVal(),sigma.getError());
  pt->AddText(tmpchr);
  sprintf(tmpchr,"#sigma_{2} = %1.6f #pm %1.6f",sigma2.getVal(),sigma2.getError());
  pt->AddText(tmpchr);
  sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
  pt->AddText(tmpchr);
  sprintf(tmpchr,"signal2 = %.2f #pm %.2f",signal2.getVal(),signal2.getError());
  pt->AddText(tmpchr);
  sprintf(tmpchr,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
  pt->AddText(tmpchr);
  sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
  pt->AddText(tmpchr);
  pt->Draw();
  //c1->Update();
  //xframe->SetName("mass_pi_best");
  sprintf(tmpchr,"mass_spectrum_%s",namesfx);
  c1->SetName(tmpchr);
  c1->Write();
  if (out){
    char name[500];
    sprintf(name,"%s/result.txt",outputdir.c_str());
    ofstream ofresult(name,std::ios::app);
    ofresult<<namesfx<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
    ofresult.close();
  }
  std::cout<<"aaaaaaaaaaa"<<std::endl;
  delete sum;
  delete dataset;
  delete xframe;
  delete pt;
  delete c1;
  std::cout<<"bbbbbbbbbba"<<std::endl;

  return;
}

gepep_fastpipill::gepep_fastpipill(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data_Rvalue_pipill_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("data_Rvalue_pipill_1.root");
      }
      f->GetObject("gepep_fastpipill",tree);

   }
   Init(tree);
}

gepep_fastpipill::~gepep_fastpipill()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gepep_fastpipill::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gepep_fastpipill::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void gepep_fastpipill::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("rec", &rec, &b_rec);
   fChain->SetBranchAddress("evttag", &evttag, &b_evttag);
   fChain->SetBranchAddress("indexmc", &indexmc, &b_indexmc);
   fChain->SetBranchAddress("pdgid", pdgid, &b_pdgid);
   fChain->SetBranchAddress("motheridx", motheridx, &b_motheridx);
   fChain->SetBranchAddress("mcgpx", mcgpx, &b_mcgpx);
   fChain->SetBranchAddress("mcgpy", mcgpy, &b_mcgpy);
   fChain->SetBranchAddress("mcgpz", mcgpz, &b_mcgpz);
   fChain->SetBranchAddress("mcge", mcge, &b_mcge);
   fChain->SetBranchAddress("mcgid", mcgid, &b_mcgid);
   fChain->SetBranchAddress("nGammatch", &nGammatch, &b_nGammatch);
   fChain->SetBranchAddress("ncharg", &ncharg, &b_ncharg);
   fChain->SetBranchAddress("ntot", &ntot, &b_ntot);
   fChain->SetBranchAddress("nneu", &nneu, &b_nneu);
   fChain->SetBranchAddress("ngch", &ngch, &b_ngch);
   fChain->SetBranchAddress("ngam", &ngam, &b_ngam);
   fChain->SetBranchAddress("npi0", &npi0, &b_npi0);
   fChain->SetBranchAddress("netap2", &netap2, &b_netap2);
   fChain->SetBranchAddress("delang", delang, &b_delang);
   fChain->SetBranchAddress("delphi", delphi, &b_delphi);
   fChain->SetBranchAddress("delthe", delthe, &b_delthe);
   fChain->SetBranchAddress("npart", npart, &b_npart);
   fChain->SetBranchAddress("nemchits", nemchits, &b_nemchits);
   fChain->SetBranchAddress("module", module, &b_module);
   fChain->SetBranchAddress("x", x, &b_x);
   fChain->SetBranchAddress("y", y, &b_y);
   fChain->SetBranchAddress("z", z, &b_z);
   fChain->SetBranchAddress("px", px, &b_px);
   fChain->SetBranchAddress("py", py, &b_py);
   fChain->SetBranchAddress("pz", pz, &b_pz);
   fChain->SetBranchAddress("theta", theta, &b_theta);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("dx", dx, &b_dx);
   fChain->SetBranchAddress("dy", dy, &b_dy);
   fChain->SetBranchAddress("dz", dz, &b_dz);
   fChain->SetBranchAddress("dtheta", dtheta, &b_dtheta);
   fChain->SetBranchAddress("dphi", dphi, &b_dphi);
   fChain->SetBranchAddress("energy", energy, &b_energy);
   fChain->SetBranchAddress("dE", dE, &b_dE);
   fChain->SetBranchAddress("eSeed", eSeed, &b_eSeed);
   fChain->SetBranchAddress("nSeed", nSeed, &b_nSeed);
   fChain->SetBranchAddress("e3x3", e3x3, &b_e3x3);
   fChain->SetBranchAddress("e5x5", e5x5, &b_e5x5);
   fChain->SetBranchAddress("secondMoment", secondMoment, &b_secondMoment);
   fChain->SetBranchAddress("latMoment", latMoment, &b_latMoment);
   fChain->SetBranchAddress("a20Moment", a20Moment, &b_a20Moment);
   fChain->SetBranchAddress("a42Moment", a42Moment, &b_a42Moment);
   fChain->SetBranchAddress("getTime", getTime, &b_getTime);
   fChain->SetBranchAddress("getEAll", getEAll, &b_getEAll);
   fChain->SetBranchAddress("mpi0", mpi0, &b_mpi0);
   fChain->SetBranchAddress("chisq1cpi0", chisq1cpi0, &b_chisq1cpi0);
   fChain->SetBranchAddress("ig1pi0", ig1pi0, &b_ig1pi0);
   fChain->SetBranchAddress("ig2pi0", ig2pi0, &b_ig2pi0);
   fChain->SetBranchAddress("chisq6c", &chisq6c, &b_chisq6c);
   fChain->SetBranchAddress("chisq4c", &chisq4c, &b_chisq4c);
   fChain->SetBranchAddress("chisq3c", &chisq3c, &b_chisq3c);
   fChain->SetBranchAddress("gpx4", gpx4, &b_gpx4);
   fChain->SetBranchAddress("gpy4", gpy4, &b_gpy4);
   fChain->SetBranchAddress("gpz4", gpz4, &b_gpz4);
   fChain->SetBranchAddress("ge4", ge4, &b_ge4);
   fChain->SetBranchAddress("lepx4", lepx4, &b_lepx4);
   fChain->SetBranchAddress("lepy4", lepy4, &b_lepy4);
   fChain->SetBranchAddress("lepz4", lepz4, &b_lepz4);
   fChain->SetBranchAddress("lee4", lee4, &b_lee4);
   fChain->SetBranchAddress("pipx4", pipx4, &b_pipx4);
   fChain->SetBranchAddress("pipy4", pipy4, &b_pipy4);
   fChain->SetBranchAddress("pipz4", pipz4, &b_pipz4);
   fChain->SetBranchAddress("pie4", pie4, &b_pie4);
   fChain->SetBranchAddress("ggm4", &ggm4, &b_ggm4);
   fChain->SetBranchAddress("kkm4", &kkm4, &b_kkm4);
   fChain->SetBranchAddress("kpi01m4", &kpi01m4, &b_kpi01m4);
   fChain->SetBranchAddress("kpi02m4", &kpi02m4, &b_kpi02m4);
   fChain->SetBranchAddress("kkpi0m4", &kkpi0m4, &b_kkpi0m4);
   fChain->SetBranchAddress("kkpipim4", &kkpipim4, &b_kkpipim4);
   fChain->SetBranchAddress("kkpipipi0m4", &kkpipipi0m4, &b_kkpipipi0m4);
   fChain->SetBranchAddress("Recoilmass", &Recoilmass, &b_Recoilmass);
   fChain->SetBranchAddress("llm4", &llm4, &b_llm4);
   fChain->SetBranchAddress("pipim4", &pipim4, &b_pipim4);
   fChain->SetBranchAddress("llpipi1m4", &llpipi1m4, &b_llpipi1m4);
   fChain->SetBranchAddress("llpipi2m4", &llpipi2m4, &b_llpipi2m4);
   fChain->SetBranchAddress("llpipi3m4", &llpipi3m4, &b_llpipi3m4);
   fChain->SetBranchAddress("llpipi4m4", &llpipi4m4, &b_llpipi4m4);
   fChain->SetBranchAddress("gpipim4", &gpipim4, &b_gpipim4);
   fChain->SetBranchAddress("decay_ee", &decay_ee, &b_decay_ee);
   fChain->SetBranchAddress("hepp", &hepp, &b_hepp);
   fChain->SetBranchAddress("hepm", &hepm, &b_hepm);
   fChain->SetBranchAddress("emcp", &emcp, &b_emcp);
   fChain->SetBranchAddress("emcm", &emcm, &b_emcm);
   fChain->SetBranchAddress("mucp", &mucp, &b_mucp);
   fChain->SetBranchAddress("mucm", &mucm, &b_mucm);
   fChain->SetBranchAddress("zcpm4", &zcpm4, &b_zcpm4);
   fChain->SetBranchAddress("zcmm4", &zcmm4, &b_zcmm4);
   fChain->SetBranchAddress("ecms", &ecms, &b_ecms);
   fChain->SetBranchAddress("angle1", &angle1, &b_angle1);
   fChain->SetBranchAddress("angle2", &angle2, &b_angle2);
   fChain->SetBranchAddress("angle3", &angle3, &b_angle3);
   fChain->SetBranchAddress("angle4", &angle4, &b_angle4);
   fChain->SetBranchAddress("psipm4", &psipm4, &b_psipm4);
   Notify();
}

Bool_t gepep_fastpipill::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gepep_fastpipill::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gepep_fastpipill::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
