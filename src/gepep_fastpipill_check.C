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
namespace PSIP_ns{
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
  sprintf(fname,"%s/plot_psip_AR.root",outputdir.c_str());
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

  TH1D *hppi = new TH1D("hppi","pt_{#pi}",200,0,1);
  TH1D *hmpsip  = new TH1D("hm" ,"M(#pi#pi J/#psi)",200,3,4);
  
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
  double start=0.05;
  double stop =0.5;
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
	        lepx4[1],lepy4[1],lepz4[1]);

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
     double pt1 = evt.pip.perp();
     double pt2 = evt.pim.perp();
     //if (p1<0.10 || p1>0.4) continue;
     //if (p2<0.10 || p2>0.4) continue;
     //if (p1+p2<0.4 || p1+p2>0.6) continue;
     costheta1 = evt.GetCostheta1();
     costheta2 = evt.GetCostheta2();
     phi1 = evt.GetPhi1();
     phi2 = evt.GetPhi2();
   //if ( mass>psiplow && mass<psipup ) {
   //  vars->Fill();
   //}
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

     hmpsip->Fill(mass);
     if (mass>psiplow && mass<psipup){
       hppi->Fill(pt1);
       hppi->Fill(pt2);
     }

  }
  //vars->Write();
   
   TFile *ftmp = new TFile("P_cmp.root","update");
   ftmp->WriteTObject(hppi,"hptpi_psip");
   ftmp->WriteTObject(hmpsip,"hmpsip_psip");
   ftmp->Close();
   delete ftmp;
   f->cd();
   return true;


   
   
   for (int parti=0; parti<Npart;parti++)
   for (int partj=0; partj<Npart;partj++){
     sprintf(fname,"part%02d%02d",parti,partj);
	 if(evts[parti][partj].size()<50) continue;
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
 
  int Npart=1;
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
/* 
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
*/


 // Ks decayl/err > 2 
// pcut[0] =0.0 ;  facmap[0] =1.0;       facemap[0] =1.0;     
// pcut[1] =0.05;  facmap[1] =1.00629 ;  facemap[1] =0.0019564;     
// pcut[2] =0.10;  facmap[2] =1.00295 ;  facemap[2] =0.0002664;
// pcut[3] =0.15;  facmap[3] =1.00182 ;  facemap[3] =0.0001468;
// pcut[4] =0.20;  facmap[4] =1.00085 ;  facemap[4] =0.0001360;
// pcut[5] =0.25;  facmap[5] =1.00032 ;  facemap[5] =0.0001242;
// pcut[6] =0.30;  facmap[6] =0.999913;  facemap[6] =0.0001315;
// pcut[7] =0.35;  facmap[7] =0.999792;  facemap[7] =0.0001581;
// pcut[8] =0.40;  facmap[8] =1.00006 ;  facemap[8] =0.0001886;
// pcut[9] =0.45;  facmap[9] =0.998787;  facemap[9] =0.0002521;
// pcut[10]=0.50;  facmap[10]=0.999473;  facemap[10]=0.0002220;
// pcut[11]=0.60;  facmap[11]=0.999828;  facemap[11]=0.0003216;
// pcut[12]=0.70;  facmap[12]=1.00034 ;  facemap[12]=0.0004939;
// pcut[13]=0.80;  facmap[13]=1.00113 ;  facemap[13]=0.0008034;
// pcut[14]=0.90;  facmap[14]=1.00172 ;  facemap[14]=0.0011911;
// pcut[15]=1.00;  facmap[15]=1.00296 ;  facemap[15]=0.0015954;
// pcut[16]=1.20;  facmap[16]=0.98714 ;  facemap[16]=0.0064642;
// pcut[17]=1.40;  facmap[17]=0.987802;  facemap[17]=0.008873 ;
// pcut[18]=1.60;  facmap[18]=1.0;       facemap[18]=1.0;     
// pcut[19]=1.80;  facmap[19]=1.0;       facemap[19]=1.0;     
// pcut[20]=2.00;  //facmap[20]=2.00;


/*
// Ks decayL/err >1 & <2
   pcut[0] =0.0 ;  facmap[0] =1.0;       facemap[0] =1.0;     
   pcut[1] =0.05;  facmap[1] =1.0019  ;  facemap[1] =0.00389584 ;     
   pcut[2] =0.10;  facmap[2] =1.0023  ;  facemap[2] =0.000994288;
   pcut[3] =0.15;  facmap[3] =1.00196 ;  facemap[3] =0.000598345;
   pcut[4] =0.20;  facmap[4] =1.00055 ;  facemap[4] =0.000473927;
   pcut[5] =0.25;  facmap[5] =1.0007  ;  facemap[5] =0.000399789;
   pcut[6] =0.30;  facmap[6] =0.999955;  facemap[6] =0.000536198;
   pcut[7] =0.35;  facmap[7] =1.00092 ;  facemap[7] =0.000681956;
   pcut[8] =0.40;  facmap[8] =1.00312 ;  facemap[8] =0.00109478 ;
   pcut[9] =0.45;  facmap[9] =1.00278 ;  facemap[9] =0.00125143 ;
   pcut[10]=0.50;  facmap[10]=0.99954 ;  facemap[10]=0.00178119 ;
   pcut[11]=0.60;  facmap[11]=1.00235 ;  facemap[11]=0.00226029 ;
   pcut[12]=0.70;  facmap[12]=1.00148 ;  facemap[12]=0.000701396;
   pcut[13]=0.80;  facmap[13]=0.998524;  facemap[13]=0.00405132 ;
   pcut[14]=0.90;  facmap[14]=1.00337 ;  facemap[14]=0.00478848 ;
   pcut[15]=1.00;  facmap[15]=1.04187 ;  facemap[15]=0.020354   ;
   pcut[16]=1.20;  facmap[16]=0.993981;  facemap[16]=0.00342047 ;
   pcut[17]=1.40;  facmap[17]=1.01747 ;  facemap[17]=0.0153603  ;
   pcut[18]=1.60;  facmap[18]=1.0;       facemap[18]=1.0;     
   pcut[19]=1.80;  facmap[19]=1.0;       facemap[19]=1.0;     
   pcut[20]=2.00;  //facmap[20]=2.00;
*/

   pcut[0] =0.0 ;  facmap[0] =1.000610;       facemap[0] =1.0;     
   pcut[1] =0.5;      
  
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
     TF1 ff("ff","1.00065+0.000630*x",0,2);
     factori = ff.Eval(p1);
     factorj = ff.Eval(p2);
     // for average correction factor
   //for (int i=0;i<Npart;i++){
   //  if (p1>=pcut[i]&&p1<pcut[i+1]){
   //     factori = facmap[i];
   //     break;
   //  }
   //}
   //for (int i=0;i<Npart;i++){
   //  if (p2>=pcut[i]&&p2<pcut[i+1]){
   //    factorj = facmap[i];
   //    break;
   //  }
   //}
     if (factori==0 || factorj==0){
       std::cout<<"Waring: factor is 0, id "<<jentry<<", factori "<<factori<<", factorj "<<factorj<<std::endl;
       std::cout<<"p1 "<<p1<<", p2 "<<p2<<", cos1 "<<costheta1<<std::endl;
       continue;
     }
     mass = evts.at(jentry).InvMass(factori,factorj);
     //mass = evts.at(jentry).InvMass(1.001,1.001);
     if (mass>psiplow-0.001 && mass<psipup+0.001)
       dataraw->Fill();
     continue; // other part ignore
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
  PSIP_ns::FitSpectrum(datarawo,tmpchr,true);
  std::cout<<"cccccccccca"<<std::endl;

  // factor at average
  sprintf(tmpchr,"nom_%s",namesfx);
  PSIP_ns::FitSpectrum(dataraw,tmpchr,true);
  std::cout<<"dddddddddda"<<std::endl;
 
//// factor at low edge
//sprintf(tmpchr,"low_%s",namesfx);
//PSIP_ns::FitSpectrum(datarawl,tmpchr,true);

//// factor at up edge
//sprintf(tmpchr,"up_%s",namesfx);
//PSIP_ns::FitSpectrum(datarawu,tmpchr,true);

//~~~~~~~~~~pion part end~~~~~~~~
  return;
}

void PSIP_ns::FitSpectrum(TTree *&dataraw, char* namesfx, bool out)
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
  RooRealVar sigma("sigma","width of gaussian",0.0018,0.001,0.003);
  RooRealVar sigma2("sigma2","width of gaussian",0.005,0.003,0.01);
  RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
  RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
  
  RooRealVar signal("signal"," ",8500,0,1000000);//event number
  RooRealVar sigfra("sigfra"," ",0.5,0.2,1.0);//event number
  //RooRealVar signal2("signal2"," ",100,0,1000000);//event number
  RooRealVar background("background"," ",3000,0,1000000);
  
  RooRealVar a0("a0","coefficient #0",100,-100000,100000);
  RooRealVar a1("a1","coefficient #1",-1,-100000,100000);
  RooPolynomial ground("ground","ground",x,RooArgList(a0,a1));
  
  RooAddPdf sig("sig","signal",RooArgList(gaus,gaus2),sigfra);
  
  //mean.setVal(3.68611);
  //mean.setError(0.0005);
  //sigma.setError(0.0048);
  //sigma2.setError(0.0034);
  //signal.setError(sqrt(signal.getVal()));
  //signal2.setError(sqrt(signal2.getVal()));
  //background.setError(sqrt(background.getVal()));
  
//mean.setVal(3.68611);
//mean.setError(0.000035);
//mean.setError(0.000035);
//sigma.setError(0.000048);
//sigma2.setError(0.000034);
//signal.setError(sqrt(signal.getVal()));
//signal2.setError(sqrt(signal2.getVal()));
//background.setError(sqrt(background.getVal()));
  // Roofit part
  RooAddPdf *sum;
  RooDataSet *dataset;
  RooPlot *xframe;  
  
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
  
  // factor at average
  xframe = x.frame(Title("fit #psi'"));
  sprintf(tmpchr,"data_pi_%s",namesfx);
  dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
  sum = new RooAddPdf("sum","sum",RooArgList(sig,ground),RooArgList(signal,background));
 // sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,ground),RooArgList(signal,signal2,background));
  Npar=8;
  sum->fitTo(*dataset,Range(psiplow,psipup));
  dataset->plotOn(xframe,Binning(nBins));
  sum->plotOn(xframe,Components(sig),LineStyle(2),LineColor(2));
  //sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
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
//pt->AddText(tmpchr);
//sprintf(tmpchr,"#sigma_{2} = %1.6f #pm %1.6f",sigma2.getVal(),sigma2.getError());
  pt->AddText(tmpchr);
  sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
//pt->AddText(tmpchr);
//sprintf(tmpchr,"signal2 = %.2f #pm %.2f",signal2.getVal(),signal2.getError());
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
    sprintf(name,"%s/result.txt",".");//outputdir.c_str());
    ofstream ofresult(name,std::ios::app);
    ofresult<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
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
