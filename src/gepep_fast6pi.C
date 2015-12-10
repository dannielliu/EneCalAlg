#define gepep_fast6pi_cxx
#include "gepep_fast6pi.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include "TPaveText.h"
#include "function.h"
#include <fstream>
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooCBShape.h"
#include "RooArgusBG.h"
#include "RooChebychev.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooMsgService.h"
#include <vector>
#include <iomanip>
#include <sstream>
extern std::string outputdir;
using namespace RooFit;
using namespace std;

namespace EEto6PI{
  void FitSpe(std::vector<EEto6pi> evts, double beame,  char* namesfx=0);
  void FitSpectrum(TTree *&dataraw, double beame, char* namesfx=0);
  void FitSpectrum2(const TH1* h, double beame, char* namesfx=0);
  void FitSpectrum3(const TH1* h, double beame, char* namesfx=0); // crystal ball
  void FitSpectrum4(const TH1* h, double beame, char* namesfx=0); // gaus
  //double GetEnergy(int runNo);
}

void gepep_fast6pi::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L gepep_fast6pi.C
//      Root > gepep_fast6pi t
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
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();

   fChain->GetEntry(1);
   double beamene = GetEnergy(run);
   //beamene = 4.415;
   ifstream enefile("enetxt");
   if(enefile.is_open()) {
      enefile>>beamene;
      beamene = beamene/1000;
      enefile.close();
   }
   std::cout<<"current beam energy is "<<beamene<<", run id "<<run<<std::endl;
   if (beamene < 0.1){
     std::cout<<"can not get a suitable beam energy!"<<std::endl;
	 return;
   }
   std::cout<<"Toral entry is "<<nentries<<std::endl;
   int nBins=100;
   double peakvalue=beamene;// mbeam
   double beamlow=beamene-0.25;
   double beamup=beamene+0.25;
   double mpi=0.13957;
    
  char fname[1000];
  sprintf(fname,"%s/plot_6pi_xyz.root",outputdir.c_str());
  TFile *f=new TFile(fname,"update");

  char name[100];
  TCanvas *c1=new TCanvas("c1","",800,600);
  TTree *vars = new TTree("vars","vars");
  double mass;
  //double phi,phi1,phi2;
  //double costheta,costheta1,costheta2;
  double p[6];
   vars->Branch("indexmc", &indexmc, "indexmc/I");
   vars->Branch("pdgid", pdgid, "pdgid[indexmc]/I");
   vars->Branch("motheridx", motheridx, "motheridx[indexmc]/I");
  vars->Branch("p1",&p[0],"p1/D");
  vars->Branch("p2",&p[1],"p2/D");
  vars->Branch("p3",&p[2],"p3/D");
  vars->Branch("p4",&p[3],"p4/D");
  vars->Branch("p5",&p[4],"p5/D");
  vars->Branch("p6",&p[5],"p6/D");
  vars->Branch("mass",&mass,"mass/D");
  
   EEto6pi evt;
   std::vector<EEto6pi> evts;

   // select useful events
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;	  
	 	  
      if (beamene == 4.26){
   	if (run < 30368) continue;
      }

      evt.Setval(pippx,pippy,pippz,pimpx,pimpy,pimpz);
      mass = evt.InvMass();
      bool good=true;;
      for (int i=0; i<6;i++) {
         p[i] = evt.GetP(i);
         //if (p[i]<0.1 || p[i]>1.2) good = false;
      }
      //check error from p distribution
      if (!good) continue;
      
      if (mass>beamlow && mass<beamup){
          evts.push_back(evt);
      }
      if (mass>peakvalue-0.1 && mass<peakvalue+0.1)    vars->Fill();
   }//select end
   vars->Write();
   sprintf(name,"%f",peakvalue);
   EEto6PI::FitSpe(evts,peakvalue,name);
   return;
}

void EEto6PI::FitSpe(std::vector<EEto6pi> evts, double beame, char *namesfx)
{
  double beamlow=beame-0.20;
  double beamup=beame+0.20;
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

  TH1D *hm = new TH1D("hm","hm",100,beamlow,beamup);
  TH1D *hmo = new TH1D("hmo","hm",100,beamlow,beamup);
  TH1D *hml = new TH1D("hml","hm",100,beamlow,beamup);
  TH1D *hmu = new TH1D("hmu","hm",100,beamlow,beamup);

  // try to correct the spectrum
  // iniialize the fit function
  //double factor4,factor4err;// for pi
 
  int Npart=2;
  double pcut[Npart+1];//={0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,
		  // 0.60,0.70,0.80,0.90,1.00,1.20,1.40,1.60,1.80,2.00};//={0.0,0.5,1.0,1.5,2.0};
   //pcut[0]=0.1;
   //pcut[1]=0.4;
   //pcut[2]=0.9;
   double facmap[Npart];
   double facemap[Npart];

   //  set normal factor in (0.2, 0.3) to 1.00061, get factor in different range
//  only for r value data, combine both part
// factor from Ks->pipi, pi corrected with vertex fit , low p range use factor from pipill

 
// factor from Ks dl/dle >2
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
   pcut[0] = 0.0;  facmap[0] = 1.000902; facemap[0] = 0.000077;
   pcut[1] = 0.4;  facmap[1] = 1.000769; facemap[1] = 0.000117;
   pcut[2] = 2.0;


/*
// factor from Ks dl/dle <2 & >1
   pcut[0] =0.0 ;  facmap[0] =1.0;       facemap[0] =1.0;     
   pcut[1] =0.05;  facmap[1] =1.0019  ;  facemap[1] =0.0038958;
   pcut[2] =0.10;  facmap[2] =1.0023  ;  facemap[2] =0.0009942;
   pcut[3] =0.15;  facmap[3] =1.00196 ;  facemap[3] =0.0005983;
   pcut[4] =0.20;  facmap[4] =1.00055 ;  facemap[4] =0.0004739;
   pcut[5] =0.25;  facmap[5] =1.0007  ;  facemap[5] =0.0003997;
   pcut[6] =0.30;  facmap[6] =0.999955;  facemap[6] =0.0005361;
   pcut[7] =0.35;  facmap[7] =1.00092 ;  facemap[7] =0.0006819;
   pcut[8] =0.40;  facmap[8] =1.00312 ;  facemap[8] =0.0010947;
   pcut[9] =0.45;  facmap[9] =1.00278 ;  facemap[9] =0.0012514;
   pcut[10]=0.50;  facmap[10]=0.99954 ;  facemap[10]=0.0017811;
   pcut[11]=0.60;  facmap[11]=1.00235 ;  facemap[11]=0.0022602;
   pcut[12]=0.70;  facmap[12]=1.00148 ;  facemap[12]=0.0007013;
   pcut[13]=0.80;  facmap[13]=0.998524;  facemap[13]=0.0040513;
   pcut[14]=0.90;  facmap[14]=1.00337 ;  facemap[14]=0.0047884;
   pcut[15]=1.00;  facmap[15]=1.04187 ;  facemap[15]=0.020354 ;
   pcut[16]=1.20;  facmap[16]=0.993981;  facemap[16]=0.0034204;
   pcut[17]=1.40;  facmap[17]=1.01747 ;  facemap[17]=0.0153603;
   pcut[18]=1.60;  facmap[18]=1.0;       facemap[18]=1.0;     
   pcut[19]=1.80;  facmap[19]=1.0;       facemap[19]=1.0;     
   pcut[20]=2.00;  //facmap[20]=2.00;
*/


  // psi 3.097
     
  
  char tmpchr[100];
  double fpilow=0, fpih=0, fpier = 0;
  ifstream inpar("CORF");
  if (inpar.is_open()){
    string line;
    istringstream isstream;
    while (!inpar.eof()){
      getline(inpar,line);
      if (line[0]=='#' || line[0]==' ') continue;
      isstream.str(line.c_str());
      sprintf(tmpchr,"%s",line.c_str());
      if (strncmp(tmpchr,"f pi",4) == 0){
        isstream.read(tmpchr,4);
	while (!isstream.eof()){
	  isstream >> line;
	  if (line == ":") {
	    isstream >> fpih;
	    isstream >> line;
	    isstream >> fpier;
	    break;
	  }
	}
      }
    }
    if (fpih ==0 ){
      std::cout<<" Can not find suitable correction factor in file: CORF"<<std::endl;
    }
    else{
       facmap[1] = fpih; facemap[1] = fpier;
    }
  }
  std::cout << "fpih = " << fpih << ", fer = " << fpier << std::endl;
  inpar.close();

  //~~~~~~~~~~part start~~~~~~~~
  //TF1 ff("ff","1.00065+0.000630*x",0,2);
  TF1 ff("ff","1.00113-0.00091738*x",0,2);

  for (Long64_t jentry=0; jentry<evts.size();jentry++) {
     
     double p[6];
     double f[6];
     double fl[6];
     double fu[6];
     // total invariant mass
	 // without correction
     mass = evts.at(jentry).InvMass();
     if (mass>beamlow-0.001 && mass<beamup+0.001) {datarawo->Fill();hmo->Fill(mass);}
     for (int i=0;i<6;i++) p[i] = evts.at(jentry).GetP(i);

     short flag = 0x0;
     for (int pid=0;pid<6;pid++){
       f[pid] = ff.Eval(p[pid]);
       fl[pid] = ff.Eval(p[pid]);
       fu[pid] = ff.Eval(p[pid]);
     //for (int i=0;i<Npart;i++){
     //  if (p[pid]>=pcut[i]&&p[pid]<pcut[i+1]){
     //    f[pid]  = facmap[i];
     //    fl[pid] = facmap[i]-facemap[i];
     //    fu[pid] = facmap[i]+facemap[i];
     //    flag = (flag<<1) +1; //if factor find, the flag should be 00111111
     //    break;
     //  }
     //}
     }
   //if (flag!=0x3f){
   //   std::cout<<"Waring: not good event, factor map is "<<flag<<std::endl;
   //   continue;
   //}
	 // for average correction factor
     mass = evts.at(jentry).InvMass(f);
     if (mass>beamlow-0.001 && mass<beamup+0.001) {dataraw->Fill();hm->Fill(mass);}
	 // factor at low edge
     mass = evts.at(jentry).InvMass(fl);
     if (mass>beamlow-0.001 && mass<beamup+0.001) {datarawl->Fill();hml->Fill(mass);}
	 // factor at up edge
     mass = evts.at(jentry).InvMass(fu);
     if (mass>beamlow-0.001 && mass<beamup+0.001) {datarawu->Fill();hmu->Fill(mass);}
  }
  //dataraw->Write();
  // no correction
  sprintf(tmpchr,"raw_%s",namesfx);
  EEto6PI::FitSpectrum(datarawo,beame,tmpchr);
  //EEto6PI::FitSpectrum2(hmo,beame,tmpchr);
  //EEto6PI::FitSpectrum3(hmo,beame,tmpchr);
  //EEto6PI::FitSpectrum4(hmo,beame,tmpchr);

  // for mc data
  // factor at average
  sprintf(tmpchr,"nom_%s",namesfx);
  EEto6PI::FitSpectrum(dataraw,beame,tmpchr);
  //EEto6PI::FitSpectrum2(hm,beame,tmpchr);
  //EEto6PI::FitSpectrum3(hm,beame,tmpchr);
  //EEto6PI::FitSpectrum4(hm,beame,tmpchr);
 
//// factor at low edge
//sprintf(tmpchr,"low_%s",namesfx);
////EEto6PI::FitSpectrum(datarawl,beame,tmpchr);
//EEto6PI::FitSpectrum2(hml,beame,tmpchr);

//// factor at up edge
//sprintf(tmpchr,"upv_%s",namesfx);
////EEto6PI::FitSpectrum(datarawu,beame,tmpchr);
//EEto6PI::FitSpectrum2(hmu,beame,tmpchr);

  //~~~~~~~~~~ part end~~~~~~~~
  return;
}

void EEto6PI::FitSpectrum(TTree *&dataraw, double beame, char* namesfx)
{
   int nBins=100;
   int Npar;
   double peakvalue = beame;
   //double beamlow=beame-0.15;
   double beamlow=beame-0.10;
   double beamup=beame+0.05;
   // try to use roofit
   RooRealVar x("x","energy",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue+0.003,beamlow,beamup);
   //RooRealVar mean2("mean2","mean of gaussian",peakvalue /*-0.002*/ ,beamlow,beamup);
   RooRealVar sigma("sigma","width of gaussian",0.010,0.005,0.013);
   //RooRealVar sigma2("sigma2","width of gaussian",0.03,0.02,0.05);
   //RooRealVar sigma3("sigma3","width of gaussian",0.005,0.002,0.008);
   //if (strncmp(namesfx,"raw",3) == 0) mean.setVal(peakvalue - 0.005);

   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   //RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean2,sigma2);
   //RooGaussian gaus3("gaus3","gauss(x,m,s)",x,mean,sigma3);
     RooRealVar co1("co1","coefficient #1",-0.8,-100.,100.);
     RooRealVar co2("co2","coefficient #1",-0.2,-100.,100.);
     RooRealVar co3("co3","coefficient #1",0.1,-100.,100.);
     //RooRealVar co4("co4","coefficient #4",0);
     //RooChebychev ground("ground","background",x,RooArgList(co1,co2,co3));
   
   RooRealVar am("am","argus m", beame+0.05, 0, 10);
   RooRealVar ac("ac","argus c", 0, 0, 10);
   RooRealVar ap("ap","argus p", 0.2, 0, 10);
   RooArgusBG ground("ground","background",x,am,ac,ap);
   
   RooRealVar signal("signal"," ",300,0,1000000);//event number
   //RooRealVar signal1("signal1"," ",20,0,1000000);//event number
   RooRealVar background("background"," ",200,0,100000);
 //RooRealVar a0("a0","coefficient #0",100,-100000,100000);
 //RooRealVar a1("a1","coefficient #1",-10,-100000,100000);
 //RooRealVar a2("a2","coefficient #2",0,-100000,100000);
   //RooRealVar a3("a3","coefficient #2",0,-100000,100000);
 //RooPolynomial ground("ground","ground",x,RooArgList(a0,a1,a2));
     
   RooRealVar alpha1("alpha1","#alpha",1.,-5,5);
   RooRealVar nnn1("n1","n",100,1,200);
   RooCBShape cbshape1("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);

   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   //RooDataHist *data_6pi;
 
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   TCanvas *c1 = new TCanvas("","",800,600);

   char tmpchr[100];
   sprintf(tmpchr,"data_6pi_%s",namesfx);
   //data_6pi = new RooDataHist(tmpchr,"data_6pi",x,h);
   xframe = x.frame(Title("fit 6 pi"));
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   //sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2),RooArgList(signal,background));
   //sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,gaus3),RooArgList(signal,background,signal1));
   sum = new RooAddPdf("sum","sum",RooArgList(cbshape1,ground),RooArgList(signal,background));
   Npar = 8;
   //sigma.setVal(0.035);
   //signal.setVal(1200);
   //background.setVal(200);
   //co1.setVal(0);
   sum->fitTo(*dataset,Range(beamlow,beamup));
   dataset->plotOn(xframe);
   sum->plotOn(xframe,Components(cbshape1),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(ground),LineStyle(2),LineColor(3));
   //sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
   //sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(3));
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
  sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
  pt->AddText(tmpchr);
  sprintf(tmpchr,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
  pt->AddText(tmpchr);
  sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
  pt->AddText(tmpchr);
  pt->Draw();
  sprintf(tmpchr,"mass_spectrum_%s",namesfx);
  c1->SetName(tmpchr);
  c1->Write();

   ofstream outf("f6pi",std::ios::app);
   outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   
   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   delete dataset;
   delete sum;
   return;
}

void EEto6PI::FitSpectrum2(const TH1* hm, double beame, char* namesfx)
{
   int nBins=100;
   int Npar;
   double peakvalue = beame;
   //double beamlow=beame-0.15;
   //double beamlow=beame-0.02;
   //double beamup=beame+0.02;
   double beamlow=beame-0.15;
   double beamup=beame+0.10;
   // try to use roofit
   RooRealVar x("x","energy",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue,beamlow,beamup);
   //RooRealVar mean2("mean2","mean of gaussian",peakvalue /*-0.002*/ ,beamlow,beamup);
   RooRealVar sigma("sigma","width of gaussian",0.012,0.008,0.015);
   //RooRealVar sigma2("sigma2","width of gaussian",0.03,0.02,0.05);
   //RooRealVar sigma3("sigma3","width of gaussian",0.005,0.002,0.008);
   if (strncmp(namesfx,"raw",3) == 0) mean.setVal(peakvalue - 0.003);

   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   //RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean2,sigma2);
   //RooGaussian gaus3("gaus3","gauss(x,m,s)",x,mean,sigma3);
     RooRealVar co1("co1","coefficient #1",-0.8,-100.,100.);
     RooRealVar co2("co2","coefficient #1",-0.2,-100.,100.);
     RooRealVar co3("co3","coefficient #1",0.1,-100.,100.);
     //RooRealVar co4("co4","coefficient #4",0);
     RooChebychev ground("ground","background",x,RooArgList(co1,co2,co3));
   RooRealVar signal("signal"," ",300,0,1000000);//event number
   //RooRealVar signal1("signal1"," ",20,0,1000000);//event number
   RooRealVar background("background"," ",200,0,100000);
 //RooRealVar a0("a0","coefficient #0",100,-100000,100000);
 //RooRealVar a1("a1","coefficient #1",-10,-100000,100000);
 //RooRealVar a2("a2","coefficient #2",0,-100000,100000);
   //RooRealVar a3("a3","coefficient #2",0,-100000,100000);
 //RooPolynomial ground("ground","ground",x,RooArgList(a0,a1,a2));
     
   RooRealVar alpha1("alpha1","#alpha",1.,-5,5);
   RooRealVar nnn1("n1","n",100,1,200);
   RooCBShape cbshape1("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);

   RooAddPdf *sum;
   //RooDataSet *dataset;
   RooDataHist *datahist;
   RooPlot *xframe;
   //RooDataHist *data_6pi;
 
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   TCanvas *c1 = new TCanvas("","",800,600);

   char tmpchr[100];
   sprintf(tmpchr,"data_6pi_%s",namesfx);
   //data_6pi = new RooDataHist(tmpchr,"data_6pi",x,h);
   xframe = x.frame(Title("fit 6 pi"));
   //dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   datahist = new RooDataHist(tmpchr,"data",RooArgSet(x),hm);
   //sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2),RooArgList(signal,background));
   //sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,gaus3),RooArgList(signal,background,signal1));
   sum = new RooAddPdf("sum","sum",RooArgList(cbshape1,ground),RooArgList(signal,background));
   Npar = 8;
   //sigma.setVal(0.035);
   //signal.setVal(1200);
   //background.setVal(200);
   //co1.setVal(0);
   //sum->fitTo(*datahist,Range(beamlow,beamup));
   gaus.fitTo(*datahist,Range(beame-0.012,beame+0.012));
   datahist->plotOn(xframe);
   gaus.plotOn(xframe);
   //sum->plotOn(xframe,Components(cbshape1),LineStyle(2),LineColor(2));
   //sum->plotOn(xframe,Components(ground),LineStyle(2),LineColor(3));
   //sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
   //sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(3));
   //sum->plotOn(xframe);
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
//sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
//pt->AddText(tmpchr);
//sprintf(tmpchr,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
//pt->AddText(tmpchr);
  sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
  pt->AddText(tmpchr);
  pt->Draw();
  sprintf(tmpchr,"mass_spectrum_%s",namesfx);
  c1->SetName(tmpchr);
  c1->Write();

   ofstream outf("f6pi",std::ios::app);
   outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   
   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   delete datahist;
   delete sum;
   return;
}

void EEto6PI::FitSpectrum3(const TH1* hm, double beame, char* namesfx)
{
   int nBins=100;
   int Npar;
   double peakvalue = beame;
   //double beamlow=beame-0.15;
   //double beamlow=beame-0.02;
   //double beamup=beame+0.02;
   double beamlow=beame-0.15;
   double beamup=beame+0.10;
   // try to use roofit
   RooRealVar x("x","energy",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue,beamlow,beamup);
   //RooRealVar mean2("mean2","mean of gaussian",peakvalue /*-0.002*/ ,beamlow,beamup);
   RooRealVar sigma("sigma","width of gaussian",0.012,0.008,0.025);
   RooRealVar sigma2("sigma2","width of gaussian",0.03,0.02,0.05);
   //RooRealVar sigma3("sigma3","width of gaussian",0.005,0.002,0.008);
   if (strncmp(namesfx,"raw",3) == 0) mean.setVal(peakvalue - 0.003);

   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
   //RooGaussian gaus3("gaus3","gauss(x,m,s)",x,mean,sigma3);
     RooRealVar co1("co1","coefficient #1",-0.8,-100.,100.);
     RooRealVar co2("co2","coefficient #1",-0.2,-100.,100.);
     RooRealVar co3("co3","coefficient #1",0.1,-100.,100.);
     //RooRealVar co4("co4","coefficient #4",0);
     RooChebychev ground("ground","background",x,RooArgList(co1,co2,co3));
   RooRealVar signal("signal"," ",300,0,1000000);//event number
   //RooRealVar signal1("signal1"," ",20,0,1000000);//event number
   RooRealVar background("background"," ",200,0,100000);
 //RooRealVar a0("a0","coefficient #0",100,-100000,100000);
 //RooRealVar a1("a1","coefficient #1",-10,-100000,100000);
 //RooRealVar a2("a2","coefficient #2",0,-100000,100000);
   //RooRealVar a3("a3","coefficient #2",0,-100000,100000);
 //RooPolynomial ground("ground","ground",x,RooArgList(a0,a1,a2));
     
   RooRealVar alpha1("alpha1","#alpha",1.,-5,5);
   RooRealVar nnn1("n1","n",100,1,200);
   RooCBShape cbshape1("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);

   RooAddPdf *sum;
   //RooDataSet *dataset;
   RooDataHist *datahist;
   RooPlot *xframe;
   //RooDataHist *data_6pi;
 
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   TCanvas *c1 = new TCanvas("","",800,600);

   char tmpchr[100];
   sprintf(tmpchr,"data_6pi_%s",namesfx);
   //data_6pi = new RooDataHist(tmpchr,"data_6pi",x,h);
   xframe = x.frame(Title("fit 6 pi"));
   //dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   datahist = new RooDataHist(tmpchr,"data",RooArgSet(x),hm);
   //sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2),RooArgList(signal,background));
   //sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,gaus3),RooArgList(signal,background,signal1));
   sum = new RooAddPdf("sum","sum",RooArgList(cbshape1,ground),RooArgList(signal,background));
   Npar = 8;
   //sigma.setVal(0.035);
   //signal.setVal(1200);
   //background.setVal(200);
   //co1.setVal(0);
   sum->fitTo(*datahist,Range(beamlow,beamup));
   //gaus.fitTo(*datahist,Range(beame-0.012,beame+0.012));
   datahist->plotOn(xframe);
   //gaus.plotOn(xframe);
   sum->plotOn(xframe,Components(cbshape1),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(ground),LineStyle(2),LineColor(3));
   //sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
   //sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(3));
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
  sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
  pt->AddText(tmpchr);
  sprintf(tmpchr,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
  pt->AddText(tmpchr);
  sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
  pt->AddText(tmpchr);
  pt->Draw();
  sprintf(tmpchr,"mass_spectrum_%s",namesfx);
  c1->SetName(tmpchr);
  c1->Write();

   ofstream outf("f6pi",std::ios::app);
   outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   
   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   delete datahist;
   delete sum;
   return;
}

void EEto6PI::FitSpectrum4(const TH1* hm, double beame, char* namesfx)
{
   int nBins=100;
   int Npar;
   double peakvalue = beame;
   //double beamlow=beame-0.15;
   //double beamlow=beame-0.02;
   //double beamup=beame+0.02;
   double beamlow=beame-0.15;
   double beamup=beame+0.10;
   // try to use roofit
   RooRealVar x("x","energy",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue,beamlow,beamup);
   RooRealVar mean2("mean2","mean of gaussian",peakvalue /*-0.002*/ ,beamlow,beamup);
   RooRealVar sigma("sigma","width of gaussian",0.012,0.008,0.015);
   RooRealVar sigma2("sigma2","width of gaussian",0.03,0.02,0.05);
   //RooRealVar sigma3("sigma3","width of gaussian",0.005,0.002,0.008);
   if (strncmp(namesfx,"raw",3) == 0) mean.setVal(peakvalue - 0.003);

   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean2,sigma2);
   //RooGaussian gaus3("gaus3","gauss(x,m,s)",x,mean,sigma3);
     RooRealVar co1("co1","coefficient #1",-0.8,-100.,100.);
     RooRealVar co2("co2","coefficient #1",-0.2,-100.,100.);
     RooRealVar co3("co3","coefficient #1",0.1,-100.,100.);
     //RooRealVar co4("co4","coefficient #4",0);
     RooChebychev ground("ground","background",x,RooArgList(co1,co2,co3));
   RooRealVar signal("signal"," ",300,0,1000000);//event number
   //RooRealVar signal1("signal1"," ",20,0,1000000);//event number
   RooRealVar background("background"," ",200,0,100000);
 //RooRealVar a0("a0","coefficient #0",100,-100000,100000);
 //RooRealVar a1("a1","coefficient #1",-10,-100000,100000);
 //RooRealVar a2("a2","coefficient #2",0,-100000,100000);
   //RooRealVar a3("a3","coefficient #2",0,-100000,100000);
 //RooPolynomial ground("ground","ground",x,RooArgList(a0,a1,a2));
     
   RooRealVar alpha1("alpha1","#alpha",1.,-5,5);
   RooRealVar nnn1("n1","n",100,1,200);
   RooCBShape cbshape1("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);

   RooAddPdf *sum;
   //RooDataSet *dataset;
   RooDataHist *datahist;
   RooPlot *xframe;
   //RooDataHist *data_6pi;
 
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   TCanvas *c1 = new TCanvas("","",800,600);

   char tmpchr[100];
   sprintf(tmpchr,"data_6pi_%s",namesfx);
   //data_6pi = new RooDataHist(tmpchr,"data_6pi",x,h);
   xframe = x.frame(Title("fit 6 pi"));
   //dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   datahist = new RooDataHist(tmpchr,"data",RooArgSet(x),hm);
   sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2),RooArgList(signal,background));
   //sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,gaus3),RooArgList(signal,background,signal1));
   //sum = new RooAddPdf("sum","sum",RooArgList(cbshape1,ground),RooArgList(signal,background));
   Npar = 8;
   //sigma.setVal(0.035);
   //signal.setVal(1200);
   //background.setVal(200);
   //co1.setVal(0);
   sum->fitTo(*datahist,Range(beamlow,beamup));
   //gaus.fitTo(*datahist,Range(beame-0.012,beame+0.012));
   datahist->plotOn(xframe);
   //gaus.plotOn(xframe);
   //sum->plotOn(xframe,Components(cbshape1),LineStyle(2),LineColor(2));
   //sum->plotOn(xframe,Components(ground),LineStyle(2),LineColor(3));
   sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(3));
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
  sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
  pt->AddText(tmpchr);
  sprintf(tmpchr,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
  pt->AddText(tmpchr);
  sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
  pt->AddText(tmpchr);
  pt->Draw();
  sprintf(tmpchr,"mass_spectrum_%s",namesfx);
  c1->SetName(tmpchr);
  c1->Write();

   ofstream outf("f6pi",std::ios::app);
   outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   
   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   delete datahist;
   delete sum;
   return;
}





#ifdef gepep_fast6pi_cxx
gepep_fast6pi::gepep_fast6pi(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data_Rvalue_f6pi_e3850.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("data_Rvalue_f6pi_e3850.root");
      }
      f->GetObject("gepep_fast6pi",tree);

   }
   Init(tree);
}

gepep_fast6pi::~gepep_fast6pi()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gepep_fast6pi::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gepep_fast6pi::LoadTree(Long64_t entry)
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

void gepep_fast6pi::Init(TTree *tree)
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
   fChain->SetBranchAddress("ngch", &ngch, &b_ngch);
   fChain->SetBranchAddress("ncharg", &ncharg, &b_ncharg);
   fChain->SetBranchAddress("nneu", &nneu, &b_nneu);
   fChain->SetBranchAddress("pippx", pippx, &b_pippx);
   fChain->SetBranchAddress("pippy", pippy, &b_pippy);
   fChain->SetBranchAddress("pippz", pippz, &b_pippz);
   fChain->SetBranchAddress("pipe", pipe, &b_pipe);
   fChain->SetBranchAddress("pimpx", pimpx, &b_pimpx);
   fChain->SetBranchAddress("pimpy", pimpy, &b_pimpy);
   fChain->SetBranchAddress("pimpz", pimpz, &b_pimpz);
   fChain->SetBranchAddress("pime", pime, &b_pime);
   fChain->SetBranchAddress("kkm4", &kkm4, &b_kkm4);
   Notify();
}

Bool_t gepep_fast6pi::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gepep_fast6pi::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gepep_fast6pi::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif 
