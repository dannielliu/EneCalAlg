#define gepep_fast4pi_cxx
#include "gepep_fast4pi.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include "function.h"
#include "EventClass.h"
#include <fstream>
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"
//#include <iostream>
extern std::string outputdir;
using namespace RooFit;

namespace PIPIKK{
  void FitSpe(std::vector<FourPi> &evts, double beame,  const char* namesfx);
  void FitSpectrum(TTree *&dataraw,double beame, const char* namesfx);
  //double GetEnergy(int runNo);
}

void gepep_fast4pi::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L gepep_fast4pi.C
//      Root > gepep_fast4pi t
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

   //TH1D *h = new TH1D("h","invariant mass",100,3.,4.);
   double mpi = 0.13957;
   double mka = 0.493677;
  
  char fname[1000];
  sprintf(fname,"%s/plot_4pi_1.root",outputdir.c_str());
  TFile *f=new TFile(fname,"RECREATE");
  
  TTree *vars = new TTree("vars","vars");
  double mass;
  //double phi,phi1,phi2;
  //double costheta,costheta1,costheta2;
  double ppi1,ppi2,ppi3,ppi4;
  double cospi1,cospi2,cospi3,cospi4;
  vars->Branch("ppi1",&ppi1,"ppi1/D");
  vars->Branch("ppi2",&ppi2,"ppi2/D");
  vars->Branch("ppi3",&ppi3,"ppi3/D");
  vars->Branch("ppi4",&ppi4,"ppi4/D");
  vars->Branch("cospi1",&cospi1,"cospi1/D");
  vars->Branch("cospi2",&cospi2,"cospi2/D");
  vars->Branch("cospi3",&cospi3,"cospi3/D");
  vars->Branch("cospi4",&cospi4,"cospi4/D");
  vars->Branch("mass",&mass,"mass/D");


   std::vector<FourPi> evts;
   std::vector<FourPi> evts2;
   //std::cout<<"factor is "<<factorpi<<std::endl;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      FourPi evt,evt2;
      HepLorentzVector pip,pim,kap,kam;
      //double deltamass1=10, deltamass2=10, tmpmass;

      pip.setVectM(Hep3Vector(pippx[0],pippy[0],pippz[0]),mpi);
      kap.setVectM(Hep3Vector(pippx[1],pippy[1],pippz[1]),mpi);
      pim.setVectM(Hep3Vector(pimpx[0],pimpy[0],pimpz[0]),mpi);
      kam.setVectM(Hep3Vector(pimpx[1],pimpy[1],pimpz[1]),mpi);
      evt.set(pip,pim,kap,kam);

/*
   // pip.setVectM(Hep3Vector(pippx[0],pippy[0],pippz[0]),mpi);
   // kap.setVectM(Hep3Vector(pippx[1],pippy[1],pippz[1]),mka);
   // pim.setVectM(Hep3Vector(pimpx[0],pimpy[0],pimpz[0]),mpi);
   // kam.setVectM(Hep3Vector(pimpx[1],pimpy[1],pimpz[1]),mka);
   // evt.set(pip,pim,kap,kam);
   // evt2.set(pip,pim,kap,kam);
   // evts2.push_back(evt2);
   // deltamass1 = evt.m() - beamene;
   
   // pip.setVectM(Hep3Vector(pippx[1],pippy[1],pippz[1]),mpi);
   // kap.setVectM(Hep3Vector(pippx[0],pippy[0],pippz[0]),mka);
   // pim.setVectM(Hep3Vector(pimpx[0],pimpy[0],pimpz[0]),mpi);
   // kam.setVectM(Hep3Vector(pimpx[1],pimpy[1],pimpz[1]),mka);
   // evt2.set(pip,pim,kap,kam);
   // evts2.push_back(evt2);
   // tmpmass = (pip+pim+kap+kam).m() - beamene;
   // if (fabs(tmpmass) < fabs(deltamass1)){
   //     evt.set(pip,pim,kap,kam);
   //     deltamass1 = tmpmass;
   // } 
   // 

   // pip.setVectM(Hep3Vector(pippx[1],pippy[1],pippz[1]),mpi);
   // kap.setVectM(Hep3Vector(pippx[0],pippy[0],pippz[0]),mka);
   // pim.setVectM(Hep3Vector(pimpx[1],pimpy[1],pimpz[1]),mpi);
   // kam.setVectM(Hep3Vector(pimpx[0],pimpy[0],pimpz[0]),mka);
   // evt2.set(pip,pim,kap,kam);
   // evts2.push_back(evt2);
   // tmpmass = (pip+pim+kap+kam).m() - beamene;
   // if (fabs(tmpmass) < fabs(deltamass1)){
   //     evt.set(pip,pim,kap,kam);
   //     deltamass1 = tmpmass;
   // } 
   // 

   // pip.setVectM(Hep3Vector(pippx[0],pippy[0],pippz[0]),mpi);
   // kap.setVectM(Hep3Vector(pippx[1],pippy[1],pippz[1]),mka);
   // pim.setVectM(Hep3Vector(pimpx[1],pimpy[1],pimpz[1]),mpi);
   // kam.setVectM(Hep3Vector(pimpx[0],pimpy[0],pimpz[0]),mka);
   // evt2.set(pip,pim,kap,kam);
   // evts2.push_back(evt2);
   // tmpmass = (pip+pim+kap+kam).m() - beamene;
   // if (fabs(tmpmass) < fabs(deltamass1)){
   //     evt.set(pip,pim,kap,kam);
   //     deltamass1 = tmpmass;
   // } 
  */    
      ppi1 = evt.pip1.rho();
      ppi2 = evt.pim1.rho();
      ppi3 = evt.pip2.rho();
      ppi4 = evt.pim2.rho();
      cospi1 = evt.pip1.cosTheta();
      cospi2 = evt.pim1.cosTheta();
      cospi3 = evt.pip2.cosTheta();
      cospi4 = evt.pim2.cosTheta();
      mass = evt.m();
      vars->Fill();
    //if (ppi3 > 1.5) continue;
    //if (ppi4 > 1.5) continue;

      evts.push_back(evt);

   }
   vars->Write();
   
   PIPIKK::FitSpe(evts,beamene,"nearest");
   //PIPIKK::FitSpe(evts2,beamene,"secondnear");

}

void PIPIKK::FitSpe(std::vector<FourPi> &evts, double beame, const char *namesfx)
{
  double beamlow= beame-0.5;
  double beamup = beame+0.5;
  // for factor fit
 
  TTree *datarawo = new TTree("datarawo","dataraw");
  TTree *dataraw = new TTree("dataraw","dataraw");
  TTree *datarawl = new TTree("datarawl","dataraw");
//TTree *datarawu = new TTree("datarawu","dataraw");
  double mass;
  datarawo->Branch("x",&mass,"x/D");
  dataraw->Branch("x",&mass,"x/D");
  datarawl->Branch("x",&mass,"x/D");
//datarawu->Branch("x",&mass,"x/D");
 
  // try to correct the spectrum
  // iniialize the fit function
  //double factor4,factor4err;// for pi
 
  //int Npart=20;
  //double pcut[Npart+1];//={0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,
		  // 0.60,0.70,0.80,0.90,1.00,1.20,1.40,1.60,1.80,2.00};//={0.0,0.5,1.0,1.5,2.0};
   //pcut[0]=0.1;
   //pcut[1]=0.4;
   //pcut[2]=0.9;
   //double facmap[Npart];
   //double facemap[Npart];

   //  set normal factor in (0.2, 0.3) to 1.00061, get factor in different range
//  only for r value data, combine both part
// factor from Ks->pipi, pi corrected with vertex fit , low p range use factor from pipill
     
  char tmpchr[100];

  //~~~~~~~~~~part start~~~~~~~~

  for (Long64_t jentry=0; jentry<evts.size();jentry++) {
     // total invariant mass
	 // without correction
     FourPi evt = evts.at(jentry);

     mass = evt.m();
     if (mass>beamlow-0.001 && mass<beamup+0.001) datarawo->Fill();
     double fpi1 = 1.000785;
     //if (evt.pip.rho()<0.4) fpi1 = 1.000902;
     //if (evt.pim.rho()<0.4) fpi2 = 1.000902;
     //if (evt.kap.rho()<0.4) fka1 = 1.000902;
     //if (evt.kam.rho()<0.4) fka2 = 1.000902;
     
  // double fpi1 = 1.00093;
  // double fpi2 = 1.00093;
  // double fka1 = 0.9992;
  // double fka2 = 0.9992;
     evt.setCorrectionFactors(fpi1);
     mass = evt.m();
     if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
   //fpi = 1.00061; fk = 1.00061  ;
   //mass = evts.at(jentry).InvMass(fpi,fpi,fk,fk);
   //if (mass>beamlow-0.001 && mass<beamup+0.001) datarawl->Fill();
  }
  //dataraw->Write();
  // no correction
  sprintf(tmpchr,"raw_%s",namesfx);
  PIPIKK::FitSpectrum(datarawo,beame,tmpchr);
  
  sprintf(tmpchr,"cor_%s",namesfx);
  PIPIKK::FitSpectrum(dataraw,beame,tmpchr);
  
  //~~~~~~~~~~ part end~~~~~~~~
  return;
}

void PIPIKK::FitSpectrum(TTree *&dataraw, double beame, const char* namesfx)
{
   int nBins=100;
   bool largesample = false;
   if (dataraw->GetEntries()>10000) largesample = true;
   int Npar;
   double peakvalue = beame;
   double beamlow   = beame - 0.15;
   double beamup    = beame + 0.15;
   // try to use roofit
   RooRealVar x("x","energy",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue-0.005,beamlow,beamup);
   //RooRealVar mean2("mean2","mean of gaussian",peakvalue-0.005,beamlow,beamup);
   RooRealVar sigma("sigma","width of gaussian",0.015,0.01,0.02);
   RooRealVar sigma2("sigma2","width of gaussian",0.025,0.02,0.05);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
     RooRealVar co1("co1","coefficient #1",0.5,-100.,100.);
     RooRealVar co2("co2","coefficient #2",0.3,-100.,100.);
     RooRealVar co3("co3","coefficient #3",0.2,-100.,100.);
     RooRealVar co4("co4","coefficient #4",0);
     RooChebychev ground("bkg","background",x,RooArgList(co1,co2,co3));
   RooRealVar signal("signal"," ",200,0,10000000);//event number
   RooRealVar signal2("signal2"," ",50,0,10000000);//event number
   RooRealVar background("background"," ",2000,0,1000000);
// RooRealVar a0("a0","coefficient #0",100,-100000,100000);
// RooRealVar a1("a1","coefficient #1",100,-100000,100000);
// RooRealVar a2("a2","coefficient #2",10,-100000,100000);
// RooRealVar a3("a3","coefficient #3",8,-100000,100000);
// RooPolynomial ground("ground","ground",x,RooArgList(a0,a1,a2,a3));
     
   RooRealVar alpha1("alpha1","#alpha",1.32,-5,5);
   RooRealVar nnn1("n1","n",100,1,200);
   RooCBShape cbshape1("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);
   
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   //RooDataHist *data_6pi;
 
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   TCanvas *c1=new TCanvas("","",800,600);
	 
   char tmpchr[100];
   sprintf(tmpchr,"data_4pi_%s",namesfx);
   //data_6pi = new RooDataHist(tmpchr,"data_6pi",x,h);
   xframe = x.frame(Title("fit kkpipi"));
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   //if (!largesample) {
     //sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,ground),RooArgList(signal,signal2,background));
     sum = new RooAddPdf("sum","sum",RooArgList(cbshape1,ground),RooArgList(signal,background));
     Npar = 9;
   //}
   //else {
     //sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,ground),RooArgList(signal,signal2,background));
     //sum = new RooAddPdf("sum","sum",RooArgList(cbshape1,gaus2,ground),RooArgList(signal,signal2,background));
     //Npar=10;
   //}
   //sigma.setVal(0.035);
   //signal.setVal(1200);
   //background.setVal(200);
   //co1.setVal(0);
   sum->fitTo(*dataset,Range(beamlow,beamup));
   //sum->fitTo(*dataset,Range(beame-0.05,beame+0.03));
   //sum->fitTo(*dataset,Range(4.22,4.28));
   dataset->plotOn(xframe);
   sum->plotOn(xframe,Components(cbshape1),LineStyle(2),LineColor(2));
   //sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(4));
   sum->plotOn(xframe,Components(ground),LineStyle(2),LineColor(3));
   //if (dataraw->GetEntries()>2000) sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(4));
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
//if (largesample){
    sprintf(tmpchr,"#sigma_{2} = %1.6f #pm %1.6f",sigma2.getVal(),sigma2.getError());
    pt->AddText(tmpchr);
//}
  sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
  pt->AddText(tmpchr);
 // if (largesample){
    //sprintf(tmpchr,"signal2 = %.2f #pm %.2f",signal2.getVal(),signal2.getError());
    //pt->AddText(tmpchr);
  //}
  sprintf(tmpchr,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
  pt->AddText(tmpchr);
  sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
  pt->AddText(tmpchr);
  pt->Draw();
  sprintf(tmpchr,"mass_spectrum_%s",namesfx);
  c1->SetName(tmpchr);
  c1->Write();

   ofstream outf("par4pi",std::ios::app);
   outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   
   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   delete dataset;
   delete sum;
   return;
}







#ifdef gepep_fast4pi_cxx
gepep_fast4pi::gepep_fast4pi(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data_Rvalue_f4pi_e3850.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("data_Rvalue_f4pi_e3850.root");
      }
      f->GetObject("gepep_fast4pi",tree);

   }
   Init(tree);
}

gepep_fast4pi::~gepep_fast4pi()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gepep_fast4pi::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gepep_fast4pi::LoadTree(Long64_t entry)
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

void gepep_fast4pi::Init(TTree *tree)
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

Bool_t gepep_fast4pi::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gepep_fast4pi::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gepep_fast4pi::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif 
