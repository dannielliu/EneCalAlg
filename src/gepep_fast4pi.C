#define gepep_fast4pi_cxx
#include "gepep_fast4pi.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include "function.h"
#include "Pars.h"
#include <fstream>
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"
//#include <iostream>
extern std::string outputdir;
using RooFit::Title;
using RooFit::Components;
using RooFit::LineStyle;
using RooFit::LineColor;
using RooFit::Range;

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
   double factore,factoreerr;
   double factorpi,factorpierr;
   double me=0.000511;
   double mmu=0.105658;
   //double mpi=0.13957;
    
   ofstream ofpar;
   ofpar.open("parf4pi.txt",std::ios::app);
   ofpar<<"\n";
   ofstream detail;
   detail.open("detailf4pi.txt",std::ios::app);
   detail<<"fast4pi will give fitting results for cms energy"<<std::endl;
   ofstream purpar("par");


   //ifstream f;
   //f.open("par.txt");
   //f >> factore >> factoreerr;
   //f >> factorpi>> factorpierr;
   //f.close();
   //factorpi =1.0;
   ParMap parmap("pion.par");
   fChain->GetEntry(1);//make run a valid number
   factorpi=parmap.GetPar(run);
   double peakvalue=parmap.GetEnergy(run);// mbeam
   //factorpi=1.00134;
   std::cout<<"factor is "<<factorpi<<std::endl;
   Long64_t nentries = fChain->GetEntriesFast();
   std::cout<<"Toral entry is "<<nentries<<std::endl;
   int nBins=40;
   double mk=0.493677;
   double mpi=0.13957018;
   //double peakvalue=3.85;// cms energy
   double Dlow=peakvalue-0.07;
   double Dup=peakvalue+0.05;
   //double factor,factorlow,factorup;
   //double minimum;
   //double minx,miny;
   
   // try to use roofit
   RooRealVar x("x","energy",peakvalue,Dlow,Dup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue,Dlow,Dup);
   RooRealVar sigma1("sigma1","width of gaussian",0.0068,0.004,0.07);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma1);
   RooRealVar co1("co1","coefficient #1",0,-100000.,100000.);
   RooRealVar co2("co2","coefficient #2",0,-100000.,100000.);
   RooChebychev bkg("bkg","background",x,RooArgList(co1,co2));
   RooRealVar signal("signal"," ",120,10,1000000);//event number
   RooRealVar background("background"," ",20,0,100000);
   RooPlot *xframe;
   RooDataHist *data_f4pi;
   RooAddPdf *sum;
 
   TH1D *h = new TH1D("h","invariant mass",200,0.,4.);
   TH1D *hene = new TH1D("hene","cms energy",50,Dlow,Dup);
   TH1D *hp = new TH1D("hp","momentum of pi",80,0,4);
   
   // no correction~~~~~~~~~~~~~~
   Long64_t nbytes = 0, nb = 0;
   ofstream ofevt("evtmark.txt");
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      double totpx,totpy,totpz,tote;
      double p[4];
      double pie[4];
      double mass;
      totpx = (pippx[0]+pippx[1]+pimpx[0]+pimpx[1]);
      totpy = (pippy[0]+pippy[1]+pimpy[0]+pimpy[1]);
      totpz = (pippz[0]+pippz[1]+pimpz[0]+pimpz[1]);
      p[0]  = TMath::Sqrt(pippx[0]*pippx[0]+pippy[0]*pippy[0]+pippz[0]*pippz[0]);
      p[1]  = TMath::Sqrt(pippx[1]*pippx[1]+pippy[1]*pippy[1]+pippz[1]*pippz[1]);
      p[2]  = TMath::Sqrt(pimpx[0]*pimpx[0]+pimpy[0]*pimpy[0]+pimpz[0]*pimpz[0]);
      p[3]  = TMath::Sqrt(pimpx[1]*pimpx[1]+pimpy[1]*pimpy[1]+pimpz[1]*pimpz[1]);
      pie[0]=TMath::Sqrt(mpi*mpi+ p[0]*p[0]);
      pie[1]=TMath::Sqrt(mpi*mpi+ p[1]*p[1]);
      pie[2]=TMath::Sqrt(mpi*mpi+ p[2]*p[2]);
      pie[3]=TMath::Sqrt(mpi*mpi+ p[3]*p[3]);
      tote  = pie[0] +pie[1] +pie[2] +pie[3];
      mass = TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
      h->Fill(mass);
      if(mass>Dlow&&mass<Dup){
        hp->Fill(p[0]);
        hp->Fill(p[1]);
        hp->Fill(p[2]);
        hp->Fill(p[3]);
        ofevt<<run<<"\t"<<rec<<"\n";
      }
      hene->Fill(mass);
   }
   ofevt.close();
   
   TCanvas *c1= new TCanvas("","",800,600);
   hp->Draw();
   char tmpchr[100];
   sprintf(tmpchr,"%s/momentumpi_f4pi.eps",outputdir.c_str());
   c1->Print(tmpchr);
   
   h->Draw();
   //char tmpchr[100];
   sprintf(tmpchr,"%s/fit4pi_pre1.eps",outputdir.c_str());
   c1->Print(tmpchr);

   xframe = x.frame(Title("fit f4pi"));
   sprintf(tmpchr,"data_4pi");
   data_f4pi = new RooDataHist(tmpchr,"data_4pi",x,hene);
   sum = new RooAddPdf("sum","sum",RooArgList(gaus,bkg),RooArgList(signal,background));
   mean.setVal(peakvalue/factorpi);
   //sigma.setVal(0.035);
   signal.setVal(120);
   background.setVal(50);
   co1.setVal(0);
   sum->fitTo(*data_f4pi,Range(Dlow,Dup));
   data_f4pi->plotOn(xframe);
   sum->plotOn(xframe);
   sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(bkg),LineStyle(2),LineColor(3));
   xframe->Draw();
   sprintf(tmpchr,"%s/fit4pi_pre2.eps",outputdir.c_str());
   c1->Print(tmpchr);
   ofpar<<run<<"\n";
   ofpar<<"pre\t"<<mean.getVal()<<"\t"<<mean.getError()<<"\t"<<sigma1.getVal()<<"\t"<<sigma1.getError()<<std::endl;
   ofpar<<"\t"<<signal.getVal()<<"\t"<<signal.getError()<<"\t"<<background.getVal()<<"\t"<<background.getError();
   ofpar<<"\t"<<signal.getVal()/(signal.getVal()+background.getVal())<<std::endl;
   purpar<<mean.getVal()<<"\t"<<mean.getError()<<"\t";
   delete data_f4pi;
   delete xframe;
   delete sum;
   
   // use the factor to refit
   h->Reset();
   hene->Reset();
   nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      double totpx,totpy,totpz,tote;
      double pie[4];
      double mass;
      totpx = factorpi*(pippx[0]+pippx[1]+pimpx[0]+pimpx[1]);
      totpy = factorpi*(pippy[0]+pippy[1]+pimpy[0]+pimpy[1]);
      totpz = factorpi*(pippz[0]+pippz[1]+pimpz[0]+pimpz[1]);
      pie[0]=TMath::Sqrt(mpi*mpi+
             factorpi*factorpi*(pippx[0]*pippx[0]+pippy[0]*pippy[0]+pippz[0]*pippz[0]) );
      pie[1]=TMath::Sqrt(mpi*mpi+
             factorpi*factorpi*(pippx[1]*pippx[1]+pippy[1]*pippy[1]+pippz[1]*pippz[1]) );
      pie[2]=TMath::Sqrt(mpi*mpi+
             factorpi*factorpi*(pimpx[0]*pimpx[0]+pimpy[0]*pimpy[0]+pimpz[0]*pimpz[0]) );
      pie[3]=TMath::Sqrt(mpi*mpi+
             factorpi*factorpi*(pimpx[1]*pimpx[1]+pimpy[1]*pimpy[1]+pimpz[1]*pimpz[1]) );
      tote  = pie[0] +pie[1] +pie[2] +pie[3];
      mass = TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
      h->Fill(mass);
      hene->Fill(mass);
   }
   h->Draw();
   sprintf(tmpchr,"%s/fit4pi_re1.eps",outputdir.c_str());
   c1->Print(tmpchr);

   xframe = x.frame(Title("fit f4pi"));
   sprintf(tmpchr,"data_f4pi");
   data_f4pi = new RooDataHist(tmpchr,"data_f4pi",x,hene);
   sum = new RooAddPdf("sum","sum",RooArgList(gaus,bkg),RooArgList(signal,background));
   mean.setVal(peakvalue);
   //sigma.setVal(0.035);
   signal.setVal(120);
   background.setVal(50);
   co1.setVal(0);
   sum->fitTo(*data_f4pi,Range(Dlow,Dup));
   data_f4pi->plotOn(xframe);
   sum->plotOn(xframe);
   sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(bkg),LineStyle(2),LineColor(3));
   xframe->Draw();
   sprintf(tmpchr,"%s/fit4pi_re2.eps",outputdir.c_str());
   c1->Print(tmpchr);
   ofpar<<"re\t"<<mean.getVal()<<"\t"<<mean.getError()<<"\t"<<sigma1.getVal()<<"\t"<<sigma1.getError()<<std::endl;
   ofpar<<"\t"<<signal.getVal()<<"\t"<<signal.getError()<<"\t"<<background.getVal()<<"\t"<<background.getError();
   ofpar<<"\t"<<signal.getVal()/(signal.getVal()+background.getVal())<<std::endl;
   purpar<<mean.getVal()<<"\t"<<mean.getError()<<"\t";
   delete data_f4pi;
   delete xframe;
   delete sum;
 

   //TCanvas *c1=new TCanvas("","",800,600);
   //gStyle->SetOptStat(0);
   //gStyle->SetOptFit(1111);
   //TF1 *fit=new TF1("fit",GausLineBack,3.75,4.0);
   //fit->SetParameters(1,3.8,0.02,10,1);
   //h->Fit(fit,"","",3.75,4.0);
   //h->Draw();
   //h->Fit("gaus","","",3.75,4.0);
   //c1->Print("f4pi.eps");
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
