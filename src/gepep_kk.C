#define gepep_kk_cxx
#include "gepep_kk.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "function.h"
#include "TF1.h"
#include "TGraphErrors.h"
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


void gepep_kk::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L gepep_kk.C
//      Root > gepep_kk t
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

   std::cout<<"Toral entry is "<<nentries<<std::endl;
   int nBins=50;
   double factorstart=0.99;
   double philow=1.0;
   double phiup=1.05;
   double mk=0.493677;
   double peakvalue=1.019455;// mphi
   int pointNo=20;
   double factors[pointNo];
   double factorserr[pointNo];
   double deltapeaks[pointNo];
   double deltapeakserr[pointNo];
   double factor=factorstart;
   double factorstep=(1.-factor)*2/pointNo;
   
   // try to use roofit
   RooRealVar x("x","energy",1.020,philow,phiup,"GeV");
   RooRealVar mean("mean","mean of gaussian",1.020,philow,phiup);
   RooRealVar sigma("sigma","width of gaussian",0.003,0.001,0.005);
   RooRealVar sigma2("sigma2","width of gaussian",0.02,0.005,0.05);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
   RooRealVar co1("co1","coefficient #1",0,-1000.,1000.);
   //RooRealVar co4("co4","coefficient #4",0);
   RooChebychev bkg("bkg","background",x,RooArgList(co1));
   RooRealVar signal("signal"," ",1200,10,1000000);//event number
   RooRealVar signal2("signal2"," ",1200,10,1000000);//event number
   RooRealVar background("background"," ",200,0,100000);
   RooPlot *xframe;
   RooDataHist *data_k;
   RooAddPdf *sum;
   
   ofstream ofpar;
   ofpar.open("parkk.txt",std::ios::app);
   ofpar<<"kk algrithm: will give factors for kaon"<<std::endl;
   ofstream ofpardetail;
   ofpardetail.open("detail.txt",std::ios::app);

   TF1 *facfit = new TF1("facfit",line2,0.9,1.1,2);
   TH1D *h1   = new TH1D("h1","2 kaon invariant mass",nBins,philow,phiup);
   TCanvas *c1= new TCanvas("","",800,600);

   // for saving the fit result
   std::string fitepsname  = outputdir+"/fitkk.eps";
   std::string fiteps_start=fitepsname+"[";
   std::string fiteps_stop =fitepsname+"]";
   c1->Print(fiteps_start.c_str());
   
   Long64_t nbytes = 0, nb = 0;
   int fittimes = 0;
   for (int i=0;i<pointNo;i++){
	  //xframe->Clear();
	  xframe = x.frame(Title("fit kaon"));

      h1->Reset();
      std::cout<<"factor is "<<factor<<std::endl;
      for (Long64_t jentry=0; jentry<nentries;jentry++) {
         Long64_t ientry = LoadTree(jentry);
         if (ientry < 0) break;
         nb = fChain->GetEntry(jentry);   nbytes += nb;
         
		 //if(ngam>0) continue;
	     double mass;
	     double totpx,totpy,totpz,tote;
		 double ke[2];
	     // total invariant mass
	     totpx=factor*(kappx+kampx);
	     totpy=factor*(kappy+kampy);
	     totpz=factor*(kappz+kampz);
         ke[0]=TMath::Sqrt(mk*mk + 
		       factor*factor*(kappx*kappx+kappy*kappy+kappz*kappz));
		 ke[1]=TMath::Sqrt(mk*mk + 
		       factor*factor*(kampx*kampx+kampy*kampy+kampz*kampz));
	     tote=ke[0]+ke[1];
	     mass=TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
	     h1->Fill(mass);
         // if (Cut(ientry) < 0) continue;
      }

	  char tmpchr[100];
	  sprintf(tmpchr,"data_k_%02d",fittimes);
      data_k = new RooDataHist(tmpchr,"data_k",x,h1);
      sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,bkg),RooArgList(signal,signal2,background));
      mean.setVal(peakvalue+0.05*(factor-1.0));
	  //sigma.setVal(0.035);
	  signal.setVal(12000);
	  background.setVal(2000);
	  co1.setVal(0);
      sum->fitTo(*data_k,Range(philow,phiup));
	  data_k->plotOn(xframe);
	  sum->plotOn(xframe);
	  sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
	  sum->plotOn(xframe,Components(gaus2),LineStyle(4),LineColor(4));
	  sum->plotOn(xframe,Components(bkg),LineStyle(3),LineColor(3));
      xframe->Draw();
	  c1->Print(fitepsname.c_str());
	  delete data_k;
	  delete xframe;
	  delete sum;

	  //sprintf(tmpchr,"data_k_%d.eps",fittimes);
	  //h1->Draw();
	  //c1->Print(tmpchr);
	  
	  // save pars
	  factors[i]=factor;
	  factorserr[i]=0;
	  deltapeaks[i] = mean.getValV() - peakvalue;
	  deltapeakserr[i] = mean.getError();

	  fittimes++;
	  factor += factorstep;
   }
   std::cout<<"entry is "<<nentries<<std::endl;
   c1->Print(fiteps_stop.c_str());
   c1->Clear();
   
   TGraphErrors *graph1 = new TGraphErrors(pointNo,factors,deltapeaks,factorserr,deltapeakserr);
   graph1->SetTitle("delta peak");
   graph1->SetMarkerStyle(5);
   graph1->Draw("AP");
   gStyle->SetOptFit(1111);
   facfit->SetParameters(1,0.3);
   facfit->SetParNames("factor","slope");
   graph1->Fit(facfit,"","",factors[0],factors[pointNo-1]);
   //factor1=facfit->GetParameter(0);
   //factor1err=facfit->GetParError(0);
   ofpar<<facfit->GetParameter(0)<<"\t"<<facfit->GetParError(0)<<std::endl;
   ofpar<<signal.getValV()<<"\t"<<signal.getError()<<std::endl;
   //std::cout<<"fit factor: "<<factor1<<", error is "<<factor1err<<std::endl;
   std::string tmpstr=outputdir+"/factork.eps";
   c1->Print(tmpstr.c_str());

   // draw the best fitting
	  xframe = x.frame(Title("fit kaon"));
      h1->Reset();
      std::cout<<"factor is "<<factor<<std::endl;
      for (Long64_t jentry=0; jentry<nentries;jentry++) {
         Long64_t ientry = LoadTree(jentry);
         if (ientry < 0) break;
         nb = fChain->GetEntry(jentry);   nbytes += nb;
         
		 //if(ngam>0) continue;
	     double mass;
	     double totpx,totpy,totpz,tote;
		 double ke[2];
	     // total invariant mass
	     totpx=factor*(kappx+kampx);
	     totpy=factor*(kappy+kampy);
	     totpz=factor*(kappz+kampz);
         ke[0]=TMath::Sqrt(mk*mk + 
		       factor*factor*(kappx*kappx+kappy*kappy+kappz*kappz));
		 ke[1]=TMath::Sqrt(mk*mk + 
		       factor*factor*(kampx*kampx+kampy*kampy+kampz*kampz));
	     tote=ke[0]+ke[1];
	     mass=TMath::Sqrt(tote*tote-totpx*totpx-totpy*totpy-totpz*totpz);
	     h1->Fill(mass);
         // if (Cut(ientry) < 0) continue;
      }

	  char tmpchr[100];
	  sprintf(tmpchr,"data_k_%02d",fittimes);
      data_k = new RooDataHist(tmpchr,"data_k",x,h1);
      sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,bkg),RooArgList(signal,signal2,background));
      mean.setVal(peakvalue+0.05*(factor-1.0));
	  //sigma.setVal(0.035);
	  signal.setVal(12000);
	  background.setVal(2000);
	  co1.setVal(0);
      sum->fitTo(*data_k,Range(philow,phiup));
	  data_k->plotOn(xframe);
	  sum->plotOn(xframe);
	  sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
	  sum->plotOn(xframe,Components(gaus2),LineStyle(4),LineColor(4));
	  sum->plotOn(xframe,Components(bkg),LineStyle(3),LineColor(3));
      xframe->Draw();
	  c1->Print("fitkk_best.eps");
	  delete data_k;
	  delete xframe;
	  delete sum;

  
   ofpar.close();
   ofpardetail.close();

}

#ifdef gepep_kk_cxx
gepep_kk::gepep_kk(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data/RValue_kk_3850.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("data/RValue_kk_3850.root");
      }
      f->GetObject("gepep_kk",tree);

   }
   Init(tree);
}

gepep_kk::~gepep_kk()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gepep_kk::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gepep_kk::LoadTree(Long64_t entry)
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

void gepep_kk::Init(TTree *tree)
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
   fChain->SetBranchAddress("kappx", &kappx, &b_kappx);
   fChain->SetBranchAddress("kappy", &kappy, &b_kappy);
   fChain->SetBranchAddress("kappz", &kappz, &b_kappz);
   fChain->SetBranchAddress("kape", &kape, &b_kape);
   fChain->SetBranchAddress("kampx", &kampx, &b_kampx);
   fChain->SetBranchAddress("kampy", &kampy, &b_kampy);
   fChain->SetBranchAddress("kampz", &kampz, &b_kampz);
   fChain->SetBranchAddress("kame", &kame, &b_kame);
   fChain->SetBranchAddress("mphi", &mphi, &b_mphi);
   fChain->SetBranchAddress("kkm4", &kkm4, &b_kkm4);
   Notify();
}

Bool_t gepep_kk::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gepep_kk::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gepep_kk::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gepep_kk_cxx
