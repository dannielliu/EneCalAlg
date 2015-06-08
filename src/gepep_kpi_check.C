#define gepep_kpi_cxx
#include "gepep_kpi.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "function.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include <fstream>
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooCBShape.h"
//#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "EventClass.h"
//#include <iostream>
extern std::string outputdir;
using namespace RooFit;


void gepep_kpi::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L gepep_kpi.C
//      Root > gepep_kpi t
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
   
   double beamene ;
   fChain->GetEntry(1);
   beamene = GetEnergy(run);
   std::cout<<"Toral entry is "<<nentries<<" @"<<beamene<<" GeV."<<std::endl;
   ofstream ofresult("result_kpi.txt",std::ios::app);
   ofresult << beamene << '\t';

   int nBins=40;
   double D0low=1.835;
   double D0up=1.895;
   double mk=0.493677;
   double mpi=0.13957018;
   double peakvalue=1.86486;// mD0
   
   // try to use roofit
   RooRealVar x("x","energy",1.865,D0low,D0up,"GeV");
   RooRealVar mean("mean","mean of gaussian",1.865,D0low,D0up);
   RooRealVar sigma("sigma","width of gaussian",0.005,0.003,0.006);
   RooRealVar sigma2("sigma2","width of gaussian",0.0075,0.006,0.009);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
   
   RooRealVar alpha("alpha","#alpha",1.32,-5,5);
   RooRealVar nnn("nnn","n",5,1,200);
   RooCBShape cbshape("cbshape","cbshape",x,mean,sigma,alpha,nnn);

   RooRealVar co1("co1","coefficient #1",0,-1000.,1000.);
   RooRealVar co2("co2","coefficient #2",0,-1000.,0.1);
   RooChebychev bkg("bkg","background",x,RooArgList(co1,co2));
   
   RooRealVar a0("a0","coefficient #0",100,-100000,100000);
   RooRealVar a1("a1","coefficient #1",-1,-100000,100000);
   RooPolynomial ground("ground","ground",x,RooArgList(a0,a1));
   
   RooRealVar signal("signal"," ",12000,10,1000000);//event number
   RooRealVar signal2("signal2"," ",12000,10,1000000);//event number
   RooRealVar background("background"," ",4000,0,1000000);
   RooPlot *xframe;
   //RooDataHist *data_kpi;
   RooDataSet *dataset;
   RooAddPdf *sum;
   int Npar;
 
   ofstream ofpar;
   ofpar.open("parkpi.txt",std::ios::app);
   ofstream ofpardetail;
   ofpardetail.open("detail.txt",std::ios::app);
   ofstream purepar;
   purepar.open("par");
   char fname[100];
   sprintf(fname,"%s/plot_kpi.root",outputdir.c_str());
   TFile *f = new TFile(fname,"RECREATE");

   TF1 *facfit = new TF1("facfit",line2,0.9,1.1,2);
   TH1D *h1    = new TH1D("h1","k- pi+ invariant mass",nBins,D0low,D0up);
   TCanvas *c1 = new TCanvas("","",800,600);
   TTree *dataraw = new TTree("dataraw","dataraw");
   double mass;
   dataraw->Branch("x",&mass,"x/D");

   Long64_t nbytes = 0, nb = 0;

   TH2D *h2p = new TH2D("h2p","#pi K momentum",200,0,2,200,0,2);
   int Npart=1;

   double m0=peakvalue;
   double sigma_m=0.0068;//0.0024 for phi,
   double width = 16.*sigma_m;
   //double mparticle=mk;

   char name[100];
   // ~~~~~~~~ draw nxn histogram, m distribution in different p range
   
   std::vector<D0KPI> evt_set;
   h2p->Reset();
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     
     D0KPI evt;
     HepLorentzVector kam, pip;
     if (nkam ==1)
        kam.setVectM(Hep3Vector(kampx[0],kampy[0],kampz[0]),mk);
     else continue;

     int besidx=0;
     double tmpdeltaold=100;
     double tmpmass;
     if (npip<1) continue;
     for(int i=0; i<npip;i++){
	   pip.setVectM(Hep3Vector(pippx[i],pippy[i],pippz[i]),mpi);
       tmpmass = (kam+pip).m();
       if(fabs(tmpmass-peakvalue)<tmpdeltaold){
         tmpdeltaold = fabs(tmpmass-peakvalue);
         besidx = i;
       }
     }
	 pip.setVectM(Hep3Vector(pippx[besidx],pippy[besidx],pippz[besidx]),mpi);
     mass = (kam+pip).m();
     
     // total invariant mass, D0 -> k- pi+
     //p1 = CalMom(pippx[besidx],pippy[besidx],pippz[besidx]);
	 //p2 = CalMom(kampx[0],kampy[0],kampz[0]);
     if (mass>m0-width/2. && mass<m0+width/2.)
     {
	   D0KPI evt;
	   evt.set(kam,pip);
	   evt_set.push_back(evt);
     }
     // if (Cut(ientry) < 0) continue;
   }
   

   // for saving the fit result
 
	  xframe = x.frame(Title("fit k pi"));

      h1->Reset();
	  dataraw->Reset();
      for (Long64_t jentry=0; jentry<evt_set.size();jentry++) {
		 mass = evt_set.at(jentry).m();
		 if (mass>D0low && mass<D0up){
		   dataraw->Fill();
		   h1->Fill(mass);
		 }
      }

	  char tmpchr[100];
	  sprintf(tmpchr,"data_kpi");
      dataset = new RooDataSet(tmpchr,"dataset",dataraw,x);
	  sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,ground),RooArgList(signal,signal2,background));
      mean.setVal(peakvalue);
	  Npar=8;
	  //sigma.setVal(0.035);
	  //signal.setVal(1200);
	  //background.setVal(500);
	  //co1.setVal(0);
	  //co2.setVal(-0.2);
      sum->fitTo(*dataset,Range(D0low,D0up));
	  dataset->plotOn(xframe);
	  sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
	  sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
	  sum->plotOn(xframe,Components(bkg),LineStyle(2),LineColor(3));
	  sum->plotOn(xframe);
      xframe->Draw();
	  TPaveText *pt = new TPaveText(0.12,0.50,0.5,0.90,"BRNDC");
    //pt->SetBorderSize(0);
    //pt->SetFillStyle(4000);
    //pt->SetTextAlign(12);
    //pt->SetTextFont(42);
    //pt->SetTextSize(0.035);
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
      sprintf(tmpchr,"#chi^{2}/(100-%d) = %5.6f",Npar,xframe->chiSquare(Npar));
      pt->AddText(tmpchr);
      pt->Draw();
      c1->Write("mSpec_raw");
	  delete dataset;
	  delete xframe;
	  delete sum;

      ofresult<<mean.getVal()<<'\t'<<mean.getError()<<'\t';

   // draw the best fit
   double fk  = 1.000257;
   double fpi = 1.000769;
   xframe = x.frame(Title("fit k pi"));
   h1->Reset();
   dataraw->Reset();
   //std::cout<<"factor is "<<factor<<std::endl;
   for (Long64_t jentry=0; jentry<evt_set.size();jentry++) {
	 D0KPI evt = evt_set.at(jentry);
	 double p1 = evt.pip.rho();
	 p1 < 0.4 ? fpi = 1.000902: fpi = 1.000769;
	 evt.setCorrectionFactors(fk,fpi);
	 mass = evt.m();
	 if (mass>D0low && mass<D0up){
	   dataraw->Fill();
	 }
     // if (Cut(ientry) < 0) continue;
   }
	//char tmpchr[100];
	sprintf(tmpchr,"data_kpi");
    dataset = new RooDataSet(tmpchr,"dataset",dataraw,x);
    sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,bkg),RooArgList(signal,signal2,background));
    mean.setVal(peakvalue);
	Npar=8;
	//sigma.setVal(0.035);
	//signal.setVal(1200);
	//background.setVal(200);
	//co1.setVal(0);
	//co2.setVal(-0.2);
    sum->fitTo(*dataset,Range(D0low,D0up));
	dataset->plotOn(xframe);
	sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
	sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
	sum->plotOn(xframe,Components(bkg),LineStyle(2),LineColor(3));
	sum->plotOn(xframe);
    xframe->Draw();
    //TPaveText *
	pt = new TPaveText(0.12,0.50,0.5,0.90,"BRNDC");
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
    sprintf(tmpchr,"#chi^{2}/(100-%d) = %5.6f",Npar,xframe->chiSquare(Npar));
    pt->AddText(tmpchr);
    pt->Draw();
	c1->Write("mSpec_cor");
	//xframe->Write(name);
	delete dataset;
	delete xframe;
	delete sum;
    
	ofresult<<mean.getVal()<<'\t'<<mean.getError()<<'\t';
   
   ofresult<<std::endl;
   f->Close();

}

#ifdef gepep_kpi_cxx
gepep_kpi::gepep_kpi(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data/RValue_kpi_3850.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("data/RValue_kpi_3850.root");
      }
      f->GetObject("gepep_kpi",tree);

   }
   Init(tree);
}

gepep_kpi::~gepep_kpi()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gepep_kpi::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gepep_kpi::LoadTree(Long64_t entry)
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

void gepep_kpi::Init(TTree *tree)
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
   fChain->SetBranchAddress("nkap", &nkap, &b_nkap);
   fChain->SetBranchAddress("kappx", kappx, &b_kappx);
   fChain->SetBranchAddress("kappy", kappy, &b_kappy);
   fChain->SetBranchAddress("kappz", kappz, &b_kappz);
   fChain->SetBranchAddress("kape", kape, &b_kape);
   fChain->SetBranchAddress("nkam", &nkam, &b_nkam);
   fChain->SetBranchAddress("kampx", kampx, &b_kampx);
   fChain->SetBranchAddress("kampy", kampy, &b_kampy);
   fChain->SetBranchAddress("kampz", kampz, &b_kampz);
   fChain->SetBranchAddress("kame", kame, &b_kame);
   fChain->SetBranchAddress("npip", &npip, &b_npip);
   fChain->SetBranchAddress("pippx", pippx, &b_pippx);
   fChain->SetBranchAddress("pippy", pippy, &b_pippy);
   fChain->SetBranchAddress("pippz", pippz, &b_pippz);
   fChain->SetBranchAddress("pipe", pipe, &b_pipe);
   fChain->SetBranchAddress("npim", &npim, &b_npim);
   fChain->SetBranchAddress("pimpx", pimpx, &b_pimpx);
   fChain->SetBranchAddress("pimpy", pimpy, &b_pimpy);
   fChain->SetBranchAddress("pimpz", pimpz, &b_pimpz);
   fChain->SetBranchAddress("pime", pime, &b_pime);
   fChain->SetBranchAddress("kkm4", &kkm4, &b_kkm4);
   Notify();
}

Bool_t gepep_kpi::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gepep_kpi::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gepep_kpi::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gepep_kpi_cxx
