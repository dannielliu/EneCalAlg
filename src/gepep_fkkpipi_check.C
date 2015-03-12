#define gepep_fkkpipi_cxx
#include "gepep_fkkpipi.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "function.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include <fstream>
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
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
namespace KKPIPI{
  void FitSpe(std::vector<KKpipi> evts, double beame,  char* namesfx);
  void FitSpectrum(TTree *&dataraw,double beame, char* namesfx);
  double GetEnergy(int runNo);
}

void gepep_fkkpipi::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L gepep_fkkpipi.C
//      Root > gepep_fkkpipi t
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
//
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();

   fChain->GetEntry(1);
   double beamene = KKPIPI::GetEnergy(run);
   beamene = 4.26;
   std::cout<<"current beam energy is "<<beamene<<", run id "<<run<<std::endl;
   if (beamene < 0.1){
     std::cout<<"can not get a suitable beam energy!"<<std::endl;
	 return;
   }

   std::cout<<"Toral entry is "<<nentries<<std::endl;
   int nBins=100;
   double peakvalue=beamene;// mbeam
   double beamlow=beamene-0.1;
   double beamup=beamene+0.1;
   double mpi=0.13957;
   double mk = 0.493677;
    
  char fname[1000];
  sprintf(fname,"%s/plot_kkpipi.root",outputdir.c_str());
  TFile *f=new TFile(fname,"RECREATE");

  char name[100];
  TCanvas *c1=new TCanvas("c1","",800,600);
  TTree *vars = new TTree("vars","vars");
  double mass;
  //double phi,phi1,phi2;
  //double costheta,costheta1,costheta2;
  double ppi1,ppi2,pk1,pk2;
  vars->Branch("ppi1",&ppi1,"ppi1/D");
  vars->Branch("ppi2",&ppi2,"ppi2/D");
  vars->Branch("pk1",&pk1,"pk1/D");
  vars->Branch("pk2",&pk2,"pk2/D");
  vars->Branch("mass",&mass,"mass/D");
  double chisq;
  vars->Branch("chi2",&chisq,"chi2/D");
   KKpipi evt;
   std::vector<KKpipi> evts;

   // select useful events
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;	  
	  
	  evt.SetVal(pippx,pippy,pippz,pimpx,pimpy,pimpz,
	             kappx,kappy,kappz,kampx,kampy,kampz);
      mass = evt.InvMass();
	  ppi1 = evt.GetP1pi();
	  ppi2 = evt.GetP2pi();
	  pk1 = evt.GetP1k();
	  pk2 = evt.GetP2k();
	  chisq = kkm4;
	  if (mass>beamlow && mass<beamup){
	    vars->Fill();
		evts.push_back(evt);
	  }
   }//select end
   vars->Write();
   sprintf(name,"%f",peakvalue);
   KKPIPI::FitSpe(evts,peakvalue,name);
   return;

}

void KKPIPI::FitSpe(std::vector<KKpipi> evts, double beame, char *namesfx)
{
  double beamlow=beame-0.1;
  double beamup=beame+0.1;
  // for factor fit
 
  TTree *datarawo = new TTree("datarawo","dataraw");
//TTree *dataraw = new TTree("dataraw","dataraw");
//TTree *datarawl = new TTree("datarawl","dataraw");
//TTree *datarawu = new TTree("datarawu","dataraw");
  double mass;
  datarawo->Branch("x",&mass,"x/D");
//dataraw->Branch("x",&mass,"x/D");
//datarawl->Branch("x",&mass,"x/D");
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
     mass = evts.at(jentry).InvMass();
     if (mass>beamlow-0.001 && mass<beamup+0.001) datarawo->Fill();
  }
  //dataraw->Write();
  // no correction
  sprintf(tmpchr,"raw_%s",namesfx);
  KKPIPI::FitSpectrum(datarawo,beame,tmpchr);
  
  //~~~~~~~~~~ part end~~~~~~~~
  return;
}

void KKPIPI::FitSpectrum(TTree *&dataraw, double beame, char* namesfx)
{
   int nBins=100;
   bool largesample = false;
   if (dataraw->GetEntries()>10000) largesample = true;
   int Npar;
   double peakvalue = beame;
   double beamlow=beame-0.1;
   double beamup=beame+0.1;
   // try to use roofit
   RooRealVar x("x","energy",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue,beamlow,beamup);
   RooRealVar sigma("sigma","width of gaussian",0.003,0.0001,0.02);
   RooRealVar sigma2("sigma2","width of gaussian",0.022,0.02,0.05);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
   //RooRealVar co1("co1","coefficient #1",0,-100.,100.);
   //RooRealVar co4("co4","coefficient #4",0);
   //RooChebychev bkg("bkg","background",x,RooArgList(co1));
   RooRealVar signal("signal"," ",200,0,10000000);//event number
   RooRealVar signal2("signal2"," ",200,0,10000000);//event number
   RooRealVar background("background"," ",20,0,1000000);
   RooRealVar a0("a0","coefficient #0",100,-100000,100000);
   RooRealVar a1("a1","coefficient #1",-1,-100000,100000);
   RooPolynomial ground("ground","ground",x,RooArgList(a0,a1));
   
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   //RooDataHist *data_6pi;
 
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   TCanvas *c1=new TCanvas("","",800,600);
	 
   char tmpchr[100];
   sprintf(tmpchr,"data_kkpipi_%s",namesfx);
   //data_6pi = new RooDataHist(tmpchr,"data_6pi",x,h);
   xframe = x.frame(Title("fit kkpipi"));
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   if (!largesample) {
     sum = new RooAddPdf("sum","sum",RooArgList(gaus,ground),RooArgList(signal,background));
     Npar = 6;
   }
   else {
     sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,ground),RooArgList(signal,signal2,background));
	 Npar=8;
   }
   //sigma.setVal(0.035);
   //signal.setVal(1200);
   //background.setVal(200);
   //co1.setVal(0);
   sum->fitTo(*dataset,Range(beamlow,beamup));
   //sum->fitTo(*dataset,Range(beame-0.05,beame+0.03));
   //sum->fitTo(*dataset,Range(4.22,4.28));
   dataset->plotOn(xframe);
   sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(ground),LineStyle(2),LineColor(3));
   if (dataraw->GetEntries()>2000) sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(4));
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
  if (largesample){
    sprintf(tmpchr,"#sigma_{2} = %1.6f #pm %1.6f",sigma2.getVal(),sigma2.getError());
    pt->AddText(tmpchr);
  }
  sprintf(tmpchr,"signal1 = %.2f #pm %.2f",signal.getVal(),signal.getError());
  pt->AddText(tmpchr);
  if (largesample){
    sprintf(tmpchr,"signal2 = %.2f #pm %.2f",signal2.getVal(),signal2.getError());
    pt->AddText(tmpchr);
  }
  sprintf(tmpchr,"backNo = %.2f #pm %.2f",background.getVal(),background.getError());
  pt->AddText(tmpchr);
  sprintf(tmpchr,"#chi^{2}/(%d-%d) = %5.6f",nBins,Npar,xframe->chiSquare(Npar));
  pt->AddText(tmpchr);
  pt->Draw();
  sprintf(tmpchr,"mass_spectrum_%s",namesfx);
  c1->SetName(tmpchr);
  c1->Write();

   ofstream outf("fkkpipi",std::ios::app);
   outf<<beame<<"\t"<<namesfx<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   
   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   delete dataset;
   delete sum;
   return;
}


double KKPIPI::GetEnergy(int runNo)
{
  runNo=runNo%10000;
  //int a[2] = {1, 2};
  int runNoLow[108] = {
   4011, 4028, 4037, 4046, 4058, 4069, 4077, 4084, 4091, 4097,
   4105, 4118, 4128, 4135, 4142, 4151, 4161, 4175, 4184, 4191,
   4197, 4203, 4211, 4221, 4231, 4240, 4246, 4253, 4258, 4266,
   4272, 4277, 4282, 4298, 4314, 4321, 4328, 4339, 4346, 4351,
   4359, 4369, 4374, 4382, 4390, 4397, 4404, 4412, 4418, 4428,
   4437, 4447, 4461, 4478, 4486, 4494, 4503, 4512, 4527, 4541,
   4555, 4564, 4574, 4585, 4593, 4603, 4613, 4623, 4634, 4642,
   4652, 4661, 4674, 4685, 4695, 4705, 4719, 4729, 4740, 4754,
   4763, 4777, 4785, 4794, 4804, 4812, 4825, 4837, 4848, 4861,
   4869, 4882, 4891, 4900, 4913, 4926, 4936, 4947, 4958, 4968,
   4982, 5010, 5027, 5041, 5060, 5082, 5099, 5119 };// #the last number is not used!
  double energy[108] = {
   3850, 3890, 3895, 3900, 3905, 3910, 3915, 3920, 3925, 3930,
   3935, 3940, 3945, 3950, 3955, 3960, 3965, 3970, 3975, 3980,
   3985, 3990, 3995, 4000, 4005, 4010, 4012, 4014, 4016, 4018,
   4020, 4025, 4030, 0000, 4035, 4040, 4050, 4055, 4060, 4065,
   4070, 4080, 4090, 4100, 4110, 4120, 4130, 4140, 4145, 4150,
   4160, 4170, 0000, 4180, 4190, 4195, 4200, 4203, 4206, 4210,
   4215, 4220, 4225, 4230, 4235, 4240, 4243, 4245, 4248, 4250,
   4255, 4260, 4265, 4270, 4275, 4280, 4285, 4290, 4300, 4310,
   4320, 4330, 4340, 4350, 4360, 4370, 4380, 4390, 4395, 4400,
   4410, 4420, 4425, 4430, 4440, 4450, 4460, 4480, 4500, 4520,
   4540, 4550, 4560, 4570, 4580, 0000, 4590, 4600 };// # the last one is skipped
  for(int i=0;i<107;i++){
    if(runNo>=runNoLow[i]&&runNo<runNoLow[i+1])
      return energy[i]/1000;
  }
  return -1;
}








#ifdef gepep_fkkpipi_cxx
gepep_fkkpipi::gepep_fkkpipi(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data/RValue_fkkpipi_3850.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("data/RValue_fkkpipi_3850.root");
      }
      f->GetObject("gepep_fastkkpipi",tree);

   }
   Init(tree);
}

gepep_fkkpipi::~gepep_fkkpipi()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gepep_fkkpipi::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gepep_fkkpipi::LoadTree(Long64_t entry)
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

void gepep_fkkpipi::Init(TTree *tree)
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
   fChain->SetBranchAddress("pippx", &pippx, &b_pippx);
   fChain->SetBranchAddress("pippy", &pippy, &b_pippy);
   fChain->SetBranchAddress("pippz", &pippz, &b_pippz);
   fChain->SetBranchAddress("pipe", &pipe, &b_pipe);
   fChain->SetBranchAddress("pimpx", &pimpx, &b_pimpx);
   fChain->SetBranchAddress("pimpy", &pimpy, &b_pimpy);
   fChain->SetBranchAddress("pimpz", &pimpz, &b_pimpz);
   fChain->SetBranchAddress("pime", &pime, &b_pime);
   fChain->SetBranchAddress("kappx", &kappx, &b_kappx);
   fChain->SetBranchAddress("kappy", &kappy, &b_kappy);
   fChain->SetBranchAddress("kappz", &kappz, &b_kappz);
   fChain->SetBranchAddress("kape", &kape, &b_kape);
   fChain->SetBranchAddress("kampx", &kampx, &b_kampx);
   fChain->SetBranchAddress("kampy", &kampy, &b_kampy);
   fChain->SetBranchAddress("kampz", &kampz, &b_kampz);
   fChain->SetBranchAddress("kame", &kame, &b_kame);
   fChain->SetBranchAddress("kkm4", &kkm4, &b_kkm4);
   Notify();
}

Bool_t gepep_fkkpipi::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gepep_fkkpipi::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gepep_fkkpipi::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gepep_fkkpipi_cxx
