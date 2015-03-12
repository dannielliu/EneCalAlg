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
//#include "RooChebychev.h"
//#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooMsgService.h"
#include <vector>
extern std::string outputdir;
using namespace RooFit;
namespace EEto6PI{
  void FitSpe(std::vector<EEto6pi> evts, double beame,  char* namesfx=0);
  void FitSpectrum(TTree *&dataraw, double beame, char* namesfx=0);
  double GetEnergy(int runNo);
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
   double beamene = EEto6PI::GetEnergy(run);
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
    
  char fname[1000];
  sprintf(fname,"%s/plot_6pi_fKsnovtx_4260.root",outputdir.c_str());
  TFile *f=new TFile(fname,"RECREATE");

  char name[100];
  TCanvas *c1=new TCanvas("c1","",800,600);
  TTree *vars = new TTree("vars","vars");
  double mass;
  //double phi,phi1,phi2;
  //double costheta,costheta1,costheta2;
  double p[6];
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
	  //if (run > 30367) break;
	  evt.Setval(pippx,pippy,pippz,pimpx,pimpy,pimpz);
      mass = evt.InvMass();
	  bool good=true;;
      for (int i=0; i<6;i++) {
	    p[i] = evt.GetP(i);
		if (p[i]<0.1 || p[i]>1.0) good = false;
	  }
	  if (mass>beamlow && mass<beamup){
	    vars->Fill();
	    if (!good) continue;
		evts.push_back(evt);
	  }
   }//select end
   vars->Write();
   sprintf(name,"%f",peakvalue);
   EEto6PI::FitSpe(evts,peakvalue,name);
   return;
}

void EEto6PI::FitSpe(std::vector<EEto6pi> evts, double beame, char *namesfx)
{
  double beamlow=beame-0.1;
  double beamup=beame+0.1;
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
 
  int Npart=20;
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
/*
   pcut[0] =0.0 ;  facmap[0] =1.0;       facemap[0] =1.0;     
   pcut[1] =0.05;  facmap[1] =1.00594 ;  facemap[1] =0.00146173 ;
   pcut[2] =0.10;  facmap[2] =1.00345 ;  facemap[2] =0.000209543;
   pcut[3] =0.15;  facmap[3] =1.00191 ;  facemap[3] =9.96687e-05;
   pcut[4] =0.20;  facmap[4] =1.00074 ;  facemap[4] =7.83398e-05;
   pcut[5] =0.25;  facmap[5] =0.99993 ;  facemap[5] =6.57175e-05;
   pcut[6] =0.30;  facmap[6] =0.999425;  facemap[6] =6.42244e-05;
   pcut[7] =0.35;  facmap[7] =0.999187;  facemap[7] =6.86471e-05;
   pcut[8] =0.40;  facmap[8] =0.998956;  facemap[8] =7.52353e-05;
   pcut[9] =0.45;  facmap[9] =0.998686;  facemap[9] =8.37434e-05;
   pcut[10]=0.50;  facmap[10]=0.998865;  facemap[10]=6.97743e-05;
   pcut[11]=0.60;  facmap[11]=0.99871 ;  facemap[11]=8.90833e-05;
   pcut[12]=0.70;  facmap[12]=0.998473;  facemap[12]=0.000118579;
   pcut[13]=0.80;  facmap[13]=0.99859 ;  facemap[13]=0.000153497;
   pcut[14]=0.90;  facmap[14]=0.998797;  facemap[14]=0.000208148;
   pcut[15]=1.00;  facmap[15]=0.999028;  facemap[15]=0.000232254;
   pcut[16]=1.20;  facmap[16]=0.999368;  facemap[16]=0.000521802;
   pcut[17]=1.40;  facmap[17]=0.999502;  facemap[17]=0.00190491 ;
   pcut[18]=1.60;  facmap[18]=1.0;       facemap[18]=1.0;     
   pcut[19]=1.80;  facmap[19]=1.0;       facemap[19]=1.0;     
   pcut[20]=2.00;  //facmap[20]=2.00;
 */

   pcut[0] =0.0 ;  facmap[0] =1.0;       facemap[0] =1.0;     
   pcut[1] =0.05;  facmap[1] =1.0103  ;  facemap[1] =0.00271272 ;
   pcut[2] =0.10;  facmap[2] =1.00242 ;  facemap[2] =0.000294394;
   pcut[3] =0.15;  facmap[3] =1.00154 ;  facemap[3] =0.000168453;
   pcut[4] =0.20;  facmap[4] =1.00072 ;  facemap[4] =0.000152402;
   pcut[5] =0.25;  facmap[5] =1.00036 ;  facemap[5] =0.000135882;
   pcut[6] =0.30;  facmap[6] =0.999792;  facemap[6] =0.000146793;
   pcut[7] =0.35;  facmap[7] =1       ;  facemap[7] =0.000178329;
   pcut[8] =0.40;  facmap[8] =0.999605;  facemap[8] =0.000256467;
   pcut[9] =0.45;  facmap[9] =0.999899;  facemap[9] =0.000262561;
   pcut[10]=0.50;  facmap[10]=0.999495;  facemap[10]=0.000247957;
   pcut[11]=0.60;  facmap[11]=1.00016 ;  facemap[11]=0.000364348;
   pcut[12]=0.70;  facmap[12]=0.999504;  facemap[12]=0.000537846;
   pcut[13]=0.80;  facmap[13]=1.00029 ;  facemap[13]=0.00077103 ;
   pcut[14]=0.90;  facmap[14]=1.00031 ;  facemap[14]=0.00123488 ;
   pcut[15]=1.00;  facmap[15]=1.00313 ;  facemap[15]=0.00150405 ;
   pcut[16]=1.20;  facmap[16]=0.986371;  facemap[16]=0.00626053 ;
   pcut[17]=1.40;  facmap[17]=0.938064;  facemap[17]=0.0999831  ;
   pcut[18]=1.60;  facmap[18]=1.0;       facemap[18]=1.0;     
   pcut[19]=1.80;  facmap[19]=1.0;       facemap[19]=1.0;     
   pcut[20]=2.00;  //facmap[20]=2.00;
 




  // psi 3.097
     
  char tmpchr[100];

  //~~~~~~~~~~part start~~~~~~~~

  for (Long64_t jentry=0; jentry<evts.size();jentry++) {
     
     double p[6];
	 double f[6];
	 double fl[6];
	 double fu[6];
     // total invariant mass
	 // without correction
     mass = evts.at(jentry).InvMass();
     if (mass>beamlow-0.001 && mass<beamup+0.001) datarawo->Fill();
	 for (int i=0;i<6;i++) p[i] = evts.at(jentry).GetP(i);

     short flag = 0x0;
	 for (int pid=0;pid<6;pid++){
       for (int i=0;i<Npart;i++){
         if (p[pid]>=pcut[i]&&p[pid]<pcut[i+1]){
           f[pid]  = facmap[i];
           fl[pid] = facmap[i]-facemap[i];
           fu[pid] = facmap[i]+facemap[i];
		   flag = (flag<<1) +1; //if factor find, the flag should be 00111111
           break;
         }
       }
	 }
	 if (flag!=0x3f){
	   std::cout<<"Waring: not good event, factor map is "<<flag<<std::endl;
	   continue;
	 }
	 // for average correction factor
     mass = evts.at(jentry).InvMass(f);
     if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
	 // factor at low edge
     mass = evts.at(jentry).InvMass(fl);
     if (mass>beamlow-0.001 && mass<beamup+0.001) datarawl->Fill();
	 // factor at up edge
     mass = evts.at(jentry).InvMass(fu);
     if (mass>beamlow-0.001 && mass<beamup+0.001) datarawu->Fill();
  }
  //dataraw->Write();
  // no correction
  sprintf(tmpchr,"raw_%s",namesfx);
  EEto6PI::FitSpectrum(datarawo,beame,tmpchr);

  // factor at average
  sprintf(tmpchr,"nom_%s",namesfx);
  EEto6PI::FitSpectrum(dataraw,beame,tmpchr);
 
  // factor at low edge
  sprintf(tmpchr,"low_%s",namesfx);
  EEto6PI::FitSpectrum(datarawl,beame,tmpchr);

  // factor at up edge
  sprintf(tmpchr,"upv_%s",namesfx);
  EEto6PI::FitSpectrum(datarawu,beame,tmpchr);

  //~~~~~~~~~~ part end~~~~~~~~
  return;
}

void EEto6PI::FitSpectrum(TTree *&dataraw, double beame, char* namesfx)
{
   int nBins=100;
   int Npar;
   double peakvalue = beame;
   double beamlow=beame-0.06;
   double beamup=beame+0.04;
   // try to use roofit
   RooRealVar x("x","energy",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue-0.0005,beamlow,beamup);
   RooRealVar sigma("sigma","width of gaussian",0.003,0.0001,0.02);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   //RooRealVar co1("co1","coefficient #1",0,-100.,100.);
   //RooRealVar co4("co4","coefficient #4",0);
   //RooChebychev bkg("bkg","background",x,RooArgList(co1));
   RooRealVar signal("signal"," ",200,0,1000000);//event number
   RooRealVar background("background"," ",20,0,100000);
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
   sprintf(tmpchr,"data_6pi_%s",namesfx);
   //data_6pi = new RooDataHist(tmpchr,"data_6pi",x,h);
   xframe = x.frame(Title("fit 6 pi"));
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(gaus,ground),RooArgList(signal,background));
   Npar = 6;
   //sigma.setVal(0.035);
   //signal.setVal(1200);
   //background.setVal(200);
   //co1.setVal(0);
   sum->fitTo(*dataset,Range(beamlow,beamup));
   dataset->plotOn(xframe);
   sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(ground),LineStyle(2),LineColor(3));
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
   outf<<beame<<"\t"<<namesfx<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   
   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   delete dataset;
   delete sum;
   return;
}


double EEto6PI::GetEnergy(int runNo)
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
