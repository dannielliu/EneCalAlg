// use D0 D0bar to reconstruct cms energy
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
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "EventClass.h"
//#include <iostream>
extern std::string outputdir;
using namespace RooFit;

namespace KKPIPI{
  double FitSpe(std::vector<D0KPI> &evts, double beame,  const char* namesfx);
  void FitSpe(std::vector<APair> &evts, double beame, const char* namesfx, double f1, double f2);
  void FitSpectrum(TTree *&dataraw,double beame, const char* namesfx, double &peak, double &peakerror);
  void FitSpectrum(TTree *&dataraw,double beame, const char* namesfx);
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
   double beamene = GetEnergy(run);
   //beamene = 4.26;
   std::cout<<"current beam energy is "<<beamene<<", run id "<<run<<std::endl;
   if (beamene < 0.1){
     std::cout<<"can not get a suitable beam energy!"<<std::endl;
	 return;
   }

   std::cout<<"Toral entry is "<<nentries<<std::endl;
   int nBins=100;
   double peakvalue=beamene;// mbeam
   double beamlow=beamene-0.2;
   double beamup=beamene+0.2;
   double mpi=0.13957;
   double mka = 0.493677;
   double D0peak = 1.86484;
   double D0low  = D0peak - 0.1;
   double D0up   = D0peak + 0.1;
    
  char fname[1000];
  sprintf(fname,"%s/plot_DDbar.root",outputdir.c_str());
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
   
   D0KPI evtD;
   std::vector<D0KPI> evtsD;
   D0KPI evtDbar;
   std::vector<D0KPI> evtsDbar;
   APair evt;
   std::vector<APair> evts;

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

	  HepLorentzVector kap,kam,pip,pim;
	  kap.setVectM(Hep3Vector(kappx,kappy,kappz),mka);
	  kam.setVectM(Hep3Vector(kampx,kampy,kampz),mka);
	  pip.setVectM(Hep3Vector(pippx,pippy,pippz),mpi);
	  pim.setVectM(Hep3Vector(pimpx,pimpy,pimpz),mpi);

	  evtD.set(kam,pip);
	  evtDbar.set(kap,pim);
	  evt.set(evtD.Get4P(),evtDbar.Get4P());

	  ppi1 = pip.rho();
	  ppi2 = pim.rho();
	  pk1 = kap.rho();
	  pk2 = kam.rho();
	  chisq = kkm4;
	  //if (pk1>1.3) continue;
	  //if (pk2>1.3) continue;
          mass = evtD.m();
	  char dd = 0;
	  if (mass>D0low && mass<D0up){
		evtsD.push_back(evtD);
		dd |= 0x01;
	  }
          mass = evtDbar.m();
	  if (mass>D0low && mass<D0up){
		evtsDbar.push_back(evtDbar);
		dd |= 0x02;
	  }
          mass = evt.m();
	  if (mass>beamlow && mass<beamup && dd == 0x03){
	    	vars->Fill();
		evts.push_back(evt);
	  }
   }//select end
   vars->Write();
   double f1, f2;
   sprintf(name,"D_%.3f",beamene);
   f1 = KKPIPI::FitSpe(evtsD   ,D0peak,name);
   sprintf(name,"Dbar_%.3f",beamene);
   f2 = KKPIPI::FitSpe(evtsDbar,D0peak,name);
   
   sprintf(name,"CM_%.3f",beamene);
   KKPIPI::FitSpe(evts ,beamene, name,f1,f2);
   return;

}

double KKPIPI::FitSpe(std::vector<D0KPI> &evts, double beame, const char *namesfx)
{
  double peak = beame;
  double beamlow=peak-0.1;
  double beamup=peak+0.1;
  // for factor fit
 
  TTree *dataraw = new TTree("dataraw","dataraw");
  double mass;
  dataraw->Branch("x",&mass,"x/D");
 
  // try to correct the spectrum
  // iniialize the fit function
 
  char tmpchr[100];

  int np = 20;
  double factors[np];
  double factorserr[np];
  for (int i=0; i<np; i++){
  	factors[i] = (1.005-0.995)/np*i+0.995;
	factorserr[i] = 0;
  }
  double peaks[np];
  double deltapeaks[np];
  double peakerrors[np];

  // dm Vs fk
  for (int i=0; i<np; i++){
    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
	  D0KPI evt = evts.at(evtid);
	  evt.setCorrectionFactors(factors[i]);
	  mass = evt.m();
      	  if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
	}
    	sprintf(tmpchr,"factor_%.4f",factors[i]);
    	KKPIPI::FitSpectrum(dataraw,beame,tmpchr,peaks[i],peakerrors[i]);
	deltapeaks[i] = peaks[i] - peak;
	if (peakerrors[i]<2e-5) peakerrors[i] = 1e-3;
  }

   TCanvas *c1 = new TCanvas("c1_1","c1");
   TGraphErrors *graph1 = new TGraphErrors(np,factors,deltapeaks,factorserr,peakerrors);
   graph1->SetTitle("delta peak");
   graph1->SetMarkerStyle(5);
   graph1->Draw("AP");
   gStyle->SetOptFit(1111);
   TF1 *facfit = new TF1("facfit","[0]*(x-[1])");
   facfit->SetParameters(0.5,1);
   facfit->SetParNames("slope","factor");
   graph1->Fit(facfit,"","",factors[0],factors[np-1]);
   double a = facfit->GetParameter(0);
   double ae= facfit->GetParError(0);
   double b = facfit->GetParameter(1);
   double be= facfit->GetParError(1);
   double factor = b;
   double factorerr = sqrt(TMath::Power(peakerrors[(np-1)/2]/a,2)+TMath::Power(be,2));

   graph1->SetName("PsVF");
   graph1->Draw();

  TPaveText *pt = new TPaveText(0.12,0.7,0.50,0.80,"BRNDC");
  sprintf(tmpchr,"factor = %1.6f #pm %1.6f",factor,factorerr);
  pt->AddText(tmpchr);
  sprintf(tmpchr,"interupt = %1.6f #pm %1.6f",b,be);
  pt->AddText(tmpchr);
  pt->Draw();
  sprintf(tmpchr,"PsVsF_%s",namesfx);
  c1->Write(tmpchr);

    dataraw->Reset();
    for (int evtid=0; evtid<evts.size();evtid++){
	  D0KPI evt = evts.at(evtid);
	  evt.setCorrectionFactors(factor);
	  mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"factor_%.6f",factor);
    double peakt,errt;
    KKPIPI::FitSpectrum(dataraw,beame,tmpchr,peakt,errt);

  //~~~~~~~~~~ part end~~~~~~~~
  delete c1;
  return factor;
}

void KKPIPI::FitSpe(std::vector<APair> &evts, double beame, const char *namesfx, double f1, double f2)
{
  double peak = beame;
  double beamlow=peak-0.15;
  double beamup=peak+0.15;
  // for factor fit
 
  TTree *dataraw = new TTree("dataraw","dataraw");
  double mass;
  dataraw->Branch("x",&mass,"x/D");
 
  // try to correct the spectrum
  // iniialize the fit function
 
  char tmpchr[100];

   for (int evtid=0; evtid<evts.size();evtid++){
	  APair evt = evts.at(evtid);
	  evt.setCorrectionFactors(f1,f2);
	  mass = evt.m();
      if (mass>beamlow-0.001 && mass<beamup+0.001) dataraw->Fill();
    }
    sprintf(tmpchr,"CMScor");
    KKPIPI::FitSpectrum(dataraw,beame,tmpchr);

  //~~~~~~~~~~ part end~~~~~~~~
}


void KKPIPI::FitSpectrum(TTree *&dataraw, double beame, const char* namesfx, double &peak, double &peakerror)
{
   int nBins=100;
   bool largesample = false;
   //if (dataraw->GetEntries()>10000) largesample = true;
   int Npar;
   double peakvalue = beame;
   double beamlow=beame-0.05;
   double beamup=beame+0.05;

   // try to use roofit
   RooRealVar x("x","energy",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue,beamlow,beamup);
   RooRealVar sigma("sigma","width of gaussian",0.0056,0.005,0.007);
   RooRealVar sigma2("sigma2","width of gaussian",0.0075,0.007,0.010);
   RooGaussian gaus("gaus","gauss(x,m,s)",x,mean,sigma);
   RooGaussian gaus2("gaus2","gauss(x,m,s)",x,mean,sigma2);
   //RooRealVar co1("co1","coefficient #1",0,-100.,100.);
   //RooRealVar co4("co4","coefficient #4",0);
   //RooChebychev bkg("bkg","background",x,RooArgList(co1));
   RooRealVar signal("signal"," ",120,0,10000000);//event number
   RooRealVar signal2("signal2"," ",120,0,10000000);//event number
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
   sprintf(tmpchr,"data_kpi_%s",namesfx);
   //data_6pi = new RooDataHist(tmpchr,"data_6pi",x,h);
   xframe = x.frame(Title("fit kpi"));
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   sum = new RooAddPdf("sum","sum",RooArgList(gaus,gaus2,ground),RooArgList(signal,signal2,background));
   Npar = 8;
   double factor = atof(&namesfx[7]);
   mean.setVal(peakvalue+1.4*(factor-1.0));
   //sigma.setVal(0.035);
   //signal.setVal(1200);
   //background.setVal(200);
   //co1.setVal(0);
   sum->fitTo(*dataset,Range(beamlow,beamup));
   //sum->fitTo(*dataset,Range(beame-0.05,beame+0.03));
   //sum->fitTo(*dataset,Range(4.22,4.28));
   dataset->plotOn(xframe);
   //sum->plotOn(xframe,Components(cbshape1),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
   sum->plotOn(xframe,Components(gaus2),LineStyle(2),LineColor(2));
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

  peak = mean.getVal();
  peakerror = mean.getError();

   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   delete dataset;
   delete sum;
   return;
}

void KKPIPI::FitSpectrum(TTree *&dataraw, double beame, const char* namesfx)
{
   int nBins=100;
   bool largesample = false;
   if (dataraw->GetEntries()>10000) largesample = true;
   int Npar;
   double peakvalue = beame;
   double beamlow=beame-0.15;
   double beamup=beame+0.15;
   // try to use roofit
   RooRealVar x("x","energy",peakvalue,beamlow,beamup,"GeV");
   RooRealVar mean("mean","mean of gaussian",peakvalue,beamlow,beamup);
   RooRealVar sigma("sigma","width of gaussian",0.014,0.005,0.20);
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
     
   RooRealVar alpha1("alpha1","#alpha",1.32,-5,5);
   RooRealVar nnn1("n1","n",5,1,200);
   RooCBShape cbshape1("cbshape1","crystal ball",x,mean,sigma,alpha1,nnn1);
   
   RooAddPdf *sum;
   RooDataSet *dataset;
   RooPlot *xframe;
   //RooDataHist *data_6pi;
 
   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); // set out put message level of roofit
   TCanvas *c1=new TCanvas("","",800,600);
	 
   char tmpchr[100];
   sprintf(tmpchr,"data_DDbar_%s",namesfx);
   //data_6pi = new RooDataHist(tmpchr,"data_6pi",x,h);
   xframe = x.frame(Title("fit D^{0} #bar{D}^{0}"));
   dataset = new RooDataSet(tmpchr,"data",RooArgSet(x),Import(*dataraw));
   if (!largesample) {
     //sum = new RooAddPdf("sum","sum",RooArgList(gaus,ground),RooArgList(signal,background));
     sum = new RooAddPdf("sum","sum",RooArgList(cbshape1,ground),RooArgList(signal,background));
     Npar = 8;
   }
   else {
     sum = new RooAddPdf("sum","sum",RooArgList(cbshape1,gaus2,ground),RooArgList(signal,signal2,background));
	 Npar=10;
   }
   //sigma.setVal(0.035);
   //signal.setVal(1200);
   //background.setVal(200);
   //co1.setVal(0);
   sum->fitTo(*dataset,Range(beamlow,beamup));
   //sum->fitTo(*dataset,Range(beame-0.05,beame+0.03));
   //sum->fitTo(*dataset,Range(4.22,4.28));
   dataset->plotOn(xframe);
   sum->plotOn(xframe,Components(cbshape1),LineStyle(2),LineColor(2));
   //sum->plotOn(xframe,Components(gaus),LineStyle(2),LineColor(2));
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
   outf<<beame<<"\t"<<mean.getVal()<<"\t"<<mean.getError()<<std::endl;
   
   //c1->Print("fit6pi.eps");
   //delete data_6pi;
   delete xframe;
   delete dataset;
   delete sum;
   return;
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
